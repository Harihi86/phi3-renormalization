"""Bind a shared overlapping two-loop two-point core to every original forest.

No whole-integral accuracy claim is made by these identities.  The two
overlapping triangles are never admitted into the same forest.
"""
import sympy as S
import primitive_triangle_factor_v3 as P

IDS = ('L3E2G7', 'L3E3G38')
CORE = (0, 1, 3, 4, 6)


def verify(raw):
    gid = raw['id']; assert gid in IDS
    gd, pd = raw['gaussian_full_descriptor'], raw['positive_full_descriptor']
    n = 1 if gid == IDS[0] else 2
    assert (gd['dimension'], gd['loops'], gd['external_legs'], gd['overall_degree']) == (6, 3, n+1, 2 if n == 1 else 0)
    assert gd['normalization'] == 'd^6k/pi^3 per loop; no symmetry weight or coupling'
    assert [p['edge_ids'] for p in gd['proper_subgraphs']] == [[0,1,6],[3,4,6],list(CORE)]
    assert [(p['loops'],p['external_legs'],p['omega']) for p in gd['proper_subgraphs']] == [(1,3,0),(1,3,0),(2,2,2)]
    words = {(), (0,), (1,), (2,), (0,2), (1,2)}
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']} == words
    assert {tuple(t['word']) for t in pd['states']} == words
    assert len(gd['terms']) == len(pd['states']) == 6
    gram = S.Matrix(raw['external_gram']).applyfunc(S.Rational)
    expected_gram = S.Matrix([[1]]) if n == 1 else S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    assert gram == expected_gram and all(S.Rational(v)==1 for v in raw['positive_invariants'])
    alpha=S.symbols('a:'+str(gd['edges'])); z=S.symbols('z:3')
    a=alpha[0]+alpha[1]; b=alpha[3]+alpha[4]; c=alpha[6]
    U=a*b+c*(a+b)
    F=(alpha[0]*alpha[1]*b+alpha[3]*alpha[4]*a
       +c*(alpha[0]+alpha[3])*(alpha[1]+alpha[4]))
    polynomials={None:(U,F),0:(b*(a+c),(a+c)*alpha[3]*alpha[4]),
                 1:(a*(b+c),(b+c)*alpha[0]*alpha[1])}
    host_edges=[e for e in range(gd['edges']) if e not in CORE]
    fixed_host=None; records=[]
    for term in gd['terms']:
        word=tuple(term['subgraphs_inner_first']); state=next(s for s in pd['states'] if tuple(s['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==state['variables']==([2] if 2 in word else [])
        assert term['jet_orders']==([2] if 2 in word else [])
        old=P.scaled_routes(term,z)
        # q = original third loop minus p (E2), or minus p1+p2 (E3).
        # This is the SAME unimodular affine change for all six forest states.
        route=[S.Matrix([[r[j] if j<3 else S.expand(r[j]+r[2]) for j in range(3+n)]]) for r in old]
        for e in CORE:assert all(v==0 for v in route[e][3:])
        qscale=z[2] if 2 in word else S.Integer(1)
        tri=0 if 0 in word else 1 if 1 in word else None
        expected={None:[[1,1,qscale],[-1,-1,0],[-1,0,-qscale],[1,0,0],[0,1,0]],
                  0:[[0,1,0],[0,-1,0],[-1,0,-qscale],[1,0,0],[0,1,0]],
                  1:[[-1,0,qscale],[1,0,0],[0,1,0],[0,-1,0],[0,1,0]]}[tri]
        assert [list(route[e][:3]) for e in CORE] == expected
        child=sum((alpha[e]*route[e][:,:3].T*route[e][:,:3] for e in CORE),S.zeros(3))
        uc,fc=polynomials[tri]
        assert S.expand(child[:2,:2].det()-uc)==0
        assert S.expand(child.det()-qscale*qscale*fc)==0
        host=sum((alpha[e]*route[e][:,2:].T*route[e][:,2:] for e in host_edges),S.zeros(1+n))
        if fixed_host is None:fixed_host=host
        assert host==fixed_host, 'Subgraph Taylor changed a host propagator'
        reduced=host.copy();reduced[0,0]+=qscale*qscale*fc/uc
        before=sum((v*r.T*r for v,r in zip(alpha,old)),S.zeros(3+n))
        counts=[]
        for name,original,new in (
            ('U',before[:3,:3].det(method='domain-ge'),uc*reduced[0,0]),
            ('F',sum(gram[i,j]*before.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(n) for j in range(n)),
             uc*sum(gram[i,j]*reduced.extract([0,1+i],[0,1+j]).det() for i in range(n) for j in range(n)))):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(original-positive)==0 and S.cancel(new-positive)==0,(gid,word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        records.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
                        'routing':[[str(v) for v in r] for r in route],
                        'core_U':str(uc),'core_F_per_q2':str(fc),'U_F_monomials':counts,
                        'three_exact_U_F_constructions_equal':True})
    return {'id':gid,'states':records,'host_edges':host_edges,'shared_core_edges':list(CORE),
            'all_six_original_forests_bound':True,'common_affine_jacobian':1,
            'boundary':'q; retained host loop, not the physical external momentum',
            'core_Taylor':'total vector degree 2; after rotation, constant and q squared',
            'whole_integral_evaluated':False}


def identities():
    a,b,c,x,y,t,u,rho=S.symbols('a b c x y t u rho',positive=True)
    aa=[a*x,a*(1-x),b*y,b*(1-y),c]
    U=a*b+c*(a+b);X=x*(1-x);Y=y*(1-y)
    F=aa[0]*aa[1]*b+aa[2]*aa[3]*a+c*(aa[0]+aa[2])*(aa[1]+aa[3])
    Q=a*X+b*Y+c*a*b/U*(x-y)**2
    assert S.cancel(F/U-Q)==0
    # Cheng-Wu section c=1, alpha=rho*(a*x,a*(1-x),b*y,b*(1-y),1).
    transformation=S.Matrix([rho*v for v in aa]).subs(c,1)
    jac=transformation.jacobian([rho,a,b,x,y]).det()
    assert S.factor(jac)==-rho**4*a*b
    h=(1+t)*S.log(1+t)-t
    assert h.subs(t,0)==0 and S.diff(h,t).subs(t,0)==0
    assert S.diff(h,t,2)==1/(1+t)
    assert S.simplify(t*t*S.integrate((1-u)/(1+u*t),(u,0,1))-h)==0
    # Two differentiations justify the common scale integral, with both
    # constants fixed at t=0: integral exp(-rho*M)(exp(-rho*t*Q)-1+rho*t*Q)/rho^2.
    M,q=S.symbols('M q',positive=True)
    integrated=M*((1+t*q/M)*S.log(1+t*q/M)-t*q/M)
    assert S.simplify(S.diff(integrated,t,2)-q*q/(M+t*q))==0
    assert integrated.subs(t,0)==0 and S.diff(integrated,t).subs(t,0)==0
    return {'Cheng_Wu_section':'central edge c=1; four coordinates a,b,x,y',
            'Jacobian':'rho^4 a b',
            'common_scale_integral':'M h(t Q/M), M=1+a+b; h(v)=(1+v)log(1+v)-v',
            'Q_bare':'a*x*(1-x)+b*y*(1-y)+a*b/(a*b+a+b)*(x-y)^2',
            'Q_triangle0':'b*y*(1-y)','Q_triangle1':'a*x*(1-x)',
            'U_bare':'a*b+a+b','U_triangle0':'b*(a+1)','U_triangle1':'a*(b+1)',
            'core_integrand':'a*b*M*[h(t*Q/U_mass)/U^3 - h(t*Q0/U_mass)/U0^3 - h(t*Q1/U_mass)/U1^3], U_mass=M',
            'independent_remainder':'h(v)=v^2 integral_0^1 (1-u)/(1+u*v) du',
            'no_overlapping_double_counterterm':True,
            'endpoint_integrability_proved':False,'whole_integral_evaluated':False}


def core_density(mp,t,a,b,x,y):
    t,a,b,x,y=map(mp.mpf,(t,a,b,x,y));U=a*b+a+b;M=1+a+b
    X=x*(1-x);Y=y*(1-y);Q=a*X+b*Y+a*b/U*(x-y)**2
    h=lambda v:(1+v)*mp.log1p(v)-v
    return a*b*M*(h(t*Q/M)/U**3-h(t*b*Y/M)/(b*(a+1))**3-h(t*a*X/M)/(a*(b+1))**3)


def core_remainder_density(mp,t,a,b,x,y,u):
    t,a,b,x,y,u=map(mp.mpf,(t,a,b,x,y,u));U=a*b+a+b;M=1+a+b
    X=x*(1-x);Y=y*(1-y);Q=a*X+b*Y+a*b/U*(x-y)**2;Q0=b*Y;Q1=a*X
    return t*t*(1-u)*a*b*(Q*Q/(U**3*(M+u*t*Q))
            -Q0*Q0/((b*(a+1))**3*(M+u*t*Q0))
            -Q1*Q1/((a*(b+1))**3*(M+u*t*Q1)))
