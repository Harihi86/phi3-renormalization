"""L3E2G5: retain the overlap forest while analytically eliminating the bubble."""
import sympy as S
import primitive_triangle_factor_v3 as P


def verify(raw):
    assert raw['id']=='L3E2G5'
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,3,2,2,8)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert [p['edge_ids'] for p in gd['proper_subgraphs']]==[[3,4],[0,1,7],[0,2,3,4,5,6]]
    assert [(p['loops'],p['external_legs'],p['omega']) for p in gd['proper_subgraphs']]==[(1,2,2),(1,3,0),(2,3,0)]
    words={(),(0,),(1,),(2,),(0,1),(0,2)}
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']}=={tuple(t['word']) for t in pd['states']}==words
    assert len(gd['terms'])==len(pd['states'])==6 and raw['external_gram']==[[1.0]] and raw['positive_invariants']==[1.0]
    alpha=S.symbols('a:8');z=S.symbols('z:3');records=[]
    for term in gd['terms']:
        word=tuple(term['subgraphs_inner_first']);state=next(s for s in pd['states'] if tuple(s['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==state['variables']==([0] if 0 in word else [])
        assert term['jet_orders']==([2] if 0 in word else [])
        old=P.scaled_routes(term,z);route=[S.Matrix([[r[0],r[1],-r[2],r[3]]]) for r in old]
        b=z[0] if 0 in word else S.Integer(1)
        expected=[[0,1,-1,0],[0,0,1,-1],[0,-1,0,1],[-1,-b,0,0],[1,0,0,0],[0,1,0,0],[0,1,0,0],[0,0,-1,0]]
        if 1 in word:expected[0]=[0,0,-1,0];expected[1]=[0,0,1,0]
        if 2 in word:expected[0]=[0,1,0,0];expected[2]=[0,-1,0,0]
        assert route==[S.Matrix([r]) for r in expected]
        before=sum((a*r.T*r for a,r in zip(alpha,old)),S.zeros(4))
        matrix=sum((a*r.T*r for a,r in zip(alpha,route)),S.zeros(4))
        bubble=alpha[3]+alpha[4]
        rest=sum((alpha[e]*route[e][:,1:].T*route[e][:,1:] for e in (0,1,2,5,6,7)),S.zeros(3))
        variance=S.zeros(3);variance[0,0]=alpha[3]*alpha[4]*b*b/bubble
        schur=matrix[1:,1:]-matrix[1:,:1]*matrix[:1,1:]/bubble
        assert (schur-rest-variance).applyfunc(S.cancel)==S.zeros(3)
        counts=[]
        for name,left,right in [('U',before[:3,:3].det(method='domain-ge'),bubble*schur[:2,:2].det()),
                                ('F',before.det(method='domain-ge'),bubble*schur.det())]:
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(left-positive)==0 and S.cancel(right-positive)==0,(word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        records.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
                        'routing':[[str(v) for v in row] for row in route],'U_F_monomials':counts,
                        'three_exact_U_F_constructions_equal':True})
    return {'id':raw['id'],'states':records,'common_affine_jacobian':-1,'bubble_boundary':'k, retained loop of the outer kite',
            'weighted_edge':'W(k^2)=B(k^2)/(1+k^2)^2',
            'paired_forest_groups':[{'words':[[],[0]],'sign':1},{'words':[[1],[0,1]],'sign':-1},{'words':[[2],[0,2]],'sign':-1}],
            'both_overlapping_vertex_subtractions_retained':True,'whole_integral_evaluated':False}


def identities():
    sigma,t,M,d,z=S.symbols('sigma t M d z',positive=True)
    c=(sigma-2)/(1-sigma)**2;dd=1/(1-sigma);e=1/(1-sigma)**2
    assert S.cancel(t*t/((1+t)**2*(1+sigma*t))-c/(1+t)-dd/(1+t)**2-e/(1+sigma*t))==0
    h=lambda v:(1+v)*S.log(1+v)-v
    mass=M+d*(1/sigma-1)
    combined=c*M*h(z/M)+dd*d*(z/M-S.log(1+z/M))+e/sigma*mass*h(z/mass)
    second=(M+z-d)**2/((M+z)**2*(d+sigma*(M+z-d)))
    assert S.cancel(S.diff(combined,z,2)-second)==0
    assert combined.subs(z,0)==0 and S.simplify(S.diff(combined,z).subs(z,0))==0
    assert S.simplify(S.diff(M*h(z/M),M)-(S.log(1+z/M)-z/M))==0
    v=(M+z-d)/d
    assert S.cancel(d*v*v/((M+z)**2*(1+sigma*v))-second)==0
    return {'bubble_spectral_weight':'(1-4*sigma)^(3/2)/6, 0<sigma<1/4',
            'weighted_edge':'C/(1+k^2)+D/(1+k^2)^2+E/[sigma*(k^2+1/sigma)] integrated with the weight',
            'double_propagator':'-d/d(mass_squared) of the single weighted propagator; its Schwinger coefficient is d',
            'J_second':'d*B((M+z-d)/d)/(M+z)^2',
            'J':'d*z^2 integral_0^1(1-u)*B((M+u*z-d)/d)/(M+u*z)^2 du',
            'J_initial_conditions':['J(0)=0','J_prime(0)=0'],
            'whole_G5':'integral da db dx dy a*b*[J(Q)/U^3-J(Q0)/U0^3-J(Q1)/U1^3], external p^2=1',
            'parameters':'same kite U,Q,U0,Q0,U1,Q1 as overlap core; M=1+a+b; d=b*(1-y)',
            'endpoint_integrability_proved':False,'whole_integral_evaluated':False}


def spectral_J(mp,M,d,z):
    M,d,z=map(mp.mpf,(M,d,z));h=lambda v:(1+v)*mp.log1p(v)-v
    def f(s):
        sigma=(1-s*s)/4
        if not sigma:return s**4/12*(-2*M*h(z/M)+d*(z/M-mp.log1p(z/M))+z*z/(2*d))
        c=(sigma-2)/(1-sigma)**2;dd=1/(1-sigma);e=1/(1-sigma)**2;mass=M+d*(1/sigma-1)
        return s**4/12*(c*M*h(z/M)+dd*d*(z/M-mp.log1p(z/M))+e/sigma*mass*h(z/mass))
    return mp.quad(f,[0,mp.mpf('.5'),1])


def primitive_J(mp,M,d,z):
    M,d,z=map(mp.mpf,(M,d,z))
    def bubble(t):return t*t/60*mp.hyp2f1(1,1,mp.mpf('3.5'),-t/4)
    return d*z*z*mp.quad(lambda u:(1-u)*bubble((M+u*z-d)/d)/(M+u*z)**2,[0,mp.mpf('.5'),1])
