"""Exact L3E3G8 reduction: a bubble inside a two-loop vertex.

The bubble and the second child loop are integrated analytically. All four
original forest states and the final physical-external-momentum subtraction
remain present. This is an integrand identity, not a numerical certificate.
"""
import sympy as S
import primitive_triangle_factor_v3 as P
from disjoint_triangle_bubble_factor import convert


def verify(raw):
    assert raw['id']=='L3E3G8'
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,3,3,0,9)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    parts=gd['proper_subgraphs'];assert len(parts)==2
    assert (parts[0]['edge_ids'],parts[0]['loops'],parts[0]['external_legs'],parts[0]['omega'])==([6,7],1,2,2)
    assert (parts[1]['edge_ids'],parts[1]['loops'],parts[1]['external_legs'],parts[1]['omega'])==([0,2,3,5,6,7],2,3,0)
    assert set(parts[0]['edge_ids'])<set(parts[1]['edge_ids'])
    assert len(gd['terms'])==len(pd['states'])==4
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    assert list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    alpha=S.symbols('a:9');z=S.symbols('z:2');states=[];seen=set();fixed_host=None
    for term in gd['terms']:
        word=tuple(term['subgraphs_inner_first']);seen.add(word)
        state=next(r for r in pd['states'] if tuple(r['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==state['variables']==([0] if 0 in word else [])
        assert term['jet_orders']==([2] if 0 in word else [])
        original=P.scaled_routes(term,z);w=convert(original,(7,3,8));assert w is not None
        routes=w['routes'];b=z[0] if 0 in word else S.Integer(1);v=0 if 1 in word else 1
        expected=[[0,-1,0,0,0],[0,0,-1,0,-1],[0,1,v,0,v],
                  [0,1,0,0,0],[0,0,1,-1,0],[0,-1,-v,v,0],
                  [-1,b,0,0,0],[1,0,0,0,0],[0,0,1,0,0]]
        assert routes==[S.Matrix([row]) for row in expected]
        matrix=sum((a*r.T*r for a,r in zip(alpha,routes)),S.zeros(5))
        host=sum((alpha[e]*routes[e][:,1:].T*routes[e][:,1:] for e in (1,4,8)),S.zeros(4))
        if fixed_host is None:fixed_host=host
        assert host==fixed_host,'Child Taylor operation changed a host propagator'
        child=sum((alpha[e]*routes[e][:,1:].T*routes[e][:,1:] for e in (0,2,3,5)),S.zeros(4))
        local=alpha[6]+alpha[7];variance=S.zeros(4);variance[0,0]=alpha[6]*alpha[7]*b*b/local
        schur=matrix[1:,1:]-matrix[1:,:1]*matrix[:1,1:]/local
        assert matrix[0,0]==local and (schur-host-child-variance).applyfunc(S.cancel)==S.zeros(4)
        before=sum((a*r.T*r for a,r in zip(alpha,original)),S.zeros(5))
        counts=[]
        for name,old,new in [
            ('U',before[:3,:3].det(method='domain-ge'),local*schur[:2,:2].det(method='domain-ge')),
            ('F',sum(gram[i,j]*before.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2)),
             local*sum(gram[i,j]*schur.extract([0,1,2+i],[0,1,2+j]).det(method='domain-ge') for i in range(2) for j in range(2)))]:
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(old-positive)==0 and S.cancel(new-positive)==0,(word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
            'basis_original_edges':[7,3,8],'jacobian_determinant':int(w['det']),
            'routing':[[str(t) for t in row] for row in routes],
            'three_exact_U_F_constructions_equal':True,'U_F_monomials':counts})
    assert seen=={(),(0,),(1,),(0,1)}
    x,y,q2,qp1,qp2=S.symbols('x y t qp1 qp2')
    x1=(1-x)*y;x2=(1-x)*(1-y)
    variance=x*x1*(q2+2*qp2+1)+x*x2*(q2-2*qp1+1)+x1*x2
    v0=x*(1-x)*q2;dv=2*x*(x1*qp2-x2*qp1)+x*(1-x)+x1*x2
    assert S.expand(variance-v0-dv)==0
    return {'graph':raw['id'],'states':states,'all_four_nested_forests_bound':True,
        'bubble_boundary':'k, the retained loop inside the two-loop child vertex',
        'host_edges':[1,4,8],'outer_host':'Delta(q) Delta(q-p1) Delta(q+p2)',
        'child_local_T0':'Only the child boundaries are zeroed; the host is unchanged',
        'overall_T0':'p1=p2=0 at fixed k,q; it does not set q=0',
        'simplex_measure':'(1-x) dx dy','whole_graph_certificate':False}


def identities():
    a,t,x,V,u,s=S.symbols('a t x V u s',positive=True)
    C=(a-2)/(1-a)**2;D=1/(1-a);E=1/(1-a)**2
    lhs=t*t/((1+t)**2*(1+a*t))
    assert S.cancel(lhs-C/(1+t)-D/(1+t)**2-E/(1+a*t))==0
    h=(1+t)*S.log(1+t)-t
    integral=S.integrate((1-u)/(1+u*t),(u,0,1))
    assert S.simplify(h-t*t*integral)==0
    # Integrate x(1-x)-a over the exact allowed x interval, then a=(1-s^2)/4.
    aa=(1-s*s)/4
    weight=S.integrate(x*(1-x)-aa,(x,(1-s)/2,(1+s)/2))
    assert S.expand(weight-s**3/6)==0
    assert S.integrate(s**4/12,(s,0,1))==S.Rational(1,60)
    mass=1-x+x/a
    combined=-C*S.log(1+V)-D*x*V/(1+V)-E/a*S.log(1+V/mass)
    derivative=-(1+V-x)**2/((1+V)**2*(x+a*(1-x+V)))
    assert S.cancel(S.diff(combined,V)-derivative)==0
    assert combined.subs(V,0)==0
    # The spectral B identity turns the entire derivative into one B evaluation.
    b=1-x+V;z=b/x
    assert S.cancel((x/(1+V)**2)*(z*z)/(1+a*z)-b*b/((1+V)**2*(x+a*b)))==0
    # Check the unit-Jacobian primitive substitution, including the lower end.
    assert S.cancel((x/(1+V)**2)/(S.diff(z,V))-1/(1+z)**2)==0
    return {'Stieltjes_weight':'(1-4a)^(3/2)/6, 0<=a<=1/4',
        'bubble':'B(t)=t^2 integral w(a)/(1+a*t) da',
        'child_kernel':'K(x,V)=-integral_{(1-x)/x}^{(1-x+V)/x} B(t)/(1+t)^2 dt',
        'derivative':'dK/dV=-x*B((1-x+V)/x)/(1+V)^2',
        'mass_derivative_and_spectral_coefficients_exact':True,
        'bubble_positive_leading_coefficient':'1/60','whole_graph_certificate':False}


def bubble(mp,t):
    t=mp.mpf(t)
    if not t:return mp.mpf(0)
    return (t+4)**mp.mpf('1.5')/(3*mp.sqrt(t))*mp.asinh(mp.sqrt(t)/2)-4*t/9-mp.mpf(4)/3


def spectral_kernel(mp,x,V):
    x,V=mp.mpf(x),mp.mpf(V)
    def integrand(s):
        a=(1-s*s)/4
        if not a:return s**4/12*(2*mp.log1p(V)-x*V/(1+V)-V/x)
        c=(a-2)/(1-a)**2;d=1/(1-a);e=1/(1-a)**2;m=1-x+x/a
        return s**4/12*(-c*mp.log1p(V)-d*x*V/(1+V)-e/a*mp.log1p(V/m))
    return mp.quad(integrand,[0,mp.mpf('.5'),1])


def primitive_kernel(mp,x,V):
    x,V=mp.mpf(x),mp.mpf(V)
    return -V*x*mp.quad(lambda s:bubble(mp,(1-x+s*V)/x)/(1+s*V)**2,[0,mp.mpf('.5'),1])
