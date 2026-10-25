"""Single bubble in a primitive three-loop vertex: exact original-forest binding."""
from itertools import combinations
import sympy as S
import primitive_triangle_factor_v3 as P


def host_polynomials(routes):
    """Cauchy--Binet coefficients, independently of symbolic Schur determinants."""
    matrix=S.Matrix.vstack(*routes);n=len(routes);U=[];F=[]
    for subset in combinations(range(n),3):
        d=matrix.extract(subset,range(3)).det();co=S.expand(d*d)
        if co:U.append([sum(1<<i for i in subset),int(co)])
    for subset in combinations(range(n),4):
        a=matrix.extract(subset,[0,1,2,3]).det();b=matrix.extract(subset,[0,1,2,4]).det()
        co=S.expand(a*a+b*b-a*b)
        assert co.is_Integer and co>=0
        if co:F.append([sum(1<<i for i in subset),int(co)])
    assert U and F and all(co==1 for _,co in U)
    return U,F


def polynomial(rows,alpha):
    return S.expand(sum(co*S.prod(alpha[i] for i in range(len(alpha)) if mask&(1<<i)) for mask,co in rows))


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,4,3,0,12)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert len(gd['proper_subgraphs'])==1
    part=gd['proper_subgraphs'][0];assert (part['loops'],part['external_legs'],part['omega'],len(part['edge_ids']))==(1,2,2,2)
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']}=={(),(0,)}
    assert len(gd['terms'])==len(pd['states'])==2 and {tuple(s['word']) for s in pd['states']}=={(),(0,)}
    assert raw['external_gram']==[[1.,-.5],[-.5,1.]] and raw['positive_invariants']==[1.,1.,1.]
    alpha=S.symbols('a:12');z=S.symbols('z:1');gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    selected=part['edge_ids'];outside=[i for i in range(12) if i not in selected];lam=sum(alpha[e] for e in selected)
    choices,witnesses=P.basis_all(gd,0,z);records=[];states={}
    for term,witness in zip(gd['terms'],witnesses):
        word=tuple(term['subgraphs_inner_first']);state=next(s for s in pd['states'] if tuple(s['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==state['variables']==list(word) and term['jet_orders']==([2] if word else [])
        route=witness['routing'];M=sum((a*r.T*r for a,r in zip(alpha,route)),S.zeros(6))
        host=sum((alpha[e]*route[e][:,1:].T*route[e][:,1:] for e in outside),S.zeros(5))
        vectors=[route[e][:,1:]/route[e][0] for e in selected];delta=vectors[0]-vectors[1]
        variance=alpha[selected[0]]*alpha[selected[1]]*delta.T*delta
        schur=(lam*M[1:,1:]-M[1:,:1]*M[:1,1:]).applyfunc(S.expand)
        assert M[0,0]==lam and (schur-lam*host-variance).applyfunc(S.expand)==S.zeros(5)
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(6))
        U=original[:4,:4].det(method='domain-ge')
        F=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        reduced_u=S.cancel(schur[:3,:3].det(method='domain-ge')/lam**2)
        reduced_f=S.cancel(sum(gram[i,j]*schur.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2))/lam**3)
        counts=[]
        for name,original_value,reduced in (('U',U,reduced_u),('F',F,reduced_f)):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(original_value-positive)==0 and S.expand(reduced-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states[word]=(host,delta,route)
        records.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
            'basis_edges':list(choices),'jacobian_determinant':int(witness['det']),
            'loop_change':[[str(v) for v in witness['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in witness['offset'].row(i)] for i in range(4)],
            'routing':[[str(v) for v in r] for r in route],'U_F_monomials':counts,
            'three_exact_U_F_constructions_equal':True})
    host,delta,route=states[()];projected,scaled,_=states[(0,)]
    assert (host-projected).applyfunc(S.expand)==S.zeros(5)
    assert (scaled-z[0]*delta).applyfunc(S.expand)==S.zeros(1,5)
    duplicate=[e for e in outside if route[e][:,1:] in (delta,-delta)]
    assert len(duplicate)==2
    kept=[e for e in outside if e not in duplicate]+[duplicate[0]]
    assert len(kept)==9
    rows=[route[e][:,1:] for e in kept];U,F=host_polynomials(rows);beta=S.symbols('b:9')
    reduced=sum((b*r.T*r for b,r in zip(beta,rows)),S.zeros(5))
    assert S.expand(reduced[:3,:3].det()-polynomial(U,beta))==0
    ff=sum(gram[i,j]*reduced.extract([0,1,2,3+i],[0,1,2,3+j]).det() for i in range(2) for j in range(2))
    assert S.expand(ff-polynomial(F,beta))==0
    cap=max((r[:,3:]*gram*r[:,3:].T)[0] for r in rows);assert cap>0
    return {'id':raw['id'],'states':records,'bubble_edges':selected,'paired_host_propagators':duplicate,
        'host_edges':kept,'weighted_host_edge':8,'host_routes':[[int(v) for v in r] for r in rows],
        'host_U':U,'host_F':F,'external_energy_trial_bound':int(cap),
        'host_unchanged_under_bubble_Taylor':True,'bubble_boundary_scaled_exactly_once':True,
        'weighted_edge':'B(k^2)/(1+k^2)^2','forest_pairing':'bare minus vector-total-degree-two bubble Taylor',
        'whole_integral_evaluated':False}


def scalar_identity():
    M,d,z,s=S.symbols('M d z s',positive=True)
    C=(s-2)/(1-s)**2;D=1/(1-s);E=1/(1-s)**2;mass=M+d*(1/s-1)
    value=-C*S.log(1+z/M)+D*d*(1/(M+z)-1/M)-E/s*S.log(1+z/mass)
    second=-(M+z-d)**2/((M+z)**2*(d+s*(M+z-d)))
    assert S.cancel(S.diff(value,z)-second)==0 and value.subs(z,0)==0
    t=(M+z-d)/d
    assert S.cancel(second+d*t*t/((M+z)**2*(1+s*t)))==0
    return {'J_prime':'-d*B((M+z-d)/d)/(M+z)^2','J_zero':0,
            'J':'-d*z*integral_0^1 B((M+u*z-d)/d)/(M+u*z)^2 du',
            'whole_integral':'projective_integral J(F/U)/U^3 over nine host Schwinger parameters',
            'full_external_momentum_dependence_retained':True,'integrand_nonpositive':True}
