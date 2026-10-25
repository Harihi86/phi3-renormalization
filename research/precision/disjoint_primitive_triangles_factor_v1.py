"""Exact simultaneous elimination of two or three disjoint triangle children."""
from itertools import combinations,product
import sympy as S
import primitive_triangle_factor_v3 as P


def basis(gd):
    children=[p['edge_ids'] for p in gd['proper_subgraphs']];k=len(children)
    outside=[e for e in range(12) if all(e not in edges for edges in children)]
    routed=[P.scaled_routes(t,[]) for t in gd['terms']]
    for own in product(*children):
        for rest in combinations(outside,4-k):
            choices=(*own,*rest);records=[]
            for rows in routed:
                change=S.Matrix.vstack(*(rows[e][:,:4] for e in choices));det=change.det()
                if det not in (-1,1):break
                offset=S.Matrix.vstack(*(rows[e][:,4:] for e in choices));inverse=change.inv();converted=[]
                for row in rows:
                    loop=row[:,:4]*inverse;external=row[:,4:]-loop*offset
                    converted.append(S.Matrix([[S.expand(v) for v in list(loop)+list(external)]]))
                if any(v for e in outside for v in converted[e][:,:k]):break
                if any(converted[e][i] not in (-1,1) or any(converted[e][j] for j in range(k) if j!=i)
                       for i,edges in enumerate(children) for e in edges):break
                records.append({'change':change,'offset':offset,'det':det,'routing':converted})
            if len(records)==len(routed):return choices,records
    raise AssertionError('No common independent triangle-loop basis')


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert raw['external_gram']==[[1.,-.5],[-.5,1.]] and raw['positive_invariants']==[1.,1.,1.]
    parts=gd['proper_subgraphs'];k=len(parts);assert k in (2,3)
    children=[p['edge_ids'] for p in parts]
    assert all((p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==(1,3,0,3) for p in parts)
    assert all(not set(a)&set(b) for a,b in combinations(children,2))
    wanted={frozenset(i for i in range(k) if mask&(1<<i)) for mask in range(1<<k)}
    assert len(gd['terms'])==len(pd['states'])==1<<k
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==wanted
    for term in gd['terms']:
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if set(s['word'])==set(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==term['jet_orders']==state['variables']==[]
    choices,witnesses=basis(gd);alpha=S.symbols('a:12');L=4-k;n=L+2
    outside=[e for e in range(12) if all(e not in edges for edges in children)]
    lambdas=[sum(alpha[e] for e in edges) for edges in children];combined=S.prod(lambdas)
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]);records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=tuple(term['subgraphs_inner_first']);state=next(s for s in pd['states'] if set(s['word'])==set(word))
        route=w['routing'];matrix=sum((a*r.T*r for a,r in zip(alpha,route)),S.zeros(6))
        host=sum((alpha[e]*route[e][:,k:].T*route[e][:,k:] for e in outside),S.zeros(n))
        assert matrix[:k,:k]==S.diag(*lambdas)
        variances=[];deltas=[]
        for i,edges in enumerate(children):
            shifts=[route[e][:,k:]/route[e][i] for e in edges]
            differences=[shifts[a]-shifts[b] for a,b in combinations(range(3),2)]
            variance=sum((alpha[edges[a]]*alpha[edges[b]]*d.T*d for (a,b),d in zip(combinations(range(3),2),differences)),S.zeros(n))
            variances.append(variance.applyfunc(S.expand));deltas.append(differences)
            if i in word:assert variances[-1]==S.zeros(n)
        reduced=combined*host+sum((S.prod(lambdas[j] for j in range(k) if j!=i)*v for i,v in enumerate(variances)),S.zeros(n))
        direct=combined*matrix[k:,k:]-sum((S.prod(lambdas[j] for j in range(k) if j!=i)*matrix[k:,i:i+1]*matrix[i:i+1,k:] for i in range(k)),S.zeros(n))
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(n)
        reduced_u=S.Poly(reduced[:L,:L].det(method='domain-ge'),*alpha).exquo(S.Poly(combined**(L-1),*alpha)).as_expr()
        fn=sum(gram[i,j]*reduced.extract([*range(L),L+i],[*range(L),L+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        reduced_f=S.Poly(fn,*alpha).exquo(S.Poly(combined**L,*alpha)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,[]))),S.zeros(6))
        ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,value,reduced_value in (('U',ug,reduced_u),('F',fg,reduced_f)):
            positive=P.positive_polynomial(state[name],state,alpha,[])
            assert S.expand(value-positive)==0 and S.expand(reduced_value-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha).terms()))
        states[frozenset(word)]=(host,variances,deltas,route)
        records.append({'word':list(word),'sign':term['sign'],'common_scale_log_coefficient':-term['sign'],
            'jacobian_determinant':int(w['det']),'basis_edges':list(choices),
            'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],
            'routing':[[str(v) for v in row] for row in route],'U_F_monomials':counts,'three_exact_U_F_constructions_equal':True})
    host,variances,deltas,route=states[frozenset()]
    for word,(h,v,_,_) in states.items():
        assert h==host
        for i in range(k):assert v[i]==(S.zeros(n) if i in word else variances[i])
    assert (12-1)-3*4==-1 and 2-3==-1
    return {'id':raw['id'],'children':children,'host_edges':outside,'host_loops':L,'states':records,
            'host_routes':[[int(v) for v in route[e][:,k:]] for e in outside],
            'child_pair_differences':[[[int(v) for v in d] for d in ds] for ds in deltas],
            'child_pair_indices':list(combinations(range(3),2)),
            'host_unchanged_in_every_forest':True,'each_child_variance_independent_of_other_child_T0':True,
            'all_child_mass_factors_retained':True,'child_simplex_volumes':['1/2']*k,
            'common_scale_power':-1,'relative_child_scale_powers':[-1]*k,
            'complete_log_kernel':'sum_forest -sign_forest*log1p(Q_active/(Mhost+sum_r))/U_active^3',
            'complete_integral_dimension':11,'whole_integral_evaluated':False}
