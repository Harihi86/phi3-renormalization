"""Exact eight-forest factorization of two bubbles and one disjoint triangle."""
from itertools import combinations,product
import sympy as S
import primitive_triangle_factor_v3 as P


def basis(gd,z):
    children=[p['edge_ids'] for p in gd['proper_subgraphs']];outside=[e for e in range(12) if all(e not in edges for edges in children)]
    routed=[P.scaled_routes(t,z) for t in gd['terms']]
    for choices in product(*children,outside):
        records=[]
        for rows in routed:
            change=S.Matrix.vstack(*(rows[e][:,:4] for e in choices));det=S.factor(change.det())
            if det not in (-1,1):break
            offset=S.Matrix.vstack(*(rows[e][:,4:] for e in choices));inverse=change.inv();converted=[]
            for row in rows:
                loop=row[:,:4]*inverse;external=row[:,4:]-loop*offset
                converted.append(S.Matrix([[S.expand(v) for v in list(loop)+list(external)]]))
            if any(v for e in outside for v in converted[e][:,:3]):break
            if any(converted[e][3] not in (-1,1) or any(v.has(*z) for v in converted[e]) for e in outside):break
            if any(converted[e][i] not in (-1,1) or any(converted[e][j] for j in range(3) if j!=i)
                   for i,edges in enumerate(children) for e in edges):break
            records.append({'change':change,'offset':offset,'det':det,'routing':converted})
        if len(records)==len(routed):return choices,records
    raise AssertionError('No common edge basis isolates all three primitive child loops')


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    gram=S.Matrix([[S.Rational(x) for x in row] for row in raw['external_gram']]);assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    assert list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    parts=gd['proper_subgraphs'];assert len(parts)==3
    kinds=[(p['loops'],p['external_legs'],p['omega'],len(p['edge_ids'])) for p in parts]
    assert sorted(kinds)==[(1,2,2,2),(1,2,2,2),(1,3,0,3)]
    bubbles=[i for i,k in enumerate(kinds) if k[1]==2];triangle=next(i for i,k in enumerate(kinds) if k[1]==3)
    children=[p['edge_ids'] for p in parts];assert all(not set(a)&set(b) for a,b in combinations(children,2))
    wanted={frozenset(i for i in range(3) if mask&(1<<i)) for mask in range(8)}
    assert len(gd['terms'])==len(pd['states'])==8
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==wanted
    z=S.symbols('z:3');choices,witnesses=basis(gd,z);alpha=S.symbols('a:12')
    outside=[e for e in range(12) if all(e not in edges for edges in children)];assert len(outside)==5
    lambdas=[sum(alpha[e] for e in edges) for edges in children];combined=S.prod(lambdas);records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        active=[i for i in word if i in bubbles]
        assert term['jet_subgraphs']==state['variables']==active and term['jet_orders']==[2]*len(active)
        route=w['routing'];matrix=sum((a*r.T*r for a,r in zip(alpha,route)),S.zeros(6))
        host=sum((alpha[e]*route[e][:,3:].T*route[e][:,3:] for e in outside),S.zeros(3))
        assert matrix[:3,:3]==S.diag(*lambdas)
        variances=[];deltas=[]
        for i,edges in enumerate(children):
            shifts=[route[e][:,3:]/route[e][i] for e in edges];pairs=list(combinations(range(len(edges)),2))
            differences=[shifts[a]-shifts[b] for a,b in pairs]
            variance=sum((alpha[edges[a]]*alpha[edges[b]]*d.T*d for (a,b),d in zip(pairs,differences)),S.zeros(3))
            variances.append(variance.applyfunc(S.expand));deltas.append(differences)
        reduced=combined*host+sum((S.prod(lambdas[j] for j in range(3) if j!=i)*v for i,v in enumerate(variances)),S.zeros(3))
        direct=combined*matrix[3:,3:]-sum((S.prod(lambdas[j] for j in range(3) if j!=i)*matrix[3:,i:i+1]*matrix[i:i+1,3:] for i in range(3)),S.zeros(3))
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(3)
        reduced_u=reduced[0,0]
        fn=sum(gram[i,j]*reduced.extract([0,1+i],[0,1+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        reduced_f=S.Poly(fn,*alpha,*z).exquo(S.Poly(combined,*alpha,*z)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(6))
        ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,value,red in (('U',ug,reduced_u),('F',fg,reduced_f)):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(value-positive)==0 and S.expand(red-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states[frozenset(word)]=(host,variances,deltas,route)
        records.append({'word':word,'sign':term['sign'],'jet_orders':term['jet_orders'],'basis_edges':list(choices),
            'jacobian_determinant':int(w['det']),'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],'routing':[[str(v) for v in row] for row in route],
            'U_F_monomials':counts,'gaussian_positive_Schur_equal':True})
    host,variances,deltas,route=states[frozenset()]
    for word,(h,v,_,_) in states.items():
        assert h==host
        for i in range(3):
            scale=(z[i]**2 if i in word else 1) if i in bubbles else (0 if i in word else 1)
            assert (v[i]-scale*variances[i]).applyfunc(S.expand)==S.zeros(3)
    bubble_rows=[[int(v) for v in deltas[i][0]] for i in bubbles]
    triangle_rows=[[int(v) for v in row] for row in deltas[triangle]]
    assert all(abs(row[0])==1 for row in bubble_rows)
    assert sorted(row[0]*row[0] for row in triangle_rows)==[0,1,1]
    fixed=[row for row in triangle_rows if row[0]==0];assert len(fixed)==1
    assert fixed[0][1]**2+fixed[0][2]**2-fixed[0][1]*fixed[0][2]==1
    triangle_rows=[row for row in triangle_rows if row[0]!=0]+fixed
    hosts=[[int(v) for v in route[e][:,3:]] for e in outside]
    same_norm=lambda a,b:a==b or a==[-v for v in b]
    unused=list(range(5));removed=[]
    for row in bubble_rows:
        hits=[i for i in unused if same_norm(hosts[i],row)];assert hits,'Bubble requires its complete retained host denominator'
        removed.append(hits[0]);unused.remove(hits[0])
    assert len(unused)==3
    return {'id':raw['id'],'children':children,'bubble_indices':bubbles,'triangle_index':triangle,'host_loops':1,
            'states':records,'bubble_routes':bubble_rows,'triangle_pair_differences':triangle_rows,
            'all_host_edges':outside,'all_host_routes':hosts,'host_routes':[hosts[i] for i in unused],
            'absorbed_bubble_host_edges':[outside[i] for i in removed],
            'complete_bubble_factor':'J(t)=B(t)/(1+t),B=int[(1+x(1-x)t)log(1+x(1-x)t)-x(1-x)t]dx',
            'complete_triangle_factor':'-L(a,b,1), including simplex volume1/2',
            'two_local_T2_projections_retained':True,'overall_T0_retains_host_loop':True,'whole_graph_acceptance':False}
