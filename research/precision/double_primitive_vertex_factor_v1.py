"""Two disjoint primitives inside one three-loop vertex: exact eight-forest reduction."""
from itertools import combinations,product
import sympy as S
import primitive_triangle_factor_v3 as P
from nested_bubble_disjoint_triangle_factor_v1 import convert


def recognize(gd,z):
    parts=gd['proper_subgraphs'];assert len(parts)==3
    parents=[i for i,p in enumerate(parts) if (p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==(3,3,0,9)]
    assert len(parents)==1;parent=parents[0];pe=parts[parent]['edge_ids'];children=[i for i in range(3) if i!=parent]
    ce=[parts[i]['edge_ids'] for i in children]
    assert all((parts[i]['loops'],parts[i]['external_legs'],parts[i]['omega'],len(parts[i]['edge_ids'])) in ((1,2,2,2),(1,3,0,3)) for i in children)
    assert all(set(edges)<set(pe) for edges in ce) and not set(ce[0])&set(ce[1])
    inner=sorted(set(pe)-set(ce[0])-set(ce[1]));outer=sorted(set(range(12))-set(pe));assert len(outer)==3
    routed=[P.scaled_routes(t,z) for t in gd['terms']]
    for choices in product(*ce,inner,outer):
        witnesses=[]
        for term,rows in zip(gd['terms'],routed):
            w=convert(rows,choices)
            if w is None:break
            rr=w['routes'];word=term['subgraphs_inner_first']
            if any(any(rr[e][j] for j in range(3)) or rr[e][3] not in (-1,1) or any(v.has(*z) for v in rr[e]) for e in outer):break
            if any(rr[e][0] or rr[e][1] or rr[e][2] not in (-1,1) for e in inner):break
            if any(rr[e][j] not in (-1,1) or rr[e][1-j] for j,edges in enumerate(ce) for e in edges):break
            if parent in word and any(v for e in inner for v in rr[e][:,3:]):break
            witnesses.append(w)
        if len(witnesses)==len(routed):return parent,children,pe,ce,inner,outer,choices,witnesses
    raise AssertionError('No original-edge basis isolates both primitives and retains parent inner k')


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]) and list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    parts=gd['proper_subgraphs'];z=S.symbols('z:3');alpha=S.symbols('a:12')
    parent,children,pe,ce,inner,outer,choices,witnesses=recognize(gd,z)
    assert len(gd['terms'])==len(pd['states'])==8
    words={frozenset(i for i in range(3) if mask&(1<<i)) for mask in range(8)}
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==words
    lam=[sum(alpha[e] for e in edges) for edges in ce];combined=S.prod(lam)
    kept=[3,2,4,5];removed=[0,1];records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        assert all(i not in word or parent not in word or word.index(i)<word.index(parent) for i in children)
        active=[i for i in children if i in word and parts[i]['omega']==2]
        assert term['jet_subgraphs']==state['variables']==active and term['jet_orders']==[2]*len(active)
        rr=w['routes'];matrix=sum((a*r.T*r for a,r in zip(alpha,rr)),S.zeros(6))
        assert matrix.extract(removed,removed)==S.diag(*lam)
        rows=[r[:,kept] for r in rr]
        H=sum((alpha[e]*rows[e].T*rows[e] for e in outer),S.zeros(4));N=sum((alpha[e]*rows[e].T*rows[e] for e in inner),S.zeros(4))
        variances=[];deltas=[]
        for j,edges in enumerate(ce):
            shifts=[rows[e]/rr[e][j] for e in edges];pairs=list(combinations(range(len(edges)),2));diff=[shifts[a]-shifts[b] for a,b in pairs]
            variance=sum((alpha[edges[a]]*alpha[edges[b]]*d.T*d for (a,b),d in zip(pairs,diff)),S.zeros(4)).applyfunc(S.expand)
            variances.append(variance);deltas.append(diff)
        reduced=combined*(H+N)+sum((lam[1-j]*v for j,v in enumerate(variances)),S.zeros(4))
        direct=combined*matrix.extract(kept,kept)-sum((lam[1-j]*matrix.extract(kept,[j])*matrix.extract([j],kept) for j in range(2)),S.zeros(4))
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(4)
        ru=S.Poly(reduced[:2,:2].det(method='domain-ge'),*alpha,*z).exquo(S.Poly(combined,*alpha,*z)).as_expr()
        fn=sum(gram[i,j]*reduced.extract([0,1,2+i],[0,1,2+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        rf=S.Poly(fn,*alpha,*z).exquo(S.Poly(combined**2,*alpha,*z)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(6));ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,old,new in (('U',ug,ru),('F',fg,rf)):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(old-positive)==0 and S.expand(new-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states[frozenset(word)]=(H,N,variances,rows,deltas)
        records.append({'word':word,'sign':term['sign'],'jet_orders':term['jet_orders'],'basis_edges':list(choices),
            'jacobian_determinant':int(w['det']),'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],'routing':[[str(v) for v in r] for r in rr],
            'U_F_monomials':counts,'gaussian_positive_Schur_equal':True})
    H,N,variances,rows,deltas=states[frozenset()];zero=states[frozenset([parent])]
    for word,(h,n,vs,_,_) in states.items():
        assert h==H and n==(zero[1] if parent in word else N)
        for j,i in enumerate(children):
            target=zero[2][j] if parent in word else variances[j]
            scale=(z[i]**2 if i in word else 1) if parts[i]['omega']==2 else (0 if i in word else 1)
            assert (vs[j]-scale*target).applyfunc(S.expand)==S.zeros(4)
    as_int=lambda items:[[int(v) for v in row] for row in items]
    zero_inner=as_int(zero[3][e] for e in inner)
    assert all(r[0]==r[2]==r[3]==0 and abs(r[1])==1 for r in zero_inner)
    child_rows=[as_int(ds) for ds in deltas];zero_children=[as_int(ds) for ds in zero[4]]
    remaining=list(inner);absorbed=[];kinds=[]
    for j,i in enumerate(children):
        is_bubble=parts[i]['omega']==2;kinds.append('bubble' if is_bubble else 'triangle')
        assert all(r[0]==r[2]==r[3]==0 for r in zero_children[j])
        if is_bubble:
            assert len(child_rows[j])==1 and abs(zero_children[j][0][1])==1
            bubble=child_rows[j][0];hits=[e for e in remaining if list(rows[e])==bubble or list(rows[e])==[-v for v in bubble]]
            assert hits,'Each original bubble needs its own matching parent host denominator'
            absorbed.append(hits[0]);remaining.remove(hits[0])
        else:assert sorted(r[1]*r[1] for r in zero_children[j])==[0,1,1]
    assert len(remaining)==3
    return {'id':raw['id'],'parent_index':parent,'child_indices':children,'child_kinds':kinds,'states':records,
        'parent_edges':pe,'child_edges':ce,'all_inner_host_edges':inner,'all_inner_host_routes':as_int(rows[e] for e in inner),
        'inner_host_routes':as_int(rows[e] for e in remaining),'absorbed_bubble_host_edges':absorbed,
        'outer_host_edges':outer,'outer_host_routes':as_int(rows[e] for e in outer),'child_routes':child_rows,
        'parent_T0_inner_host_routes':zero_inner,'parent_T0_child_routes':zero_children,
        'parent_zero_kernel':'Delta(k)^3 * product(J(k²) for bubble, -E(k²)/2 for triangle)',
        'complete_kernel':'(1-Tphysical0) Houter * (Hinner*O1*O2-parent_zero_kernel)',
        'local_parent_T0_retains_inner_k':True,'all_eight_original_forests_retained':True,'primitive_children_integrated':2,
        'candidate_dimension':7,'integrability_and_native_tests_pending':True,'whole_graph_acceptance':False}
