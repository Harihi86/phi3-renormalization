"""Retain a nested vertex T0 while exactly eliminating both primitive children."""
from itertools import combinations,product
import sympy as S
import primitive_triangle_factor_v3 as P
from nested_bubble_disjoint_triangle_factor_v1 import convert


def recognize(gd,z):
    parts=gd['proper_subgraphs'];assert len(parts)==3
    parents=[i for i,p in enumerate(parts) if (p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==(2,3,0,6)]
    assert len(parents)==1;parent=parents[0];pe=parts[parent]['edge_ids']
    inside=[i for i,p in enumerate(parts) if set(p['edge_ids'])<set(pe)];assert len(inside)==1;child=inside[0]
    assert (parts[child]['loops'],parts[child]['external_legs'],parts[child]['omega'],len(parts[child]['edge_ids']))==(1,3,0,3)
    other=next(i for i in range(3) if i not in (parent,child));op=parts[other]
    assert (op['loops'],op['external_legs'],op['omega'],len(op['edge_ids'])) in ((1,2,2,2),(1,3,0,3))
    ce=parts[child]['edge_ids'];oe=op['edge_ids'];assert not set(pe)&set(oe)
    inner=sorted(set(pe)-set(ce));outer=sorted(set(range(12))-set(pe)-set(oe))
    assert len(inner)==3 and len(outer)==6-len(oe)
    routed=[P.scaled_routes(t,z) for t in gd['terms']]
    for choices in product(ce,inner,oe,outer):
        witnesses=[]
        for term,rows in zip(gd['terms'],routed):
            w=convert(rows,choices)
            if w is None:break
            rr=w['routes'];word=term['subgraphs_inner_first']
            if any(rr[e][0] or rr[e][1] or rr[e][2] or rr[e][3] not in (-1,1) or any(v.has(*z) for v in rr[e]) for e in outer):break
            if any(rr[e][0] or rr[e][2] or rr[e][1] not in (-1,1) for e in inner):break
            if any(rr[e][0] not in (-1,1) or rr[e][2] for e in ce):break
            if any(rr[e][0] or rr[e][1] or rr[e][2] not in (-1,1) for e in oe):break
            if parent in word and any(v for e in inner for v in rr[e][:,3:]):break
            witnesses.append(w)
        if len(witnesses)==len(routed):return parent,child,other,pe,ce,oe,inner,outer,choices,witnesses
    raise AssertionError('No common original-edge basis respects the nested parent and disjoint child')


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]) and list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    parts=gd['proper_subgraphs'];z=S.symbols('z:3');alpha=S.symbols('a:12')
    parent,child,other,pe,ce,oe,inner,outer,choices,witnesses=recognize(gd,z)
    assert len(gd['terms'])==len(pd['states'])==8
    words={frozenset(i for i in range(3) if mask&(1<<i)) for mask in range(8)}
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==words
    lc,lo=sum(alpha[e] for e in ce),sum(alpha[e] for e in oe);combined=lc*lo
    # Retained order is outer q,inner k,p1,p2. Removed loops are the primitive children.
    kept=[3,1,4,5];removed=[0,2];records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        assert child not in word or parent not in word or word.index(child)<word.index(parent)
        active=[other] if other in word and parts[other]['omega']==2 else []
        assert term['jet_subgraphs']==state['variables']==active and term['jet_orders']==[2]*len(active)
        rr=w['routes'];matrix=sum((a*r.T*r for a,r in zip(alpha,rr)),S.zeros(6))
        assert matrix.extract(removed,removed)==S.diag(lc,lo)
        rows=[r[:,kept] for r in rr]
        H=sum((alpha[e]*rows[e].T*rows[e] for e in outer),S.zeros(4))
        N=sum((alpha[e]*rows[e].T*rows[e] for e in inner),S.zeros(4))
        variances=[];deltas=[]
        for edges,column in ((ce,0),(oe,2)):
            shifts=[rows[e]/rr[e][column] for e in edges];pairs=list(combinations(range(len(edges)),2))
            diff=[shifts[a]-shifts[b] for a,b in pairs]
            variance=sum((alpha[edges[a]]*alpha[edges[b]]*d.T*d for (a,b),d in zip(pairs,diff)),S.zeros(4)).applyfunc(S.expand)
            variances.append(variance);deltas.append(diff)
        vc,vo=variances;assert all(vo[1,j]==vo[j,1]==0 for j in range(4))
        if child in word:assert vc==S.zeros(4)
        if other in word and parts[other]['omega']==0:assert vo==S.zeros(4)
        reduced=combined*(H+N)+lo*vc+lc*vo
        direct=combined*matrix.extract(kept,kept)-lo*matrix.extract(kept,[0])*matrix.extract([0],kept)-lc*matrix.extract(kept,[2])*matrix.extract([2],kept)
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(4)
        ru=S.Poly(reduced[:2,:2].det(method='domain-ge'),*alpha,*z).exquo(S.Poly(combined,*alpha,*z)).as_expr()
        fn=sum(gram[i,j]*reduced.extract([0,1,2+i],[0,1,2+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        rf=S.Poly(fn,*alpha,*z).exquo(S.Poly(combined**2,*alpha,*z)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(6))
        ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,old,new in (('U',ug,ru),('F',fg,rf)):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(old-positive)==0 and S.expand(new-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states[frozenset(word)]=(H,N,vc,vo,rows,deltas)
        records.append({'word':word,'sign':term['sign'],'jet_orders':term['jet_orders'],'basis_edges':list(choices),
            'jacobian_determinant':int(w['det']),'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],'routing':[[str(v) for v in r] for r in rr],
            'U_F_monomials':counts,'gaussian_positive_Schur_equal':True})
    H,N,vc,vo,rows,deltas=states[frozenset()];parent_zero=states[frozenset([parent])]
    for word,(h,n,c,o,_,_) in states.items():
        assert h==H and n==(parent_zero[1] if parent in word else N)
        expected_c=S.zeros(4) if child in word else parent_zero[2] if parent in word else vc
        assert c==expected_c
        scale=(z[other]**2 if other in word else 1) if parts[other]['omega']==2 else (0 if other in word else 1)
        assert (o-scale*vo).applyfunc(S.expand)==S.zeros(4)
    as_int=lambda items:[[int(v) for v in row] for row in items]
    bare_inner=as_int(rows[e] for e in inner);zero_inner=as_int(parent_zero[4][e] for e in inner)
    assert all(r[0]==r[2]==r[3]==0 and abs(r[1])==1 for r in zero_inner)
    zero_child=as_int(parent_zero[5][0]);assert sorted(r[1]*r[1] for r in zero_child)==[0,1,1]
    assert all(r[0]==r[2]==r[3]==0 for r in zero_child)
    outside=as_int(rows[e] for e in outer);other_routes=as_int(deltas[1]);absorbed=None
    if parts[other]['omega']==2:
        assert len(other_routes)==1;bubble=other_routes[0]
        hits=[i for i,r in enumerate(outside) if r==bubble or r==[-v for v in bubble]];assert hits
        absorbed=outer[hits[0]];outside=[r for i,r in enumerate(outside) if i!=hits[0]]
    assert len(outside)==3 and all(r[1]==0 and abs(r[0])==1 for r in outside)
    return {'id':raw['id'],'parent_index':parent,'inner_triangle_index':child,'disjoint_child_index':other,
            'disjoint_child_kind':'bubble' if parts[other]['omega']==2 else 'triangle','states':records,
            'parent_edges':pe,'inner_triangle_edges':ce,'disjoint_child_edges':oe,'inner_host_edges':inner,
            'all_outer_host_edges':outer,'all_outer_host_routes':as_int(rows[e] for e in outer),
            'outer_host_routes':outside,'absorbed_bubble_host_edge':absorbed,'inner_host_routes':bare_inner,
            'parent_T0_inner_host_routes':zero_inner,'inner_triangle_routes':as_int(deltas[0]),
            'parent_T0_inner_triangle_routes':zero_child,'disjoint_child_routes':other_routes,
            'retained_loop_order':['outer q','inner k'],'local_parent_T0_retains_inner_k':True,
            'local_parent_T0_kernel':'Delta(k)^3*(-E(k²)/2)',
            'complete_kernel':'(1-Tphysical0) Houter*OuterChild*(Hinner*(-L(inner invariants))-Delta(k)^3*(-E(k²)/2))',
            'physical_T0_retains_both_loops':True,'all_eight_original_forests_retained':True,
            'primitive_children_integrated':2,'candidate_dimension':7,'integrability_pending':True,'whole_graph_acceptance':False}
