"""Exact forests for a triangle inside a two-loop vertex with two retained hosts."""
from itertools import product,combinations
import sympy as S
import primitive_triangle_factor_v3 as P
from nested_bubble_disjoint_triangle_factor_v1 import convert

def recognize(gd):
    parts=gd['proper_subgraphs'];assert len(parts) in (2,3)
    def unique(kind):
        hits=[i for i,p in enumerate(parts) if (p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==kind]
        assert len(hits)==1;return hits[0]
    tri=unique((1,3,0,3));child=unique((2,3,0,6));parent=unique((3,3,0,9)) if len(parts)==3 else None
    te,ce=[parts[i]['edge_ids'] for i in (tri,child)];assert set(te)<set(ce)
    inner=sorted(set(ce)-set(te));host=sorted(set(range(12))-set(ce));assert len(inner)==3 and len(host)==6
    if parent is not None:
        pe=parts[parent]['edge_ids'];assert set(ce)<set(pe)
        middle=sorted(set(pe)-set(ce));outer=sorted(set(range(12))-set(pe));assert len(middle)==len(outer)==3
        choices_list=product(te,inner,middle,outer)
    else:
        middle=[];outer=host;choices_list=((a,b,c,d) for a,b in product(te,inner) for c,d in combinations(host,2))
    routed=[P.scaled_routes(t,[]) for t in gd['terms']]
    for choices in choices_list:
        witnesses=[]
        for term,rows in zip(gd['terms'],routed):
            w=convert(rows,choices)
            if w is None:break
            rr=w['routes'];word=term['subgraphs_inner_first']
            if any(rr[e][0] or rr[e][1] for e in host):break
            if parent is not None and any(rr[e][2] or rr[e][3] not in (-1,1) for e in outer):break
            if parent is not None and any(rr[e][2] not in (-1,1) for e in middle):break
            if parent in word and any(v for e in middle for v in rr[e][:,3:]):break
            if any(rr[e][0] or rr[e][1] not in (-1,1) for e in inner):break
            if any(rr[e][0] not in (-1,1) for e in te):break
            if child in word and any(v for e in inner for v in rr[e][:,2:]):break
            if parent in word and any(v for e in inner for v in rr[e][:,3:]):break
            witnesses.append(w)
        if len(witnesses)==len(routed):return tri,child,parent,te,ce,inner,host,middle,outer,choices,witnesses
    raise AssertionError('No original-edge basis isolates the complete nested triangle child')

def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]) and list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    parts=gd['proper_subgraphs'];alpha=S.symbols('a:12')
    tri,child,parent,te,ce,inner,host,middle,outer,choices,witnesses=recognize(gd)
    expected={frozenset(i for i in range(len(parts)) if mask&(1<<i)) for mask in range(2**len(parts))}
    assert len(gd['terms'])==len(pd['states'])==len(expected)
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==expected
    lam=sum(alpha[e] for e in te);kept=[1,3,2,4,5];records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        assert word==[i for i in (tri,child,parent) if i in word]
        assert term['jet_subgraphs']==term['jet_orders']==state['variables']==[]
        rr=w['routes'];matrix=sum((a*r.T*r for a,r in zip(alpha,rr)),S.zeros(6));assert matrix[0,0]==lam
        rows=[r[:,kept] for r in rr]
        H=sum((alpha[e]*rows[e].T*rows[e] for e in host),S.zeros(5));N=sum((alpha[e]*rows[e].T*rows[e] for e in inner),S.zeros(5))
        shifts=[rows[e]/rr[e][0] for e in te];pairs=list(combinations(range(3),2));diff=[shifts[a]-shifts[b] for a,b in pairs]
        variance=sum((alpha[te[a]]*alpha[te[b]]*d.T*d for (a,b),d in zip(pairs,diff)),S.zeros(5)).applyfunc(S.expand)
        if tri in word:assert variance==S.zeros(5)
        reduced=lam*(H+N)+variance
        direct=lam*matrix.extract(kept,kept)-matrix.extract(kept,[0])*matrix.extract([0],kept)
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(5)
        ru=S.Poly(reduced[:3,:3].det(method='domain-ge'),*alpha).exquo(S.Poly(lam**2,*alpha)).as_expr()
        rf0=sum(gram[i,j]*reduced.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        rf=S.Poly(rf0,*alpha).exquo(S.Poly(lam**3,*alpha)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,[]))),S.zeros(6));ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,old,new in (('U',ug,ru),('F',fg,rf)):
            positive=P.positive_polynomial(state[name],state,alpha,[])
            assert S.expand(old-positive)==0 and S.expand(new-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha).terms()))
        normalized=[rows[e][:,1:]/rows[e][0] for e in inner]
        states[frozenset(word)]=(H,N,variance,rows,normalized,diff)
        records.append({'word':word,'sign':term['sign'],'jet_orders':[], 'basis_edges':list(choices),
            'jacobian_determinant':int(w['det']),'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],'routing':[[str(v) for v in r] for r in rr],
            'U_F_monomials':counts,'gaussian_positive_Schur_equal':True})
    bare=states[frozenset()];zero=states[frozenset([child])]
    parent0=states[frozenset([parent])] if parent is not None else None
    for word,(h,n,v,rows,normalized,diff) in states.items():
        assert h==(parent0[0] if parent in word else bare[0])
        assert n==(zero[1] if child in word else parent0[1] if parent in word else bare[1])
        expected_v=S.zeros(5) if tri in word else zero[2] if child in word else parent0[2] if parent in word else bare[2]
        assert v==expected_v
        assert all(rows[e]==bare[3][e] for e in outer)
    assert all(row==S.zeros(1,4) for row in zero[4])
    assert all(all(v==0 for v in row[1:]) for row in zero[5]) and sorted(int(row[0]**2) for row in zero[5])==[0,1,1]
    assert [row[0] for row in bare[5]]==[row[0] for row in zero[5]]
    as_int=lambda rows:[[int(v) for v in row] for row in rows]
    return {'id':raw['id'],'triangle_index':tri,'child_index':child,'parent_index':parent,'states':records,
        'triangle_edges':te,'child_edges':ce,'inner_edges':inner,'host_edges':host,'middle_edges':middle,'outer_edges':outer,
        'child_inner_shifts':as_int(bare[4]),'child_triangle_routes':as_int(bare[5]),
        'host_routes':as_int(bare[3][e][:,1:] for e in host),
        'middle_host_routes':as_int(bare[3][e][:,1:] for e in middle),'outer_host_routes':as_int(bare[3][e][:,1:] for e in outer),
        'parent_T0_child_inner_shifts':as_int(parent0[4]) if parent0 else None,
        'parent_T0_child_triangle_routes':as_int(parent0[5]) if parent0 else None,
        'parent_T0_middle_host_routes':as_int(parent0[3][e][:,1:] for e in middle) if parent0 else None,
        'nested_triangle_boundary':'retained inner h and actual child boundaries; child T0 retains h',
        'two_inner_loops_integrated_analytically':True,'remaining_child_simplex_measure':'da0 da1 db0 db1; each simplex has volume1/2',
        'child_kernel':'R(q,ell,De,Ge), the complete common-scale and relative-scale logarithmic marginal',
        'complete_kernel':'(1-Tphysical0) Ho*(Hm*R-Delta(k)^3*Rparent0)' if parent0 else '(1-Tphysical0) Hall6*R',
        'candidate_dimension':11,'whole_graph_acceptance':False}
