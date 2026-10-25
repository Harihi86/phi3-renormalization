"""Exact bubble/two-loop-child reduction, optionally inside a three-loop parent."""
from itertools import product,combinations
import sympy as S
import primitive_triangle_factor_v3 as P
from nested_bubble_disjoint_triangle_factor_v1 import convert

def recognize(gd,z):
    parts=gd['proper_subgraphs'];assert len(parts) in (2,3)
    def unique(kind):
        hits=[i for i,p in enumerate(parts) if (p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==kind]
        assert len(hits)==1;return hits[0]
    bub=unique((1,2,2,2));child=unique((2,3,0,6));parent=unique((3,3,0,9)) if len(parts)==3 else None
    be,ce=[parts[i]['edge_ids'] for i in (bub,child)];assert set(be)<set(ce)
    inner=sorted(set(ce)-set(be));host=sorted(set(range(12))-set(ce));assert len(inner)==4 and len(host)==6
    if parent is not None:
        pe=parts[parent]['edge_ids'];assert set(ce)<set(pe)
        middle=sorted(set(pe)-set(ce));outer=sorted(set(range(12))-set(pe));assert len(middle)==len(outer)==3
        choices_list=product(be,inner,middle,outer)
    else:
        middle=[];outer=host;choices_list=((a,b,c,d) for a,b in product(be,inner) for c,d in combinations(host,2))
    routed=[P.scaled_routes(t,z) for t in gd['terms']]
    for choices in choices_list:
        witnesses=[]
        for term,rows in zip(gd['terms'],routed):
            w=convert(rows,choices)
            if w is None:break
            rr=w['routes'];word=term['subgraphs_inner_first'];scale=z[bub] if bub in word else S.Integer(1)
            if any(rr[e][0] or rr[e][1] or any(v.has(*z) for v in rr[e]) for e in host):break
            if parent is not None and any(rr[e][2] or rr[e][3] not in (-1,1) for e in outer):break
            if parent is not None and any(rr[e][2] not in (-1,1) for e in middle):break
            if parent in word and any(v for e in middle for v in rr[e][:,3:]):break
            if any(rr[e][0] or rr[e][1] not in (-1,1) for e in inner):break
            if any(rr[e][0] not in (-1,1) for e in be):break
            shifts=[rr[e][:,1:]/rr[e][0] for e in be];delta=shifts[0]-shifts[1]
            if delta[0] not in (-scale,scale) or any(delta[j] for j in range(1,5)):break
            normalized=[rr[e][:,2:]/rr[e][1] for e in inner]
            if child in word:
                if any(any(v for v in row) for row in normalized):break
            elif not word and sum(row==S.zeros(1,4) for row in normalized)!=2:break
            if parent in word and any(v for row in normalized for v in row[:,1:]):break
            witnesses.append(w)
        if len(witnesses)==len(routed):return bub,child,parent,be,ce,inner,host,middle,outer,choices,witnesses
    raise AssertionError('No original-edge basis isolates the complete nested bubble child')

def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]) and list(map(S.Rational,raw['positive_invariants']))==[1,1,1]
    parts=gd['proper_subgraphs'];z=S.symbols('z:'+str(len(parts)));alpha=S.symbols('a:12')
    bub,child,parent,be,ce,inner,host,middle,outer,choices,witnesses=recognize(gd,z)
    expected={frozenset(i for i in range(len(parts)) if mask&(1<<i)) for mask in range(2**len(parts))}
    assert len(gd['terms'])==len(pd['states'])==len(expected)
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(s['word']) for s in pd['states']}==expected
    lam=sum(alpha[e] for e in be);kept=[1,3,2,4,5];records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==state['sign']==(-1)**len(word)
        order=[i for i in (bub,child,parent) if i in word];assert word==order
        active=[bub] if bub in word else [];assert term['jet_subgraphs']==state['variables']==active and term['jet_orders']==[2]*len(active)
        rr=w['routes'];matrix=sum((a*r.T*r for a,r in zip(alpha,rr)),S.zeros(6));assert matrix[0,0]==lam
        rows=[r[:,kept] for r in rr]
        H=sum((alpha[e]*rows[e].T*rows[e] for e in host),S.zeros(5));N=sum((alpha[e]*rows[e].T*rows[e] for e in inner),S.zeros(5))
        shifts=[rows[e]/rr[e][0] for e in be];d=shifts[0]-shifts[1];variance=(alpha[be[0]]*alpha[be[1]]*d.T*d).applyfunc(S.expand)
        wanted=S.zeros(5);wanted[0,0]=alpha[be[0]]*alpha[be[1]]*(z[bub]**2 if bub in word else 1);assert variance==wanted
        reduced=lam*(H+N)+variance;direct=lam*matrix.extract(kept,kept)-matrix.extract(kept,[0])*matrix.extract([0],kept)
        assert (direct-reduced).applyfunc(S.expand)==S.zeros(5)
        ru=S.Poly(reduced[:3,:3].det(method='domain-ge'),*alpha,*z).exquo(S.Poly(lam**2,*alpha,*z)).as_expr()
        rf0=sum(gram[i,j]*reduced.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        rf=S.Poly(rf0,*alpha,*z).exquo(S.Poly(lam**3,*alpha,*z)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(6));ug=original[:4,:4].det(method='domain-ge')
        fg=sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts=[]
        for name,old,new in (('U',ug,ru),('F',fg,rf)):
            positive=P.positive_polynomial(state[name],state,alpha,z)
            assert S.expand(old-positive)==0 and S.expand(new-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        normalized=[rows[e][:,1:]/rows[e][0] for e in inner]
        states[frozenset(word)]=(H,N,variance,rows,normalized)
        records.append({'word':word,'sign':term['sign'],'jet_orders':term['jet_orders'],'basis_edges':list(choices),
            'jacobian_determinant':int(w['det']),'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],'routing':[[str(v) for v in r] for r in rr],
            'U_F_monomials':counts,'gaussian_positive_Schur_equal':True})
    bare=states[frozenset()];zero=states[frozenset([child])]
    parent0=states[frozenset([parent])] if parent is not None else None
    for word,(h,n,v,rows,normalized) in states.items():
        assert h==(parent0[0] if parent in word else bare[0])
        assert n==(zero[1] if child in word else parent0[1] if parent in word else bare[1])
        assert v==bare[2]*(z[bub]**2 if bub in word else 1)
        assert all(rows[e]==bare[3][e] for e in outer)
    indices=[i for i,row in enumerate(bare[4]) if row!=S.zeros(1,4)];assert len(indices)==2
    shifted=[bare[4][i] for i in indices];as_int=lambda rows:[[int(v) for v in row] for row in rows]
    # Exact Feynman completion of the two repeated zero-shift propagators.
    x,y=S.symbols('x y');a,b=shifted;x1=(1-x)*y;x2=(1-x)*(1-y)
    A=S.Matrix([[1,*list(a)]]);B=S.Matrix([[1,*list(b)]]);K=S.Matrix([[1,0,0,0,0]])
    center=K+S.Matrix([[0,*list(x1*a+x2*b)]])
    variance=x*x1*a.T*a+x*x2*b.T*b+x1*x2*(a-b).T*(a-b)
    check=x*K.T*K+x1*A.T*A+x2*B.T*B-center.T*center
    expected_matrix=S.zeros(5);expected_matrix[1:,1:]=variance
    assert (check-expected_matrix).applyfunc(S.expand)==S.zeros(5)
    zero_routes=as_int(parent0[4][i] for i in indices) if parent0 else None
    if zero_routes is not None:assert all(r[0]==r[2]==r[3]==0 for r in zero_routes)
    return {'id':raw['id'],'bubble_index':bub,'child_index':child,'parent_index':parent,'states':records,
        'bubble_edges':be,'child_edges':ce,'inner_edges':inner,'host_edges':host,'middle_edges':middle,'outer_edges':outer,
        'child_boundary_shifts':as_int(shifted),'host_routes':as_int(bare[3][e][:,1:] for e in host),
        'middle_host_routes':as_int(bare[3][e][:,1:] for e in middle),'outer_host_routes':as_int(bare[3][e][:,1:] for e in outer),
        'parent_T0_child_boundary_shifts':zero_routes,
        'parent_T0_middle_host_routes':as_int(parent0[3][e][:,1:] for e in middle) if parent0 else None,
        'nested_bubble_boundary':'retained child inner h,not physical momentum or either host loop',
        'K':'-integral_((1-x)/x)^((1-x+V)/x) B(t)/(1+t)^2 dt','simplex_measure':'(1-x)dxdy',
        'child_variance':'x(1-x)y |a|²+x(1-x)(1-y)|b|²+(1-x)²y(1-y)|a-b|²',
        'complete_kernel':'(1-Tphysical0) Houter*(Hmiddle*K(x,V)-Delta(k)^3*K(x,Vparent0))' if parent0 else '(1-Tphysical0) Hall6*K(x,V)',
        'two_inner_loops_integrated_analytically':True,'candidate_dimension':9,'whole_graph_acceptance':False}
