"""Exact original-state binding of an analytically integrated primitive vertex.

A selected triangle is eliminated as one Gaussian loop. Pairing its bare and
zero-degree Taylor states produces -log1p of a positive simplex quadratic.
Other forest variables remain symbolic; no outer propagator is projected away.
This module proves a kernel reduction, not a completed numerical integral.
"""
from itertools import combinations
import sympy as s


def scaled_routes(term,z):
    variables=term['jet_subgraphs']
    assert all(len(v['powers'])==len(variables) for row in term['scaled_routing'] for v in row)
    return [s.Matrix([[s.Rational(v['coefficient'])*s.prod(z[i]**p for i,p in zip(variables,v['powers'])) for v in row]])
            for row in term['scaled_routing']]


def find_triangle(gd):
    edges=gd['graph']['edges'];candidates=[]
    for i,part in enumerate(gd['proper_subgraphs']):
        if (part['loops'],part['external_legs'],len(part['edge_ids']))!=(1,3,3):continue
        selected=part['edge_ids'];vertices=set(v for e in selected for v in edges[e])
        if len(vertices)!=3 or any(sum(v in edges[e] for e in selected)!=2 for v in vertices):continue
        if any(set(selected)&set(other['edge_ids']) for j,other in enumerate(gd['proper_subgraphs']) if j!=i):continue
        candidates.append(i)
    assert len(candidates)==1,'Exactly one edge-disjoint primitive triangle required'
    return candidates[0]


def basis_all(gd,triangle,z):
    L=gd['loops'];selected=gd['proper_subgraphs'][triangle]['edge_ids']
    outside=[i for i in range(gd['edges']) if i not in selected]
    routes=[scaled_routes(term,z) for term in gd['terms']]
    for own in selected:
        for others in combinations(outside,L-1):
            choices=(own,*others);records=[]
            for rows in routes:
                change=s.Matrix.vstack(*(rows[e][:,:L] for e in choices));det=s.factor(change.det())
                if det not in (-1,1):break
                offset=s.Matrix.vstack(*(rows[e][:,L:] for e in choices));inv=change.inv();converted=[]
                for row in rows:
                    loop=row[:,:L]*inv;external=row[:,L:]-loop*offset
                    converted.append(s.Matrix([[s.expand(v) for v in list(loop)+list(external)]]))
                if any(row[0]!=0 for i,row in enumerate(converted) if i in outside):break
                if any(converted[e][0] not in (-1,1) for e in selected):break
                records.append({'change':change,'offset':offset,'det':det,'routing':converted})
            if len(records)==len(routes):return choices,records
    raise AssertionError('No common original-edge unimodular loop basis')


def positive_polynomial(rows,state,alpha,z):
    result=0
    for row in rows:
        mask,slot,channel=row[:3];co=s.Rational(row[3] if len(row)==4 else 1)
        result+=co*s.prod(alpha[e] for e in range(len(alpha)) if mask&(1<<e))*s.prod(z[v]**2 for j,v in enumerate(state['variables']) if slot&(1<<j))
    return s.expand(result)


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert gd['dimension']==6 and gd['external_legs']==3 and gd['overall_degree']==0
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert all(v==1 for v in raw['positive_invariants'])
    parts=gd['proper_subgraphs'];tri=find_triangle(gd);tri_edges=parts[tri]['edge_ids']
    assert parts[tri]['omega']==0
    L=gd['loops'];alpha=s.symbols('a:'+str(gd['edges']));z=s.symbols('z:'+str(len(parts)))
    gram=s.Matrix([[s.Rational(v) for v in row] for row in raw['external_gram']]);assert gram.shape==(2,2)
    words={tuple(term['subgraphs_inner_first']) for term in gd['terms']}
    assert len(words)==len(gd['terms'])==len(pd['states'])
    assert {tuple(state['word']) for state in pd['states']}==words
    choices,witnesses=basis_all(gd,tri,z);records=[];matrices={}
    lam=sum(alpha[e] for e in tri_edges)
    for term,witness in zip(gd['terms'],witnesses):
        word=tuple(term['subgraphs_inner_first']);state=next(v for v in pd['states'] if tuple(v['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        active={i for i in word if parts[i]['omega']>0}
        assert set(term['jet_subgraphs'])==active==set(state['variables'])
        assert term['jet_orders']==[parts[i]['omega'] for i in term['jet_subgraphs']]
        route=witness['routing'];M=s.zeros(L+2);host=s.zeros(L+1);vectors=[]
        for e,row in enumerate(route):
            M+=alpha[e]*row.T*row
            if e not in tri_edges:host+=alpha[e]*row[:,1:].T*row[:,1:]
        for e in tri_edges:vectors.append(route[e][:,1:]/route[e][0])
        variance=s.zeros(L+1)
        for i,j in combinations(range(3),2):
            delta=vectors[i]-vectors[j];variance+=alpha[tri_edges[i]]*alpha[tri_edges[j]]*delta.T*delta
        variance=variance.applyfunc(s.expand)
        schur=(lam*M[1:,1:]-M[1:,:1]*M[:1,1:]).applyfunc(s.expand)
        assert s.expand(M[0,0]-lam)==0 and (schur-lam*host-variance).applyfunc(s.expand)==s.zeros(L+1),'Gaussian completion mismatch'
        if tri in word:
            assert variance==s.zeros(L+1),'Triangle T0 retained a boundary component'
            assert not any(v.has(z[tri]) for v in host),'Triangle Taylor changed host propagators'
        original=s.zeros(L+2)
        for a,row in zip(alpha,scaled_routes(term,z)):original+=a*row.T*row
        ug=s.expand(original[:L,:L].det(method='domain-ge'))
        fg=s.expand(sum(gram[i,j]*original.extract([*range(L),L+i],[*range(L),L+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        # Fraction-free Schur determinants: exact divisibility is checked.
        reduced_u=s.cancel(schur[:L-1,:L-1].det(method='domain-ge')/lam**(L-2))
        reduced_f=s.cancel(sum(gram[i,j]*schur.extract([*range(L-1),L-1+i],[*range(L-1),L-1+j]).det(method='domain-ge') for i in range(2) for j in range(2))/lam**(L-1))
        counts=[]
        for name,gaussian,reduced in (('U',ug,reduced_u),('F',fg,reduced_f)):
            positive=positive_polynomial(state[name],state,alpha,z)
            assert s.expand(gaussian-positive)==0 and s.expand(reduced-positive)==0,(name,word)
            counts.append(len(s.Poly(positive,*alpha,*z).terms()))
        matrices[word]=(host,variance)
        records.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
                        'jacobian_determinant':int(witness['det']),'basis_edges':list(choices),
                        'loop_change':[[str(v) for v in witness['change'].row(i)] for i in range(L)],
                        'external_offset':[[str(v) for v in witness['offset'].row(i)] for i in range(L)],
                        'original_edge_routing':[[str(v) for v in row] for row in route],
                        'triangle_shift_vectors':[[str(v) for v in row] for row in vectors],
                        'U_F_monomials':counts,'all_three_polynomial_constructions_equal':True})
    pairs=[]
    for word,(host,variance) in matrices.items():
        if tri in word:continue
        partner=next((w for w in words if set(w)==set(word)|{tri}),None);assert partner is not None
        h2,v2=matrices[partner]
        assert (host-h2).applyfunc(s.expand)==s.zeros(L+1)
        assert variance!=s.zeros(L+1) and v2==s.zeros(L+1)
        pairs.append({'unsubtracted_word':list(word),'triangle_projected_word':list(partner),
                      'host_unchanged_exactly':True,'projected_boundary_variance_zero_exactly':True,
                      'unprojected_variance_nonzero':True,'triangle_T0_already_evaluated_in_original_serialization':True,
                      'primitive_integral':'-integral_simplex log(1+sum_i<j x_i*x_j*|a_i-a_j|^2)',
                      'local_mass_squared':1,'simplex_measure':'dx0 dx1; x2=1-x0-x1; volume=1/2',
                      'no_extra_factorial_or_symmetry_weight':True})
    assert 2*len(pairs)==len(records)
    return {'triangle_part_index':tri,'triangle_original_edges':tri_edges,'states':records,'paired_subtractions':pairs,
            'subgraph_loop_integrated_analytically':1,'all_external_and_outer_loop_components_retained':True,
            'whole_integral_not_computed':True,'overall_three_point_T0_still_required':True}
