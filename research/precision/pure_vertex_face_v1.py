"""Exact projective regrouping for disjoint/nested, arbitrary-loop T0 vertex parts."""
from fractions import Fraction as Q
from itertools import combinations
import nested_triangle_projective_v1 as N
import projective_pair_face_transport_v2 as V
H=N.H
minors=N.minors
positive=N.positive
polynomial=N.polynomial
gaussian_uf=N.gaussian_uf
prepare_pairs=N.prepare_pairs
verify_table=N.verify_table


def compile_original(raw,recompiled):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert raw['external_gram']==[[1,-.5],[-.5,1]] and raw['positive_invariants']==[1,1,1]
    assert recompiled['terms']==gd['terms'] and recompiled['proper_subgraphs']==gd['proper_subgraphs']
    parts=gd['proper_subgraphs'];assert 1<=len(parts)<=3
    assert all(p['external_legs']==3 and p['omega']==0 and 1<=p['loops']<=3
               and len(p['edge_ids'])==3*p['loops'] for p in parts)
    edge_sets=[set(p['edge_ids']) for p in parts]
    assert all(not a&b or a<b or b<a for a,b in combinations(edge_sets,2))
    indices=sorted(range(len(parts)),key=lambda i:(len(parts[i]['edge_ids']),i))
    words={tuple(i for i in indices if mask&(1<<i)) for mask in range(1<<len(parts))}
    assert len(gd['terms'])==len(pd['states'])==len(words)
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']}==words
    assert {tuple(s['word']) for s in pd['states']}==words
    states=[]
    for term in gd['terms']:
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if s['word']==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==term['jet_orders']==state['variables']==[]
        assert all(v['powers']==[] and type(v['coefficient']) is int for row in term['scaled_routing'] for v in row)
        routes=[[v['coefficient'] for v in row] for row in term['scaled_routing']]
        u,f=minors(routes)
        assert u==positive(state['U']) and f==positive(state['F'])
        assert all(m.bit_count()==4 for m,c in u) and all(m.bit_count()==5 for m,c in f)
        states.append({'word':word,'sign':term['sign'],'routes':routes,'U':u,'F':f})
    selected=indices[0];pairs=[]
    for index,state in enumerate(states):
        if selected in state['word']:continue
        target=sorted([*state['word'],selected],key=lambda i:indices.index(i))
        projected=next(i for i,s in enumerate(states) if s['word']==target)
        assert states[projected]['sign']==-state['sign'];pairs.append([index,projected])
    assert len(pairs)*2==len(states)
    return {'id':raw['id'],'states':states,'pairs':pairs,'pairing_part_index':selected,
        'proper_parts':[{'loops':p['loops'],'external_legs':p['external_legs'],'omega':p['omega'],'edge_ids':p['edge_ids']} for p in parts],
        'original_edges':12,'original_loops':4,'physical_dimension':11,'sampling_dimension':12,
        'complete_density_homogeneity':-12,'fresh_graph_BPH_equals_original':True,
        'gaussian_Cauchy_Binet_equals_positive_forests':True,'all_overall_and_proper_T0_retained':True,
        'added_symmetry_or_factorial':False,'radial_kernel':'-log1p(F/(M*U))/U^3; M=sum(alpha), m^2=1',
        'paired_kernel':'s*((1/Uz^3-1/Ub^3)*log1p(rz)-log1p((rb-rz)/(1+rz))/Ub^3)'}


def face_tables(binding):
    proof=V.construct(binding,False)
    assert min(proof['gamma'].values())>=1,'all coordinate-face orders must be proved positive'
    tables={};checks={}
    for cap in (2,3):
        table=H.build(12,{s:min(g,cap) for s,g in proof['gamma'].items()},2)
        assert min(table['shapes'].values())>=1 and max(table['shapes'].values())<=cap
        tables[str(cap)]=table;checks[str(cap)]=verify_table(table)
    return tables,proof,checks
