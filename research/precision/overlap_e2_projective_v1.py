"""Original two-point forests, analytic common scale, and exact paired remainders.

Only internal three-point T0 subgraphs are accepted. The overall T2 projector
is integrated analytically. Overlapping subgraphs are never multiplied unless
the independently recompiled BPH program contains that forest.
"""
from fractions import Fraction as Q
from itertools import combinations
import math
import nested_triangle_projective_v1 as N
import exact_rational_sector_transport_v1 as H

EDGES=11


def minors(routes):
    assert len(routes)==EDGES and all(len(row)==5 for row in routes)
    u=[];f=[]
    for selected in combinations(range(EDGES),4):
        d=N.determinant([[routes[e][j] for j in range(4)] for e in selected])
        if d:u.append([sum(1<<e for e in selected),d*d])
    for selected in combinations(range(EDGES),5):
        d=N.determinant([routes[e] for e in selected])
        if d:f.append([sum(1<<e for e in selected),d*d])
    assert u and f
    return sorted(u),sorted(f)


def forest_pairs(states,parts):
    """Deterministic maximum matching on actual forest toggles, with leftovers.

    Matching is an algebraic regrouping, not an extra forest enumeration.
    A left/right augmenting-path search pairs opposite parities. Every edge
    differs by exactly one allowed Taylor subgraph. Unmatched terms remain.
    """
    words=[frozenset(s['word']) for s in states]
    assert len(set(words))==len(words)
    left=[i for i,s in enumerate(states) if s['sign']==1]
    right=[i for i,s in enumerate(states) if s['sign']==-1]
    neighbors={i:sorted((j for j in right if len(words[i]^words[j])==1),
        key=lambda j:(parts[next(iter(words[i]^words[j]))]['loops'],j)) for i in left}
    owner={}
    def assign(i,seen):
        for j in neighbors[i]:
            if j in seen:continue
            seen.add(j)
            if j not in owner or assign(owner[j],seen):owner[j]=i;return True
        return False
    for i in left:assign(i,set())
    pairs=[];used=set()
    for j,i in sorted(owner.items()):
        b,z=(i,j) if words[i]<words[j] else (j,i)
        assert words[b]<words[z] and len(words[z]-words[b])==1
        pairs.append([b,z]);used.update((b,z))
    pairs.sort();unpaired=[i for i in range(len(states)) if i not in used]
    assert sorted([i for pair in pairs for i in pair]+unpaired)==list(range(len(states)))
    return pairs,unpaired


def compile_original(raw,recompiled):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,2,EDGES,2)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert raw['external_gram']==[[1]] and raw['positive_invariants']==[1]
    assert recompiled['terms']==gd['terms'] and recompiled['proper_subgraphs']==gd['proper_subgraphs']
    parts=gd['proper_subgraphs'];assert parts
    assert all(p['external_legs']==3 and p['omega']==0 and len(p['edge_ids'])==3*p['loops'] for p in parts)
    assert len(gd['terms'])==len(pd['states'])
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']}=={tuple(s['word']) for s in pd['states']}
    states=[]
    for term in gd['terms']:
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if s['word']==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==term['jet_orders']==state['variables']==[]
        for a,b in combinations(word,2):
            x,y=set(parts[a]['edge_ids']),set(parts[b]['edge_ids'])
            assert x<y or y<x or not x&y,'overlapping subgraphs cannot coexist in a forest'
        assert all(v['powers']==[] and type(v['coefficient']) is int for row in term['scaled_routing'] for v in row)
        routes=[[v['coefficient'] for v in row] for row in term['scaled_routing']]
        u,f=minors(routes)
        assert u==N.positive(state['U']) and f==N.positive(state['F'])
        assert all(m.bit_count()==4 for m,c in u) and all(m.bit_count()==5 for m,c in f)
        states.append({'word':word,'sign':term['sign'],'routes':routes,'U':u,'F':f})
    pairs,unpaired=forest_pairs(states,parts)
    return {'id':raw['id'],'states':states,'pairs':pairs,'unpaired':unpaired,
        'original_edges':EDGES,'original_loops':4,'physical_dimension':10,'sampling_dimension':11,
        'complete_density_homogeneity':-11,'fresh_graph_BPH_equals_original':True,
        'gaussian_Cauchy_Binet_equals_positive_forests':True,
        'overall_T2_and_all_internal_T0_retained':True,'overlap_products_not_invented':True,
        'added_symmetry_or_factorial':False,'radial_kernel':'sign(w)*M*phi(F/(M*U))/U^3; phi(x)=(1+x)log1p(x)-x',
        'pairing':'maximum matching of actual opposite-sign forest toggles; unmatched states retained'}


def polynomial(terms,alpha):
    return sum((co*math.prod(alpha[e] for e in range(EDGES) if mask&(1<<e)) for mask,co in terms),Q())


def gaussian_uf(state,alpha):
    """Exact Gaussian Schur complement, without positive polynomial evaluation."""
    routes=state['routes']
    matrix=[[sum((a*r[i]*r[j] for a,r in zip(alpha,routes)),Q()) for j in range(5)] for i in range(5)]
    u=Q(1)
    for k in range(4):
        pivot=matrix[k][k];assert pivot>0;u*=pivot
        for i in range(k+1,5):
            for j in range(k+1,5):matrix[i][j]-=matrix[i][k]*matrix[k][j]/pivot
    return u,u*matrix[4][4]


def validate_pairs(binding):
    states=binding['states'];used=[]
    for b,z in binding['pairs']:
        assert 0<=b<len(states) and 0<=z<len(states)
        wb,wz=set(states[b]['word']),set(states[z]['word'])
        assert wb<wz and len(wz-wb)==1 and states[b]['sign']==-states[z]['sign']
        used.extend((b,z))
    used.extend(binding['unpaired'])
    assert sorted(used)==list(range(len(states))),'every original forest must occur exactly once'


def prepare_pairs(binding,alpha,use_gaussian=False):
    assert len(alpha)==EDGES and all(type(a) in (int,Q) and a>0 for a in alpha)
    validate_pairs(binding);mass=sum(alpha);values=[]
    for state in binding['states']:
        u,f=gaussian_uf(state,alpha) if use_gaussian else (polynomial(state['U'],alpha),polynomial(state['F'],alpha))
        assert u>0 and f>0;values.append((u,f/(mass*u)))
    records=[]
    for bare,zero in binding['pairs']:
        ub,rb=values[bare];uz,rz=values[zero];scale=binding['states'][bare]['sign']*mass
        delta=rb-rz;relative=delta/(1+rz);assert relative>-1
        records.extend([
            {'function':'phi','argument':rz,'coefficient':scale*(1/ub**3-1/uz**3)},
            {'function':'log','argument':rz,'coefficient':scale*delta/ub**3},
            {'function':'phi','argument':relative,'coefficient':scale*(1+rz)/ub**3}])
    for index in binding['unpaired']:
        u,r=values[index]
        records.append({'function':'phi','argument':r,'coefficient':binding['states'][index]['sign']*mass/u**3})
    return records


def table(binding,power=2):
    gamma={}
    for mask in range(1,(1<<EDGES)-1):
        estimates=[]
        for state in binding['states']:
            nu=min((mask&m).bit_count() for m,c in state['U'])
            nf=min((mask&m).bit_count() for m,c in state['F'])
            assert nf>=nu
            estimates.append(mask.bit_count()+2*nf-5*nu)
        gamma[mask]=min(estimates)
    return H.build(EDGES,gamma,power)


def verify_table(table):
    n=table['edges'];full=(1<<n)-1;assert n==EDGES
    assert set(table['weights'])==set(range(full)) and table['weights'][0]==table['scale']
    assert set(table['shapes'])==set(range(1,full)) and min(table['weights'].values())>0
    probabilities={full:Q(1)};counts={full:1};transitions=0
    for size in range(n,0,-1):
        for mask in range(1,full+1):
            if mask.bit_count()!=size:continue
            options=[mask^(1<<e) for e in range(n) if mask&(1<<e)]
            total=sum(table['weights'][s] for s in options)
            if mask==full:assert total==table['normalization']
            else:
                shape=table['shapes'][mask];beta=Q(shape) if shape>0 else Q(1,-shape)
                assert shape>=1 or shape in (-2,-4)
                assert table['weights'][mask]*beta==total
            assert total<2**62
            for child in options:
                probabilities[child]=probabilities.get(child,Q())+probabilities[mask]*Q(table['weights'][child],total)
                counts[child]=counts.get(child,0)+counts[mask];transitions+=1
    assert probabilities[0]==1 and counts[0]==math.factorial(n)
    return {'positive_subset_partitions':full,'transitions':transitions,'complete_order_paths':counts[0],
        'path_mass':probabilities[0],'method':'exact DAG recurrence, not explicit permutation enumeration'}
