"""Use alternative exact forest matchings to bound each face of one fixed density.

All subgraph priority permutations of a deterministic augmenting-path matcher
are tried, with two left-vertex traversals. This does not enumerate or claim
all possible matchings. Each resulting representation retains every forest.
"""
from itertools import permutations
import numpy as np
import projective_pair_face_transport_v1 as A


def matchings(states):
    words=[frozenset(s['word']) for s in states]
    parts=sorted(set().union(*words));left=[i for i,s in enumerate(states) if s['sign']==1]
    right=[i for i,s in enumerate(states) if s['sign']==-1]
    assert len(parts)<=6
    edges={i:[j for j in right if len(words[i]^words[j])==1] for i in left};found=set()
    for priority in permutations(parts):
        rank={part:i for i,part in enumerate(priority)}
        neighbors={i:sorted(edges[i],key=lambda j:(rank[next(iter(words[i]^words[j]))],j)) for i in left}
        for traversal in (left,list(reversed(left))):
            owner={}
            def assign(i,seen):
                for j in neighbors[i]:
                    if j in seen:continue
                    seen.add(j)
                    if j not in owner or assign(owner[j],seen):owner[j]=i;return True
                return False
            for i in traversal:assign(i,set())
            pairs=[];used=set()
            for j,i in owner.items():
                b,z=(i,j) if words[i]<words[j] else (j,i)
                assert words[b]<words[z] and len(words[z]-words[b])==1
                pairs.append((b,z));used.update((b,z))
            unpaired=tuple(i for i in range(len(states)) if i not in used)
            found.add((tuple(sorted(pairs)),unpaired))
    return [{'pairs':[list(p) for p in pairs],'unpaired':list(unpaired)} for pairs,unpaired in sorted(found)]


def construct(binding,two_point):
    base=A.construct(binding,two_point);n=binding['original_edges'];full=(1<<n)-1
    states=binding['states'];options=matchings(states)
    original={'pairs':binding['pairs'],'unpaired':binding.get('unpaired',[])}
    if original not in options:options.append(original)
    orders=[];single=[];factor=2 if two_point else 1
    for state in states:
        u=A.all_valuations([[m,0,int(c)] for m,c in state['U']],n)
        f=A.all_valuations([[m,0,int(c)] for m,c in state['F']],n)
        orders.append((u,f));single.append(-3*u+factor*(f-u))
    edges=sorted(set(tuple(pair) for option in options for pair in option['pairs']))
    bounds={};numerators={}
    for b,z in edges:
        ub,fb=orders[b];uz,fz=orders[z];rz=fz-uz
        du=[[m,0,c] for m,c in A.difference(states[b]['U'],states[z]['U'])]
        cross=A.cross_difference(states[b]['F'],states[z]['U'],states[z]['F'],states[b]['U'])
        vu=A.all_valuations(du,n);vc=A.all_valuations(cross,n);terms=[]
        if vu is not None:terms.append(vu+2*np.minimum(ub,uz)-3*ub-3*uz+factor*rz)
        if vc is not None:
            d=vc-ub-uz;assert np.all(d>=0)
            terms.append(-3*ub+d+rz if two_point else -3*ub+d)
            if two_point:terms.append(-3*ub+2*d)
        # An identically zero pair imposes no bound on a sum containing it.
        bounds[(b,z)]=np.min(np.stack(terms),axis=0) if terms else np.full(full,1000,dtype=np.int16)
        numerators[(b,z)]={'pair':[b,z],'U_difference':du,'FU_cross_difference':cross}
    cardinality=np.array([i.bit_count() for i in range(full)],dtype=np.int16)
    prior=np.array([0]+[base['gamma'][i] for i in range(1,full)],dtype=np.int16)
    best=prior.copy();witness=np.full(full,-1,dtype=np.int16);all_bounds=[]
    for index,option in enumerate(options):
        used=[i for pair in option['pairs'] for i in pair]+option['unpaired']
        assert sorted(used)==list(range(len(states)))
        pieces=[bounds[tuple(pair)] for pair in option['pairs']]+[single[i] for i in option['unpaired']]
        value=np.min(np.stack(pieces),axis=0)+cardinality;all_bounds.append(value)
        improved=value>best;witness[improved]=index;best=np.maximum(best,value)
    assert np.all(best[1:]>=prior[1:]) and max(int(x) for x in best[1:])<1000
    gamma={i:int(best[i]) for i in range(1,full)}
    proof={'id':binding['id'],'two_point':two_point,'alternative_matchings':options,
        'witness':{i:int(witness[i]) for i in range(1,full) if witness[i]>=0},
        'gamma':gamma,'fixed_pair_gamma':base['gamma'],'original_gamma':base['old_gamma'],
        'numerators':list(numerators.values()),
        'summary':{'id':binding['id'],'candidate_matchings':len(options),'legal_pair_edges':len(edges),
            'faces_improved_over_fixed_pairing':sum(gamma[i]>base['gamma'][i] for i in gamma),
            'new_critical_to_positive_faces':sum(base['gamma'][i]<=0<gamma[i] for i in gamma),
            'total_critical_to_positive_faces':sum(base['old_gamma'][i]<=0<gamma[i] for i in gamma),
            'remaining_nonpositive_faces':sum(gamma[i]<=0 for i in gamma),
            'minimum_gamma':min(gamma.values()),'density_changed':False,
            'global_optimal_matching_claimed':False,'whole_integral_acceptance':False}}
    return proof


def transport(binding,two_point,power=2):
    result=construct(binding,two_point)
    table=A.H.build(binding['original_edges'],{s:min(g,2) for s,g in result['gamma'].items()},power)
    return table,result
