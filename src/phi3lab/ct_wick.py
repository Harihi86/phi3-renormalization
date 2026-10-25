"""Independent typed action/Wick generator. No bare graph catalogue or subgraphs.

Z[J]=exp(-sum_v lambda_v/r_v! (delta_J)^r_v) exp(J D J/2).
Vertex permutations and source-to-stub attachments are summed as exact orbits;
ALL residual labelled half-edge pairings are visited. Matching ranks and one
witness per raw adjacency permit complete deterministic replay.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import factorial, prod
from .graphs import Graph, connected, components, is_1pi
from .higher_graphs import canonical as fast_canonical
from .ct_graph import canonical_record

SPECS = {
    'phi3': (3,0,1),
    'C1': (1,1,1), 'C2m': (2,1,2), 'C2z': (2,1,2), 'C3': (3,1,3),
    'C1L2': (1,2,3), 'C2mL2': (2,2,4), 'C2zL2': (2,2,4), 'C3L2': (3,2,5),
}


def pairings(stubs):
    if not stubs:
        yield (); return
    a = stubs[0]
    for i in range(1,len(stubs)):
        for tail in pairings(stubs[1:i]+stubs[i+1:]):
            yield ((a,stubs[i]),)+tail


def attachment_orbits(kinds, labels):
    """Restricted first-use order only within a vertex species; labels fixed."""
    pools = {k: tuple(i for i,t in enumerate(kinds) if t==k) for k in set(kinds)}
    arities = [SPECS[k][0] for k in kinds]
    for attachment in product(range(len(kinds)), repeat=len(labels)):
        used=defaultdict(list); valid=True
        for v in attachment:
            k=kinds[v]
            if v not in used[k]:
                if v != pools[k][len(used[k])]: valid=False; break
                used[k].append(v)
        if not valid: continue
        counts=Counter(attachment)
        if any(counts[v]>r for v,r in enumerate(arities)): continue
        vertex_choices=prod(factorial(len(pools[k]))//factorial(len(pools[k])-len(u))
                            for k,u in used.items())
        stub_choices=prod(factorial(r)//factorial(r-counts[v]) for v,r in enumerate(arities))
        yield attachment, counts, vertex_choices, stub_choices


def source_matchings(labels):
    """Partial external-source contractions, including the empty matching."""
    if not labels:
        yield (), (); return
    a=labels[0]
    for pairs, attached in source_matchings(labels[1:]):
        yield pairs, (a,)+attached
    for i,b in enumerate(labels[1:],1):
        for pairs,attached in source_matchings(labels[1:i]+labels[i+1:]):
            yield ((a,b),)+pairs, attached


@lru_cache(maxsize=64)
def enumerate_wick(kinds, labels=(1,2,3), include_source_pairs=False):
    kinds=tuple(kinds); labels=tuple(labels); n=len(kinds)
    if not n: raise ValueError('Use at least one interaction vertex')
    species_factorial=prod(factorial(c) for c in Counter(kinds).values())
    vertex_factorial=prod(factorial(SPECS[k][0]) for k in kinds)
    denominator=species_factorial*vertex_factorial
    classes={}; raw=[]; pairings_visited=0; weighted_pairings=0
    source_cases=source_matchings(labels) if include_source_pairs else [((),labels)]
    for source_pairs, attached_labels in source_cases:
        for attachment, counts, choices, stub_choices in attachment_orbits(kinds, attached_labels):
            remaining=tuple((v,s) for v,k in enumerate(kinds) for s in range(counts[v],SPECS[k][0]))
            if len(remaining)%2: continue
            grouped={}
            for rank,matching in enumerate(pairings(remaining)):
                pairings_visited+=1
                edges=tuple(sorted(tuple(sorted((a[0],b[0]))) for a,b in matching))
                if edges not in grouped: grouped[edges]=[0,matching,[]]
                grouped[edges][0]+=1; grouped[edges][2].append(rank)
            for edges,(count,witness,ranks) in grouped.items():
                g=Graph(n,edges,tuple(zip(attached_labels,attachment)),kinds)
                fk,_=fast_canonical(g)
                key=(tuple(source_pairs),fk)
                cpts=components(range(n),edges)
                attach_vertices=set(attachment)
                vacuum=[sorted(c) for c in cpts if not c&attach_vertices]
                conn=not source_pairs and connected(range(n),edges)
                one_pi=conn and is_1pi(range(n),edges)
                multiplicity=count*choices*stub_choices
                weighted_pairings+=multiplicity
                raw_id=len(raw)
                raw.append({'raw_id':raw_id,'graph':g.to_dict(),'source_pairs':source_pairs,
                            'attachment':list(zip(attached_labels,attachment)),
                            'source_half_edges':[[label,v,j] for v in range(n)
                                                 for j,label in enumerate(l for l,w in zip(attached_labels,attachment) if w==v)],
                            'remaining_half_edges':remaining,'matching_count':count,'matching_ranks':ranks,
                            'matching_witness':witness,'vertex_attachment_choices':choices,
                            'stub_attachment_choices':stub_choices,'wick_multiplicity':multiplicity})
                if key not in classes:
                    rec=canonical_record(g)
                    # Missing external labels are paired directly by the Gaussian.
                    suffix='' if not source_pairs else '_sources_'+'_'.join(f'{a}-{b}' for a,b in source_pairs)
                    classes[key]={'id':rec['id']+suffix,'graph':rec['graph'],'source_pairs':source_pairs,
                        'raw_ids':[],'wick_multiplicity':0,'connected':conn,'one_pi':one_pi,
                        'vacuum_components':len(vacuum),'self_loops':sum(a==b for a,b in edges),
                        'accepted':one_pi,
                        'reason':('accepted_amputated_1PI' if one_pi else 'vacuum_factor_cancels_in_Z_over_Z0' if vacuum
                                  else 'disconnected_removed_by_log_Z' if not conn else 'internal_bridge_1PR'),
                        'species_factorial':species_factorial,'vertex_factorial':vertex_factorial}
                row=classes[key]; row['wick_multiplicity']+=multiplicity; row['raw_ids'].append(raw_id)
                raw[-1]['canonical_id']=row['id']; raw[-1]['mapping']=canonical_record(g)['mapping']
                bins=defaultdict(list)
                for edge_id,edge in enumerate(g.edges):bins[edge].append(edge_id)
                halfedges={}
                for left,right in witness:
                    edge_id=bins[tuple(sorted((left[0],right[0])))].pop(0)
                    orientation=int(left[0]>right[0])
                    halfedges[f'{left[0]}:{left[1]}']=raw[-1]['mapping']['half_edges'][f'{edge_id}:{orientation}']
                    halfedges[f'{right[0]}:{right[1]}']=raw[-1]['mapping']['half_edges'][f'{edge_id}:{1-orientation}']
                raw[-1]['witness_stub_to_canonical_halfedge']=halfedges
    records=[]
    for row in classes.values():
        row['positive_weight']=str(Fraction(row['wick_multiplicity'],denominator))
        records.append(row)
    return {'kinds':kinds,'labels':labels,'include_source_pairs':include_source_pairs,
            'denominator':denominator,'species_factorial':species_factorial,
            'vertex_factorial':vertex_factorial,'pairings_visited':pairings_visited,
            'full_wick_multiplicity':weighted_pairings,
            'records':sorted(records,key=lambda r:r['id']),'raw_terms':raw,
            'source_factorial':'External derivatives cancel the source 1/E!; labels remain fixed.'}


def action_sectors(loops=2, legs=3):
    """Solve the action's coupling order and loop grading, before graph generation."""
    power=2*loops+legs-2
    names=tuple(SPECS)
    out=[]
    def visit(i, remaining, kinds):
        if i==len(names):
            if remaining or not kinds:return
            r=sum(SPECS[k][0] for k in kinds); v=len(kinds)
            if r<legs or (r-legs)%2:return
            internal=(r-legs)//2; topological=internal-v+1
            ct_loops=sum(SPECS[k][1] for k in kinds)
            if topological<0 or topological+ct_loops!=loops:return
            out.append({'kinds':tuple(kinds),'coupling_power':power,'topological_loops':topological,
                        'counterterm_loops':ct_loops,'total_loops':loops})
            return
        name=names[i]; degree=SPECS[name][2]
        for number in range(remaining//degree+1):
            visit(i+1,remaining-number*degree,kinds+[name]*number)
    visit(0,power,[])
    return out
