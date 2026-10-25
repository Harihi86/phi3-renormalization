"""Independent exhaustive original E2 graph subsets; no imported list of divergences."""
from itertools import combinations
from collections import Counter

def subgraph_audit(raw):
    graph=raw['graph'];edges=graph['edges'];n=graph['vertices']
    assert n==8 and len(edges)==11 and graph['vertex_kinds']==['phi3']*n
    degree=Counter(v for e in edges for v in e)
    for label,v in graph['external']:degree[v]+=1
    assert len(graph['external'])==2 and all(degree[v]==3 for v in range(n))
    def connected(ids,vertices):
        reached={next(iter(vertices))}
        while True:
            old=len(reached)
            for i in ids:
                a,b=edges[i]
                if a in reached or b in reached:reached.update((a,b))
            if len(reached)==old:return reached==vertices
    divergent=[];onepi=0
    for size in range(1,12):
        for ids in combinations(range(11),size):
            vertices={v for i in ids for v in edges[i]};loops=size-len(vertices)+1
            if loops<=0 or not connected(ids,vertices):continue
            if any(not connected([i for i in ids if i!=cut],vertices) for cut in ids):continue
            onepi+=1;external=3*len(vertices)-2*size;omega=6*loops-2*size
            assert omega==6-2*external
            if omega>=0:divergent.append((tuple(ids),loops,external,omega))
    actual={(tuple(sorted(p['edge_ids'])),p['loops'],p['external_legs'],p['omega']) for p in raw['gaussian_full_descriptor']['proper_subgraphs']}
    expected=set(divergent)-{(tuple(range(11)),4,2,2)}
    assert actual==expected and len(actual) in (4,5) and (tuple(range(11)),4,2,2) in divergent
    assert raw['gaussian_full_descriptor']['overall_degree']==2
    return {'edge_subsets_examined':2047,'one_particle_irreducible_loop_subgraphs':onepi,
            'all_nonnegative_degree_subgraphs':divergent,'all_proper_parts_and_orders_match':True}

