"""Cycle skeletons with E2 children on any number of branches."""
from itertools import combinations
def recognize(desc):
    graph=desc['graph'];edges=graph['edges'];parts=desc['proper_subgraphs']
    assert desc['external_legs'] in (2,3) and all(s['external_legs']==2 and s['omega']==2 for s in parts)
    sets=[frozenset(s['edge_ids']) for s in parts]
    assert len(set(sets))==len(sets)
    vertices=lambda es:set(v for e in es for v in edges[e])
    for a,b in combinations(sets,2):assert a<b or b<a or not vertices(a)&vertices(b)
    def build(node):
        vs=vertices(node);inside=[s for s in sets if s<node]
        maximum=[s for s in inside if not any(s<t for t in inside)]
        child_nodes=[build(s) for s in maximum];port_vertices=[]
        for e,(u,v) in enumerate(edges):
            if e not in node:
                if u in vs:port_vertices.append(u)
                if v in vs:port_vertices.append(v)
        port_vertices.extend(v for label,v in graph['external'] if v in vs)
        assert len(port_vertices) in (2,3) and len(set(port_vertices))==len(port_vertices)
        if len(node)!=len(edges):assert len(port_vertices)==2
        reps={v:v for v in vs};child_reps=[]
        for child in maximum:
            cv=vertices(child);rep=min(cv);child_reps.append(rep)
            for v in cv:reps[v]=rep
        assert all(v not in vertices(child) for v in port_vertices for child in maximum)
        remaining=set(node)-set().union(*maximum);adj={reps[v]:[] for v in vs}
        for e in sorted(remaining):
            u,v=map(reps.get,edges[e]);assert u!=v
            adj[u].append((e,v));adj[v].append((e,u))
        assert all(len(v)==2 for v in adj.values())
        start=min(port_vertices);nodes=[start];cycle=[];edge,v=adj[start][0]
        while True:
            cycle.append(edge)
            if v==start:break
            assert v not in nodes;nodes.append(v)
            edge,v=next(pair for pair in adj[v] if pair[0]!=edge)
        assert set(cycle)==remaining and set(nodes)==set(adj)
        indices=sorted(nodes.index(v) for v in port_vertices);arcs=[]
        for i,a in enumerate(indices):
            b=indices[(i+1)%len(indices)]
            if b<=a:b+=len(nodes)
            arc=[cycle[j%len(cycle)] for j in range(a,b)]
            internal=[nodes[j%len(nodes)] for j in range(a+1,b)]
            arcs.append({'edges':arc,'child_representatives':internal})
        decorated_arcs=sum(len(a['edges'])>1 for a in arcs)
        assert sorted(v for a in arcs for v in a['child_representatives'])==sorted(child_reps)
        loops=len(node)-len(vs)+1
        assert loops==1+sum(c['loops'] for c in child_nodes)
        return {'edge_ids':sorted(node),'loops':loops,'external_legs':len(port_vertices),'external_vertices':port_vertices,
            'arcs':arcs,'decorated_arcs':decorated_arcs,'children':child_nodes,'source_part_index':sets.index(node) if node in sets else None}
    tree=build(frozenset(range(len(edges))))
    all_subsets={frozenset(c) for n in range(len(parts)+1) for c in combinations(range(len(parts)),n)}
    assert len(desc['terms'])==len(all_subsets)
    assert {frozenset(t['subgraphs_inner_first']) for t in desc['terms']}==all_subsets
    assert all(t['sign']==(-1)**len(t['subgraphs_inner_first']) for t in desc['terms'])
    return tree
