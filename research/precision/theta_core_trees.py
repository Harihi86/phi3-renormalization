"""Recognize a K2,3 core decorated by recursively nested E2 cycle children."""
from itertools import combinations

def recognize(desc):
    graph=desc['graph'];edges=graph['edges'];parts=desc['proper_subgraphs'];sets=[frozenset(p['edge_ids']) for p in parts]
    assert desc['external_legs']==3 and all(p['external_legs']==2 and p['omega']==2 for p in parts) and len(sets)==len(set(sets))
    vertices=lambda es:set(v for e in es for v in edges[e])
    for a,b in combinations(sets,2):assert a<b or b<a or not(vertices(a)&vertices(b))
    def build(node,is_root=False):
        vs=vertices(node);inside=[s for s in sets if s<node];maximum=[s for s in inside if not any(s<t for t in inside)];children=[build(s) for s in maximum]
        ports=[v for e,(a,b) in enumerate(edges) if e not in node for v in (a,b) if v in vs]+[v for label,v in graph['external'] if v in vs]
        assert len(ports)==(3 if is_root else 2) and len(set(ports))==len(ports)
        reps={v:v for v in vs};child_reps=[]
        for child in maximum:
            cv=vertices(child);rep=min(cv);child_reps.append(rep)
            for v in cv:reps[v]=rep
        assert all(v not in vertices(child) for v in ports for child in maximum)
        remaining=set(node)-set().union(*maximum);adj={reps[v]:[] for v in vs}
        for e in sorted(remaining):
            a,b=map(reps.get,edges[e]);assert a!=b;adj[a].append((e,b));adj[b].append((e,a))
        loops=len(node)-len(vs)+1
        common={'edge_ids':sorted(node),'loops':loops,'external_legs':len(ports),'external_vertices':ports,'children':children,'source_part_index':sets.index(node) if node in sets else None}
        if not is_root:
            assert all(len(v)==2 for v in adj.values());start=min(ports);nodes=[start];cycle=[];edge,v=adj[start][0]
            while True:
                cycle.append(edge)
                if v==start:break
                assert v not in nodes;nodes.append(v);edge,v=next(pair for pair in adj[v] if pair[0]!=edge)
            assert set(cycle)==remaining and set(nodes)==set(adj);cuts=sorted(nodes.index(v) for v in ports);arcs=[]
            for i,a in enumerate(cuts):
                b=cuts[(i+1)%len(cuts)]
                if b<=a:b+=len(nodes)
                arcs.append({'edges':[cycle[j%len(cycle)] for j in range(a,b)],'child_representatives':[nodes[j%len(nodes)] for j in range(a+1,b)]})
            assert sorted(v for a in arcs for v in a['child_representatives'])==sorted(child_reps) and loops==1+sum(c['loops'] for c in children)
            return {**common,'kind':'cycle','arcs':arcs}
        hubs=sorted(v for v in adj if len(adj[v])==3);assert len(hubs)==2 and all(len(adj[v])==2 for v in adj if v not in hubs) and not set(hubs)&set(ports)
        paths=[];visited=[]
        for first,next_vertex in adj[hubs[0]]:
            nodes=[hubs[0]];path=[];edge,v=first,next_vertex
            while True:
                path.append(edge);nodes.append(v)
                if v==hubs[1]:break
                assert len(adj[v])==2 and nodes.count(v)==1;edge,v=next(pair for pair in adj[v] if pair[0]!=edge)
            external=[v for v in nodes if v in ports];assert len(external)==1;cut=nodes.index(external[0]);arcs=[{'edges':path[:cut],'child_representatives':nodes[1:cut]},{'edges':path[cut:],'child_representatives':nodes[cut+1:-1]}]
            paths.append({'port':external[0],'arcs':arcs});visited+=path
        assert sorted(visited)==sorted(remaining) and {p['port'] for p in paths}==set(ports)
        assert sorted(v for p in paths for a in p['arcs'] for v in a['child_representatives'])==sorted(child_reps)
        assert loops==2+sum(c['loops'] for c in children)
        return {**common,'kind':'theta_K2_3','hubs':hubs,'paths':paths}
    tree=build(frozenset(range(len(edges))),True);words={frozenset(c) for n in range(len(parts)+1) for c in combinations(range(len(parts)),n)}
    assert len(desc['terms'])==len(words) and {frozenset(t['subgraphs_inner_first']) for t in desc['terms']}==words
    assert all(t['sign']==(-1)**len(t['subgraphs_inner_first']) for t in desc['terms'])
    return tree
