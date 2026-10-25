"""Shared exact graph identities and explicit half-edge isomorphisms (4.3 CT)."""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from .graphs import Graph, canonical, _encoding


@lru_cache(maxsize=32768)
def identity(g, fixed=True):
    key, cg = canonical(g, fixed)
    return 'g_' + sha256(json.dumps(key, separators=(',', ':')).encode()).hexdigest()[:20]


def isomorphism(g, target):
    """Return vertex, edge and half-edge maps; fixed external labels never move."""
    if g.n != target.n or len(g.edges) != len(target.edges):
        raise ValueError('Different sizes')
    wanted = _encoding(target, tuple(range(target.n)), True)
    for p in permutations(range(g.n)):
        if _encoding(g, p, True) != wanted:
            continue
        bins = defaultdict(list)
        for i, e in enumerate(target.edges):
            bins[e].append(i)
        edge_map = {}; half_map = {}
        for i, (a, b) in enumerate(g.edges):
            e = tuple(sorted((p[a], p[b])))
            j = bins[e].pop(0); edge_map[str(i)] = j
            flip = int(p[a] > p[b])
            half_map[f'{i}:0'] = [j, flip]
            half_map[f'{i}:1'] = [j, 1-flip]
        return {'vertices': list(p), 'edges': edge_map, 'half_edges': half_map,
                'external': {str(k): [k, p[v]] for k, v in g.external}}
    raise ValueError('Not isomorphic with fixed labels and vertex types')


def canonical_record(g):
    _, cg = canonical(g)
    return {'id': identity(cg), 'graph': cg.to_dict(), 'mapping': isomorphism(g, cg)}


def marked_identity(g, edge_ids):
    """Isomorphism orbit of an embedded edge set, including parallel edges."""
    chosen = set(edge_ids); best = None
    for p in permutations(range(g.n)):
        colored_edges = Counter((min(p[a],p[b]), max(p[a],p[b]), int(i in chosen))
                                for i,(a,b) in enumerate(g.edges))
        key = (_encoding(g,p,True)[:2], tuple(sorted(colored_edges.items())))
        if best is None or key < best: best = key
    return 'm_' + sha256(json.dumps(best, separators=(',', ':')).encode()).hexdigest()[:20]


def external_orbit(g):
    return sorted({identity(Graph(g.n,g.edges,tuple((p[k-1],v) for k,v in g.external),g.kinds))
                   for p in permutations(range(1,len(g.external)+1))})


def momentum_transport(g, target, mapping):
    """Exact change of loop basis induced by an oriented edge isomorphism."""
    from .routing import momentum_routing
    before=momentum_routing(g); after=momentum_routing(target)
    width=len(before['basis']); change=[]
    for i in before['chord_edges']:
        edge=mapping['edges'][str(i)]
        sign=1 if mapping['half_edges'][f'{i}:0'][1]==0 else -1
        change.append([sign*c for c in after['edge_momenta'][str(edge)]])
    for j in range(g.loops,width):change.append([int(i==j) for i in range(width)])
    checks=[]
    for i,row in before['edge_momenta'].items():
        edge=mapping['edges'][i]
        sign=1 if mapping['half_edges'][f'{i}:0'][1]==0 else -1
        transformed=[sum(row[j]*change[j][k] for j in range(width)) for k in range(width)]
        checks.append(transformed==[sign*c for c in after['edge_momenta'][str(edge)]])
    if not all(checks):raise ArithmeticError('Loop change-of-variables fails on an edge')
    return {'from_basis':before['basis'],'to_basis':after['basis'],'matrix':change,
            'equation':'from_basis = matrix * to_basis','all_edge_momenta_verified':all(checks)}
