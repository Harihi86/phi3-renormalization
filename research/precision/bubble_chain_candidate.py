"""Positive radial candidate for graph-certified disjoint bubble chains.

No numerical result from this module promotes an original campaign case.
The normalization is one d^6 k / pi^3 per loop, without symmetry/coupling.
"""
from fractions import Fraction as Q
from itertools import combinations, product
import sympy as sp


def recognize(descriptor):
    graph = descriptor['graph']; edges = graph['edges']; n = graph['vertices']
    assert descriptor['dimension'] == 6 and descriptor['external_legs'] == 2
    assert descriptor['overall_degree'] == 2
    proper = descriptor['proper_subgraphs']
    assert proper and all(s['loops'] == 1 and s['external_legs'] == 2 and s['omega'] == 2 and len(s['edge_ids']) == 2 for s in proper)
    parts = [set(s['vertices']) for s in proper]
    assert all(len(v) == 2 for v in parts)
    assert all(not a & b for a, b in combinations(parts, 2))
    for s, vertices in zip(proper, parts):
        assert all(set(edges[e]) == vertices for e in s['edge_ids'])
    k = len(parts); excluded = {e for s in proper for e in s['edge_ids']}
    rep = {i: i for i in range(n)}
    for vs in parts:
        for v in vs: rep[v] = min(vs)
    adj = {rep[v]: [] for v in range(n)}
    for i, (u, v) in enumerate(edges):
        if i not in excluded:
            u, v = rep[u], rep[v]; assert u != v
            adj[u].append((i, v)); adj[v].append((i, u))
    assert all(len(v) == 2 for v in adj.values())
    external = graph['external']; assert len(external) == 2
    a, b = [rep[e[1]] for e in external]; assert a != b
    assert a not in {min(v) for v in parts} and b not in {min(v) for v in parts}
    paths = []
    for first in adj[a]:
        path = []; nodes = [a]; edge, v = first
        while True:
            assert v not in nodes
            path.append(edge); nodes.append(v)
            if v == b: break
            edge, v = next(x for x in adj[v] if x[0] != edge)
        paths.append((path, nodes))
    paths.sort(key=lambda p: len(p[0]))
    assert len(paths[0][0]) == 1 and len(paths[1][0]) == k + 1
    assert set(paths[0][0] + paths[1][0]) == set(range(len(edges))) - excluded
    assert set(paths[1][1][1:-1]) == {min(v) for v in parts}
    words = [tuple(t['subgraphs_inner_first']) for t in descriptor['terms']]
    expected = {tuple(c) for r in range(k + 1) for c in combinations(range(k), r)}
    assert set(words) == expected and len(words) == len(expected)
    witnesses = []
    for term in descriptor['terms']:
        word = term['subgraphs_inner_first']; assert term['sign'] == (-1)**len(word)
        raw = term['unscaled_routing']; scaled = term['scaled_routing']
        assert len(raw[0]) == k + 2
        qrow = raw[paths[1][0][0]]
        found = None
        for local_edges in product(*(s['edge_ids'] for s in proper)):
            selected = [raw[e] for e in local_edges] + [qrow]
            matrix = sp.Matrix([r[:-1] for r in selected])
            det = matrix.det()
            if det not in (-1, 1): continue
            inverse = matrix.inv(); offsets = sp.Matrix([r[-1] for r in selected])
            transformed = []
            for row in scaled:
                powers = {tuple(v['powers']) for v in row if v['coefficient']}
                coefficients = {}
                for power in powers:
                    values = [Q(v['coefficient']) if tuple(v['powers']) == power else Q() for v in row]
                    loop = sp.Matrix([values[:-1]]) * inverse
                    ext = sp.Rational(values[-1]) - (loop * offsets)[0]
                    for j, c in enumerate(list(loop) + [ext]):
                        if c: coefficients[(j, power)] = Q(c)
                transformed.append(coefficients)
            zero = (0,)*len(word)
            ok = True
            for j, s in enumerate(proper):
                qpower = list(zero)
                if j in word: qpower[word.index(j)] = 1
                allowed = {(j, zero), (k, tuple(qpower))}
                rows = [transformed[e] for e in s['edge_ids']]
                ok &= all(set(r) <= allowed and abs(r.get((j, zero), 0)) == 1 for r in rows)
                ok &= sorted(abs(r.get((k, tuple(qpower)), 0)) for r in rows) == [0, 1]
            for e in paths[1][0]:
                ok &= set(transformed[e]) == {(k, zero)} and abs(transformed[e].get((k, zero), 0)) == 1
            outer = transformed[paths[0][0][0]]
            ok &= set(outer) == {(k, zero), (k+1, zero)} and all(abs(c) == 1 for c in outer.values())
            if ok:
                found = {'loop_edges': list(local_edges), 'chain_momentum_edge': paths[1][0][0],
                         'jacobian_determinant': int(det), 'loop_affine_rows': selected,
                         'transformed_rows': [[{'index': j, 'powers': list(p), 'coefficient': str(c)} for (j,p),c in sorted(r.items())] for r in transformed]}
                break
        assert found is not None, ('No exact local Taylor routing factorization', word)
        witnesses.append({'word': word, 'sign': term['sign'], **found})
    return {'bubble_count': k, 'loops': k+1, 'bubble_edges': [s['edge_ids'] for s in proper],
            'chain_edges': paths[1][0], 'outer_edge': paths[0][0][0], 'states': witnesses,
            'all_original_recursive_states_retained': True, 'all_local_Taylor_orders': 2,
            'normalization': 'd^6k/pi^3; no symmetry or coupling',
            'radial_formula': 'integral_0^infinity t^2/2 B_R(t)^k/(t+M)^(k+1) R_angle(t,P,M) dt',
            'main162_new_passes': 0}


def bubble(ctx, t, mass2=1):
    t, m = ctx.mpf(t), ctx.mpf(mass2)
    if t == 0: return ctx.mpf(0)
    z = t/(t + 4*m)
    if z < ctx.mpf('0.2'):
        value = ctx.mpf(0); power = z*z; j = 2
        while True:
            term = power/(2*j+1); value += term
            if term < ctx.eps*value: break
            j += 1; power *= z
        return (t + 4*m)*value/3
    r = ctx.sqrt(z)
    return (t+4*m)/3*(ctx.asinh(ctx.sqrt(t/(4*m)))/r-1)-t/9


def angular(ctx, t, p2=1, mass2=1):
    a = t + p2 + mass2; b2 = 4*t*p2
    h = a + ctx.sqrt(a*a-b2)
    return 2/h-2*b2/(3*h**3)


def angular_remainder(ctx, t, p2=1, mass2=1):
    t,p,m = map(ctx.mpf, (t,p2,mass2)); d = t+m
    if p == 0: return ctx.mpf(0)
    root = ctx.sqrt((d-p)**2+4*m*p)
    delta = 2*m*p/(d-p+root) if d >= p else (p-d+root)/2
    e = d+delta
    return p*delta*(3*d*d*m+6*d*m*delta+(3*m+t)*delta*delta)/(3*d**3*e**3)


def radial(ctx, t, k=2, p2=1, mass2=1):
    return t*t/2*bubble(ctx,t,mass2)**k/(t+mass2)**(k+1)*angular_remainder(ctx,t,p2,mass2)
