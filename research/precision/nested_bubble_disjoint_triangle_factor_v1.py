"""Exact eight-forest reduction of a nested bubble vertex and disjoint triangle.

All bases are discovered from original edge routes. Local Taylor projections
are compared before eliminating loop variables. This proves an integrand
identity; it does not certify a numerical value of the remaining integral.
"""
from itertools import combinations, product
import sympy as S
import primitive_triangle_factor_v3 as P


def convert(rows, choices):
    L = len(choices)
    change = S.Matrix.vstack(*(rows[e][:, :L] for e in choices))
    det = S.factor(change.det())
    if det not in (-1, 1):
        return None
    offset = S.Matrix.vstack(*(rows[e][:, L:] for e in choices))
    inverse = change.inv()
    routes = []
    for row in rows:
        loop = row[:, :L] * inverse
        external = row[:, L:] - loop * offset
        routes.append(S.Matrix([[S.expand(v) for v in list(loop) + list(external)]]))
    return dict(change=change, offset=offset, det=det, routes=routes)


def recognize(gd, z):
    parts = gd['proper_subgraphs']
    assert len(parts) == 3
    def find(kind):
        hits = [i for i, p in enumerate(parts)
                if (p['loops'], p['external_legs'], p['omega'], len(p['edge_ids'])) == kind]
        assert len(hits) == 1
        return hits[0]
    bub, tri, par = find((1, 2, 2, 2)), find((1, 3, 0, 3)), find((2, 3, 0, 6))
    be, te, pe = [parts[i]['edge_ids'] for i in (bub, tri, par)]
    assert set(be) < set(pe) and not set(pe) & set(te)
    inner = sorted(set(pe) - set(be))
    host = sorted(set(range(gd['edges'])) - set(pe) - set(te))
    assert len(inner) == 4 and len(host) == 3
    source = [P.scaled_routes(t, z) for t in gd['terms']]
    for choices in product(be, inner, te, host):
        records = []
        for term, rows in zip(gd['terms'], source):
            w = convert(rows, choices)
            if w is None:
                break
            r = w['routes']; word = term['subgraphs_inner_first']
            if any(r[e][0] or r[e][1] or r[e][2] or r[e][3] not in (-1, 1)
                   or any(v.has(*z) for v in r[e]) for e in host):
                break
            if any(r[e][0] or r[e][1] or r[e][2] not in (-1, 1) for e in te):
                break
            if any(r[e][0] or r[e][2] or r[e][1] not in (-1, 1) for e in inner):
                break
            if any(r[e][0] not in (-1, 1) or r[e][2] for e in be):
                break
            shifts = [r[e][:, 1:] / r[e][0] for e in be]
            delta = shifts[0] - shifts[1]
            b = z[bub] if bub in word else S.Integer(1)
            if delta[0] not in (-b, b) or any(delta[j] for j in range(1, 5)):
                break
            normalized = [r[e][:, 3:] / r[e][1] for e in inner]
            if par in word:
                if any(any(v for v in row) for row in normalized):
                    break
            elif sum(row == S.zeros(1, 3) for row in normalized) != 2:
                break
            records.append(w)
        if len(records) == len(source):
            return bub, tri, par, be, te, pe, inner, host, choices, records
    raise AssertionError('No common edge basis isolates the nested and disjoint loops')


def verify(raw):
    gd, pd = raw['gaussian_full_descriptor'], raw['positive_full_descriptor']
    assert (gd['dimension'], gd['loops'], gd['external_legs'], gd['overall_degree'], gd['edges']) == (6, 4, 3, 0, 12)
    assert gd['normalization'] == 'd^6k/pi^3 per loop; no symmetry weight or coupling'
    assert list(map(S.Rational, raw['positive_invariants'])) == [1, 1, 1]
    gram = S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram == S.Matrix([[1, -S.Rational(1, 2)], [-S.Rational(1, 2), 1]])
    parts = gd['proper_subgraphs']; z = S.symbols('z:3'); alpha = S.symbols('a:12')
    assert len(parts) == 3 and len(gd['terms']) == len(pd['states']) == 8
    bub, tri, par, be, te, pe, inner, host, choices, witnesses = recognize(gd, z)
    wanted = {frozenset(i for i in range(3) if mask & (1 << i)) for mask in range(8)}
    assert {frozenset(t['subgraphs_inner_first']) for t in gd['terms']} == wanted
    assert {frozenset(t['word']) for t in pd['states']} == wanted
    lb, lt = sum(alpha[e] for e in be), sum(alpha[e] for e in te)
    combined = lb * lt; kept = [1, 3, 4, 5]; removed = [0, 2]
    records = []; matrices = {}
    for term, w in zip(gd['terms'], witnesses):
        word = tuple(term['subgraphs_inner_first']); key = frozenset(word)
        state = next(s for s in pd['states'] if frozenset(s['word']) == key)
        assert term['sign'] == state['sign'] == (-1) ** len(word)
        assert bub not in word or par not in word or word.index(bub) < word.index(par)
        jets = [bub] if bub in word else []
        assert term['jet_subgraphs'] == state['variables'] == jets
        assert term['jet_orders'] == ([2] if jets else [])
        r = w['routes']; M = sum((a * row.T * row for a, row in zip(alpha, r)), S.zeros(6))
        assert M.extract(removed, removed) == S.diag(lb, lt)
        reduced_rows = [row[:, kept] for row in r]
        H = sum((alpha[e] * reduced_rows[e].T * reduced_rows[e] for e in host), S.zeros(4))
        N = sum((alpha[e] * reduced_rows[e].T * reduced_rows[e] for e in inner), S.zeros(4))
        bs = [reduced_rows[e] / r[e][0] for e in be]
        bv = (alpha[be[0]] * alpha[be[1]] * (bs[0] - bs[1]).T * (bs[0] - bs[1])).applyfunc(S.expand)
        b = z[bub] if bub in word else S.Integer(1)
        expected_bv = S.zeros(4); expected_bv[0, 0] = alpha[be[0]] * alpha[be[1]] * b * b
        assert bv == expected_bv, 'Bubble external momentum must be the retained inner loop k'
        ts = [reduced_rows[e] / r[e][2] for e in te]
        tv = sum((alpha[te[i]] * alpha[te[j]] * (ts[i] - ts[j]).T * (ts[i] - ts[j])
                  for i, j in combinations(range(3), 2)), S.zeros(4)).applyfunc(S.expand)
        assert all(tv[0, j] == tv[j, 0] == 0 for j in range(4))
        if tri in word:
            assert tv == S.zeros(4)
        reduced = (combined * (H + N) + lt * bv + lb * tv).applyfunc(S.expand)
        direct = combined * M.extract(kept, kept) - lt * M.extract(kept, [0]) * M.extract([0], kept) - lb * M.extract(kept, [2]) * M.extract([2], kept)
        assert (direct - reduced).applyfunc(S.expand) == S.zeros(4)
        # Fraction-free Schur construction is independent of the positive U/F list.
        ru = S.Poly(reduced[:2, :2].det(method='domain-ge'), *alpha, *z).exquo(S.Poly(combined, *alpha, *z)).as_expr()
        rf0 = sum(gram[i, j] * reduced.extract([0, 1, 2 + i], [0, 1, 2 + j]).det(method='domain-ge') for i in range(2) for j in range(2))
        rf = S.Poly(rf0, *alpha, *z).exquo(S.Poly(combined ** 2, *alpha, *z)).as_expr()
        original = sum((a * row.T * row for a, row in zip(alpha, P.scaled_routes(term, z))), S.zeros(6))
        u = original[:4, :4].det(method='domain-ge')
        f = sum(gram[i, j] * original.extract([0, 1, 2, 3, 4 + i], [0, 1, 2, 3, 4 + j]).det(method='domain-ge') for i in range(2) for j in range(2))
        counts = []
        for name, old, new in [('U', u, ru), ('F', f, rf)]:
            positive = P.positive_polynomial(state[name], state, alpha, z)
            assert S.expand(old - positive) == 0 and S.expand(new - positive) == 0, (word, name)
            counts.append(len(S.Poly(positive, *alpha, *z).terms()))
        matrices[key] = (H, N, bv, tv, r)
        records.append({'word':list(word), 'sign':term['sign'], 'jet_orders':term['jet_orders'],
                        'basis_original_edges':list(choices), 'jacobian_determinant':int(w['det']),
                        'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],
                        'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],
                        'routing':[[str(v) for v in row] for row in r], 'U_F_monomials':counts,
                        'three_exact_U_F_constructions_equal':True})
    H0, N0, bv0, tv0, bare = matrices[frozenset()]
    parent_zero = matrices[frozenset([par])][1]
    for word, (H, N, bv, tv, r) in matrices.items():
        assert H == H0
        assert N == (parent_zero if par in word else N0)
        assert bv == (bv0 * z[bub] ** 2 if bub in word else bv0)
        assert tv == (S.zeros(4) if tri in word else tv0)
        assert all(r[e] == bare[e] for e in host)
    shifted = [(e, bare[e][:, 3:] / bare[e][1]) for e in inner]
    shifted = [(e, row) for e, row in shifted if row != S.zeros(1, 3)]
    assert len(shifted) == 2
    tri_shift = [bare[e][:, 3:] / bare[e][2] for e in te]
    return {'id':raw['id'], 'bubble_index':bub, 'triangle_index':tri, 'parent_index':par,
            'bubble_edges':be, 'triangle_edges':te, 'parent_edges':pe, 'host_edges':host,
            'states':records, 'parent_shift_edges':[e for e, _ in shifted],
            'parent_boundary_shifts':[[int(v) for v in row] for _, row in shifted],
            'triangle_pair_differences':[[int(v) for v in tri_shift[i]-tri_shift[j]] for i, j in combinations(range(3), 2)],
            'host_routes':[[int(v) for v in bare[e][:, 3:]] for e in host],
            'bubble_boundary':'retained inner vertex loop k, not outer q or physical external p',
            'local_parent_T0':'zero parent boundary shifts, retain k and every host propagator',
            'overall_T0':'zero physical p1,p2 at fixed q; keep q in both reduced vertex kernels',
            'complete_reduced_integrand':'(1-T_physical^0)[H(q,p) * integral_simplex K(x,V_parent(q,p)) * integral_simplex -log(1+V_triangle(q,p))]',
            'kernel_K':'-integral_(1-x)/x^((1-x+V)/x) B(t)/(1+t)^2 dt',
            'simplex_measures':['(1-x) dx dy', '(1-u) du dv'],
            'no_symmetry_weight_or_coupling':True, 'all_eight_forests_bound':True,
            'three_child_loops_integrated_analytically':True, 'whole_integral_evaluated':False}
