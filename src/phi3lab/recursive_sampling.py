"""Program-derived importance sampling for recursive massive BPH integrals.

Only the probability distribution is constructed here.  Every density must
still be evaluated as the *complete* graph-compiled recursive sum.  Positive
clipped Hepp gaps define normalized proposals, not remainder certificates or
a proof of finite variance.  Spanning trees/two-trees are the graph polynomial
theorem; they are not the BPHZ forest formula.
"""
from functools import lru_cache
from itertools import combinations
from math import isfinite, lgamma, log


def _compiled(value):
    if hasattr(value, 'terms') and hasattr(value, 'graph'):
        return value
    from .automatic_bph import compile_bph
    return compile_bph(value)


def _program_key(compiled):
    parts = tuple((sum(1 << e for e in part['edge_ids']), part['loops'],
                   part['omega']) for part in compiled.parts)
    words = tuple(tuple(term['subgraphs_inner_first']) for term in compiled.terms)
    return compiled.graph, parts, words


@lru_cache(maxsize=256)
def graph_polynomial_support(graph):
    """Original-edge U/F masks, with no routing, jet or numerical kernel calls.

    U consists of complements of spanning trees. F consists of complements of
    spanning two-trees separating physical external legs; its coefficients
    are nonnegative Euclidean invariant squares. Supports are independent of
    momentum values. At exceptional kinematics some F coefficients can vanish,
    which only makes this proposal envelope more conservative.
    """
    size = len(graph.edges)
    if size > 18:
        raise ValueError('Exact recursive importance supports at most18 edges')
    full = (1 << size)-1

    def forest_components(removed, count):
        parent = list(range(graph.n))
        def root(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v
        for edge, (a, b) in enumerate(graph.edges):
            if removed & (1 << edge):
                continue
            a, b = root(a), root(b)
            if a == b:
                return None
            parent[a] = b
        labels = tuple(root(v) for v in range(graph.n))
        return labels if len(set(labels)) == count else None

    u_masks, f_masks = [], []
    for ids in combinations(range(size), graph.loops):
        mask = sum(1 << e for e in ids)
        if forest_components(mask, 1) is not None:
            u_masks.append(mask)
    if len(graph.external) > 1:
        for ids in combinations(range(size), graph.loops+1):
            mask = sum(1 << e for e in ids)
            labels = forest_components(mask, 2)
            if labels is None:
                continue
            first = labels[0]
            number = sum(labels[v] == first for _, v in graph.external)
            if 0 < number < len(graph.external):
                f_masks.append(mask)
    if not u_masks:
        raise ArithmeticError('Connected graph has no spanning tree support')
    return tuple(u_masks), tuple(f_masks)


def _word_supports(parts, word, masks):
    quadratic = [j for j in word if parts[j][2] == 2]
    if any(parts[j][2] not in (0, 2) for j in word):
        raise ValueError('Recursive sampler currently supports quadratic/logarithmic proper parts')
    positions = {j: k for k, j in enumerate(quadratic)}
    supports = [[] for _ in range(1 << len(quadratic))]
    for mask in masks:
        jet = 0
        for j in word:
            edges, loops, degree = parts[j]
            power = (mask & edges).bit_count()-loops
            if power < 0 or power > int(degree == 2):
                break
            if power:
                jet |= 1 << positions[j]
        else:
            supports[jet].append(mask)
    return supports


def _valuation(supports, faces, popcount, np):
    output = []
    for masks in supports:
        value = np.full(len(faces), np.inf)
        # Keep the largest temporary below about8 MiB at the current12 edges.
        for begin in range(0, len(masks), 64):
            selected = np.asarray(masks[begin:begin+64], dtype=faces.dtype)
            value = np.minimum(value, popcount[faces[:, None] & selected].min(axis=1))
        output.append(value)
    return output


def _subsets(mask, proper=False):
    sub = (mask-1) & mask if proper else mask
    while sub:
        yield sub
        sub = (sub-1) & mask


def _composition_values(values, np):
    """Min-plus valuation of any nonempty ordered product of positive jets."""
    result = [np.zeros_like(values[0])]
    for mask in range(1, len(values)):
        value = np.full_like(values[0], np.inf)
        for sub in _subsets(mask):
            value = np.minimum(value, values[sub]+result[mask ^ sub])
        result.append(value)
    return result


def _density_valuation(u, f, legs, np):
    """Lower valuations of all retained coefficients of U^-3 Phi(F/U).

    This is finite min-plus Taylor algebra. Cancellations may raise the true
    order and are intentionally not inferred from a few sampled points.
    """
    u0 = u[0]
    if not np.all(np.isfinite(u0)):
        raise ArithmeticError('A recursive word has a singular constant Gaussian program')
    h = []
    for mask in range(len(u)):
        value = f[mask].copy()
        for sub in _subsets(mask):
            value = np.minimum(value, u[sub]+h[mask ^ sub])
        h.append(value-u0)
    if np.any(h[0] < 0):
        raise ArithmeticError('Euclidean Schur complement lost its bounded constant valuation')
    normalized_u = [np.full_like(u0, np.inf)]+[v-u0 for v in u[1:]]
    negative_power = [value-3*u0 for value in _composition_values(normalized_u, np)]
    h_products = _composition_values(h, np)
    if legs == 3:
        phi = [h[0]]+h_products[1:]
    else:
        phi = [2*h[0]]
        for mask in range(1, len(h)):
            # Phi2'=log1p(H0/M), while all derivatives of order>=2
            # have valuation0 at bounded nonnegative H0 and positive mass.
            value = h[0]+h[mask]
            for sub in _subsets(mask, proper=True):
                value = np.minimum(value, h[sub]+h_products[mask ^ sub])
            phi.append(value)
    total = np.full_like(u0, np.inf)
    for mask in range(len(u)):
        sub = mask
        while True:
            total = np.minimum(total, negative_power[sub]+phi[mask ^ sub])
            if not sub:
                break
            sub = (sub-1) & mask
    return total


def _distribution_from_gaps(gaps, np):
    size = (len(gaps)-1).bit_length()
    full = len(gaps)-1
    values = np.zeros(full+1)
    values[0] = 1.
    cdf = np.zeros((full+1, size))
    for mask in range(1, full+1):
        edges = [e for e in range(size) if mask & (1 << e)]
        weights = np.asarray([values[mask ^ (1 << e)] for e in edges])
        subtotal = float(weights.sum())
        values[mask] = subtotal if mask == full else subtotal/gaps[mask]
        row = np.zeros(size)
        row[edges] = weights/subtotal
        cdf[mask] = np.cumsum(row)
        cdf[mask, edges[-1]:] = 1.
    return cdf, float(values[full])


@lru_cache(maxsize=128)
def _sector_data(key, gap_floor, gap_scale):
    import numpy as np
    graph, parts, words = key
    size = len(graph.edges)
    if not isfinite(gap_floor) or not .25 <= gap_floor <= 2:
        raise ValueError('Gap floor must be finite and between.25 and2')
    if not isfinite(gap_scale) or not .25 <= gap_scale <= 2:
        raise ValueError('Gap scale must be finite and between.25 and2')
    if len(graph.external) == 1:
        # The declared full one-point subtraction gives an exact zero density;
        # any valid proposal can be used. Keep a simple normalized Hepp map.
        raw = np.asarray([mask.bit_count() for mask in range(1 << size)], dtype=float)
        minimum = raw.copy()
        supports = ((), ())
    else:
        supports = graph_polynomial_support(graph)
        faces = np.arange(1 << size, dtype=np.uint32)
        popcount = np.asarray([m.bit_count() for m in range(1 << size)], dtype=np.int16)
        minimum = np.full(len(faces), np.inf)
        for word in words:
            u = _valuation(_word_supports(parts, word, supports[0]), faces, popcount, np)
            f = _valuation(_word_supports(parts, word, supports[1]), faces, popcount, np)
            minimum = np.minimum(minimum, _density_valuation(u, f, len(graph.external), np))
        raw = popcount+minimum
    gaps = np.maximum(gap_floor, gap_scale*raw)
    gaps[0] = 0.
    gaps[-1] = 0.
    cdf, normalization = _distribution_from_gaps(gaps, np)
    from .automatic_sampling import sector_distribution
    old_gaps, _, _, _ = sector_distribution(graph)
    changed = [int(mask) for mask in range(1, len(gaps)-1) if gaps[mask] < old_gaps[mask]]
    metadata = {
        'method': 'recursive Gaussian polynomial-support Hepp importance proposal',
        'gap_rule': 'max(gap_floor, gap_scale*(|A|+min_word min_retained_coefficient valuation_A))',
        'support_rule': 'tree/two-tree complements; word exponent |B intersection gamma|-loops(gamma); u_gamma=t_gamma^2',
        'gap_floor': gap_floor, 'gap_scale': gap_scale,
        'proper_subsets': len(gaps)-2, 'recursive_terms': len(words),
        'tree_monomials': len(supports[0]), 'two_tree_monomials': len(supports[1]),
        'minimum_unclipped_gap': float(np.min(raw[1:-1])) if size > 1 else None,
        'faces_more_heavily_sampled_than_original_graph_hepp': len(changed),
        'first_changed_edge_sets': [[e for e in range(size) if mask & (1 << e)] for mask in changed[:24]],
        'hepp_normalization': normalization,
        'rigorous_recursive_remainder_certificate': False,
        'finite_variance_proved': False,
        'normalization_scope': 'Exact probability change of variables; DP constants are evaluated in floating point',
    }
    for value in (gaps, cdf, raw, minimum):
        value.setflags(write=False)
    return gaps, cdf, normalization, raw, minimum, metadata


def recursive_sector_data(compiled, *, gap_floor=.5, gap_scale=1.):
    compiled = _compiled(compiled)
    return _sector_data(_program_key(compiled), float(gap_floor), float(gap_scale))


def _open_uniforms(uniforms, dimension, np):
    points = np.asarray(uniforms, dtype=float)
    if points.ndim != 2 or points.shape[1] != dimension:
        raise ValueError(f'Expected a two-dimensional array with{dimension} uniform coordinates')
    if not np.all(np.isfinite(points)) or np.any(points < 0) or np.any(points > 1):
        raise ValueError('Uniform coordinates must lie in the finite closed unit cube')
    return np.clip(points, np.nextafter(0., 1.), np.nextafter(1., 0.))


def _hepp_map(gaps, cdf, normalization, points, np):
    size = (len(gaps)-1).bit_length()
    count = len(points)
    if size == 1:
        return np.ones((count, 1)), np.zeros(count)
    mask = np.full(count, (1 << size)-1, dtype=int)
    order = np.empty((count, size), dtype=int)
    for position in range(size-1, 0, -1):
        choice = np.sum(points[:, size-1-position, None] > cdf[mask], axis=1)
        order[:, position] = choice
        mask ^= 1 << choice
    order[:, 0] = np.rint(np.log2(mask)).astype(int)
    prefix = np.zeros(count, dtype=int)
    powers = np.empty((count, size-1))
    for r in range(1, size):
        prefix |= 1 << order[:, r-1]
        powers[:, r-1] = gaps[prefix]
    logt = np.log(points[:, size-1:])/powers
    ordered_log = np.zeros((count, size))
    ordered_log[:, :-1] = np.cumsum(logt[:, ::-1], axis=1)[:, ::-1]
    maximum = ordered_log.max(axis=1, keepdims=True)
    normalized = np.exp(ordered_log-maximum)
    total = normalized.sum(axis=1)
    alpha = np.empty((count, size))
    alpha[np.arange(count)[:, None], order] = normalized/total[:, None]
    logweight = log(normalization)-size*np.log(total)+np.sum((np.arange(1,size)[None,:]-powers)*logt, axis=1)
    if not np.all(np.isfinite(alpha)) or np.any(alpha <= 0):
        raise ArithmeticError('The proposal produced an unrepresentable positive simplex point')
    return alpha, logweight


def hepp_logpdf(alpha, gaps, normalization):
    """Evaluate a Hepp proposal at arbitrary positive original-edge parameters."""
    import numpy as np
    points = np.asarray(alpha, dtype=float)
    size = (len(gaps)-1).bit_length()
    if points.ndim != 2 or points.shape[1] != size or np.any(points <= 0) or not np.all(np.isfinite(points)):
        raise ValueError('Positive finite simplex points are required')
    if not np.allclose(points.sum(axis=1), 1., rtol=0., atol=2e-12):
        raise ValueError('Simplex points must sum to one')
    if size == 1:
        return np.zeros(len(points))
    order = np.argsort(points, axis=1, kind='stable')
    ordered = np.take_along_axis(points, order, axis=1)
    logt = np.log(ordered[:, :-1])-np.log(ordered[:, 1:])
    prefixes = np.bitwise_or.accumulate(1 << order[:, :-1], axis=1)
    powers = gaps[prefixes]
    # S=sum_e(alpha_e/alpha_max)=1/alpha_max.
    return -log(normalization)-size*np.log(ordered[:, -1])+np.sum((powers-np.arange(1,size)[None,:])*logt, axis=1)


def recursive_hepp_samples(compiled, uniforms, *, gap_floor=.5, gap_scale=1.):
    import numpy as np
    compiled = _compiled(compiled)
    data = recursive_sector_data(compiled, gap_floor=gap_floor, gap_scale=gap_scale)
    points = _open_uniforms(uniforms, 2*(len(compiled.graph.edges)-1), np)
    alpha, logweight = _hepp_map(*data[:3], points, np)
    weights = np.exp(logweight)
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ArithmeticError('The proposal inverse density is unrepresentable')
    return alpha, weights, dict(data[-1])


def recursive_mixture_samples(compiled, uniforms, *, gap_floor=.5, gap_scale=1.,
                              mixture_weights=(.8,.1,.1), dirichlet_beta=.5):
    """Mixture of recursive Hepp, original-graph Hepp and symmetric Dirichlet.

    Every returned point is weighted by the *total mixture density*, not by
    its selected component. The first uniform selects the component and the
    remaining coordinates generate that component without reusing it.
    """
    import numpy as np
    from scipy.special import gammaincinv, logsumexp
    compiled = _compiled(compiled)
    size = len(compiled.graph.edges)
    dimension = 1+max(2*(size-1), size)
    points = _open_uniforms(uniforms, dimension, np)
    probabilities = np.asarray(mixture_weights, dtype=float)
    if (probabilities.shape != (3,) or not np.all(np.isfinite(probabilities))
            or np.any(probabilities < 0) or probabilities.sum() <= 0):
        raise ValueError('Use three nonnegative finite mixture probabilities')
    if not isfinite(dirichlet_beta) or not .15 <= dirichlet_beta <= 2:
        raise ValueError('Dirichlet beta must be between.15 and2')
    probabilities = probabilities/probabilities.sum()
    choice = np.searchsorted(np.cumsum(probabilities), points[:, 0], side='right')
    choice = np.minimum(choice, 2)
    recursive = recursive_sector_data(compiled, gap_floor=gap_floor, gap_scale=gap_scale)
    from .automatic_sampling import sector_distribution
    old = sector_distribution(compiled.graph)
    alpha = np.empty((len(points), size))
    for channel in range(3):
        selected = choice == channel
        if not np.any(selected):
            continue
        u = points[selected, 1:]
        if channel < 2:
            data = recursive[:3] if channel == 0 else old[:3]
            alpha[selected] = _hepp_map(*data, u[:, :2*(size-1)], np)[0]
        else:
            gamma = gammaincinv(dirichlet_beta, u[:, :size])
            alpha[selected] = gamma/gamma.sum(axis=1, keepdims=True)
    logpdfs = []
    for channel, probability in enumerate(probabilities):
        if probability == 0:
            continue
        if channel < 2:
            data = recursive if channel == 0 else old
            density = hepp_logpdf(alpha, data[0], data[2])
        else:
            density = (lgamma(size*dirichlet_beta)-size*lgamma(dirichlet_beta)
                       +(dirichlet_beta-1)*np.log(alpha).sum(axis=1))
        logpdfs.append(log(probability)+density)
    logdensity = logsumexp(np.asarray(logpdfs), axis=0)
    weights = np.exp(-logdensity)
    if not np.all(np.isfinite(alpha)) or np.any(alpha <= 0) or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ArithmeticError('The mixture point or inverse density is unrepresentable')
    metadata = dict(recursive[-1])
    metadata.update(mixture_weights=probabilities.tolist(), dirichlet_beta=dirichlet_beta,
                    mixture_components=['recursive Hepp','original graph Hepp','symmetric Dirichlet'],
                    selected_channel_counts=[int(np.sum(choice == j)) for j in range(3)],
                    weighting='inverse of complete mixture PDF; selected-channel probability included',
                    sobol_dimension=dimension)
    return alpha, weights, metadata
