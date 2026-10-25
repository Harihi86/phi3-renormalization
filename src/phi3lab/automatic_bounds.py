"""Exact graph-derived Hepp certificates for primitive three-point graphs.

The certificate controls both the absolute projective parameter integral and
the original momentum L1 norm of the overall-subtracted primitive BPH kernel.
It neither evaluates that integral nor supplies a quadrature-error estimate.
Subdivergent graphs are deliberately returned as
``needs_recursive_remainder``: a bare-period certificate cannot certify their
recursive cancellations.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import exp, factorial, isfinite, log, log1p

from .graphs import Graph, degrees, is_1pi
from .kinematics import kinematics_status, validate_three_point


def _graph(g):
    if not isinstance(g, Graph):
        g = Graph.from_dict(g)
    if (len(g.external) != 3 or g.loops < 1 or
            degrees(g) != (3,) * g.n or
            g.vertex_kinds != ('phi3',) * g.n or
            not is_1pi(range(g.n), g.edges)):
        raise ValueError('The primitive certificate requires a bare trivalent 1PI three-point graph')
    # This is a deliberate complexity limit, not a mathematical restriction.
    # It includes every graph in the current one-through-four-loop generator.
    if len(g.edges) > 18:
        raise ValueError('The exact subset certificate currently supports at most 18 edges')
    return g


def _cycle_rank(g, mask):
    """Cycle rank of an edge subset; full isolated vertices cancel from it."""
    parent = list(range(g.n))

    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    count = 0
    for e, (a, b) in enumerate(g.edges):
        if mask & (1 << e):
            count += 1
            a, b = root(a), root(b)
            if a != b:
                parent[a] = b
    return count - g.n + len({root(v) for v in range(g.n)})


@lru_cache(maxsize=256)
def _hepp_data(g):
    edge_count = len(g.edges)
    full = (1 << edge_count) - 1
    gaps = [0] * (full + 1)
    histogram = Counter()
    failing = []
    for mask in range(1, full):
        gap = mask.bit_count() - 3 * _cycle_rank(g, mask)
        gaps[mask] = gap
        histogram[gap] += 1
        if gap <= 0:
            failing.append(mask)
    if failing:
        return min(histogram), tuple(sorted(histogram.items())), tuple(failing), None

    # Each ordering of A ends in one uniquely chosen edge.  This recurrence
    # sums the products of all proper-prefix reciprocal gaps exactly.
    values = [Fraction(0)] * (full + 1)
    values[0] = Fraction(1)
    for mask in range(1, full):
        subtotal = Fraction(0)
        todo = mask
        while todo:
            bit = todo & -todo
            subtotal += values[mask ^ bit]
            todo ^= bit
        values[mask] = subtotal / gaps[mask]
    bound = sum((values[full ^ (1 << e)] for e in range(edge_count)), Fraction(0))
    return min(histogram), tuple(sorted(histogram.items())), (), bound


def primitive_uv_certificate(g):
    """Return an exact subset-DP certificate, or an explicit recursive-work status.

    For a primitive E=3 graph, ``B_exact`` is a rational upper bound for
    ``integral_(sum alpha=1) U(alpha)^(-3) d alpha``.  This is a bound on the
    period, not the value of the period.  No topology registry is consulted.
    """
    g = _graph(g)
    minimum, histogram, failing, bound = _hepp_data(g)
    result = {
        'schema': 'phi3lab.primitive-hepp-certificate.v1',
        'graph': g.to_dict(),
        'loops': g.loops,
        'edges': len(g.edges),
        'external_legs': 3,
        'dimension': 6,
        'proper_edge_subsets_checked': (1 << len(g.edges)) - 2,
        'gap_definition': 'd(A)=|A|-3 h(A), for every nonempty proper original edge subset A',
        'min_proper_d': minimum,
        'proper_gap_histogram': {str(k): v for k, v in histogram},
        'normalization': 'd^6k/pi^3 per loop; no coupling or symmetry weight',
        'certified_integral': 'absolute projective Schwinger-parameter integral after overall BPH subtraction',
        'original_momentum_L1_claimed': False,
        'quadrature_error_certified': False,
    }
    if failing:
        result.update({
            'status': 'needs_recursive_remainder',
            'primitive': False,
            'rigorous_analytic_certificate': False,
            'nonpositive_gap_subset_count': len(failing),
            'first_nonpositive_gap_edge_sets': [
                [e for e in range(len(g.edges)) if mask & (1 << e)]
                for mask in failing[:16]],
            'reason': 'A nonpositive proper UV gap prevents a bare-period bound; the recursive BPH sum must be combined and its remainders certified.',
        })
        return result
    result.update({
        'status': 'certified_primitive',
        'primitive': True,
        'rigorous_analytic_certificate': True,
        'B_exact': str(bound),
        'hepp_bound': str(bound),
        'B_numerator': bound.numerator,
        'B_denominator': bound.denominator,
        'B': float(bound),
        'period_bound_expression': 'integral_(sum alpha=1) U^(-3) <= B',
        'recurrence': 'D(empty)=1; D(A)=sum_(e in A) D(A without e)/d(A); B=sum_(e in E) D(E without e)',
        'proof_method': 'Kruskal tree monomial in each Hepp sector; exact rational subset dynamic programming',
        'original_momentum_L1_claimed': True,
        'original_momentum_L1_bound_expression': '8*L*B*P/sqrt(M), where max_i |p_i|<=P and m2>=M>0',
        'original_momentum_L1_proof': 'First external Taylor remainder; positive fractional-propagator Schwinger integral with one edge power 3/2; Gamma(1/2)/Gamma(3/2)=2.',
        'common_schwinger_cutoff_removal_proved': True,
        'floating_values_are_display_approximations': True,
    })
    return result


def _from_log(value):
    """JSON-safe display approximation; the log and analytic expression persist."""
    try:
        result = exp(value)
    except OverflowError:
        return None
    return result if isfinite(result) else None


def primitive_compact_bounds(g, invariants=(1, 1, 1), m2=1, *, derivative_order=2, scale=1):
    """Add uniform amplitude and derivative bounds to a primitive certificate.

    The nth scale derivative refers to ``J(scale*s1,scale*s2,scale*s3)`` with
    the supplied physical shape fixed.  Cartesian derivative bounds use
    momentum directions with max_i |v_i|<=1 and sum_i v_i=0.  The displayed
    real numbers are ordinary floating approximations to the stated analytic
    bounds, not outward-rounded interval enclosures.
    """
    g = _graph(g)
    shape = validate_three_point(invariants)
    if isinstance(m2, bool) or not isfinite(float(m2)) or m2 <= 0:
        raise ValueError('m2 must be positive and finite')
    if isinstance(scale, bool) or not isfinite(float(scale)) or scale < 0:
        raise ValueError('scale must be nonnegative and finite')
    if type(derivative_order) is not int or not 0 <= derivative_order <= 16:
        raise ValueError('derivative_order must be an integer from zero to sixteen')
    result = primitive_uv_certificate(g)
    result.update({'kinematics': kinematics_status(invariants), 'm2': float(m2),
                   'shape_smax': max(shape), 'scale': float(scale)})
    if not result['primitive']:
        return result
    bound = Fraction(result['B_numerator'], result['B_denominator'])
    log_b = log(bound.numerator) - log(bound.denominator)
    maximum = max(shape)
    # log(1+scale*S/M) without first forming an overflowing ratio/product.
    if maximum and scale:
        log_ratio = log(maximum) + log(scale) - log(m2)
        logarithm = max(log_ratio, 0) + log1p(exp(-abs(log_ratio)))
    else:
        logarithm = 0.0
    result.update({
        'absolute_integral_bound_expression': 'B*log(1+scale*S/M), S=max(s1,s2,s3), M=m2',
        'absolute_integral_bound': float(bound) * logarithm,
        'compact_uniformity': 'The same formulas hold uniformly for all physical momenta with max_i |p_i|^2<=S and all m2>=M.',
        'scale_derivative_bound_expression': 'B*(n-1)!*(S/(M+scale*S))^n, n>=1',
        'cartesian_direction_definition': 'sum_i v_i=0 and max_i |v_i|<=1; mixed second derivative permits two such directions',
        'cartesian_first_bound_expression': '2*B*sqrt(scale*S)/M',
        'cartesian_mixed_second_bound_expression': 'B*(2/M+4*scale*S/M^2)',
        'original_momentum_L1_bound_expression': '8*L*B*sqrt(scale*S/M)',
        'derivative_bounds_are_for_exact_integral': True,
        'floating_values_are_display_approximations': True,
    })
    scale_derivatives = {}
    if maximum:
        # log[S/(M+scale*S)] = -log(M/S+scale), evaluated by log-sum-exp.
        log_mass_ratio = log(m2) - log(maximum)
        if scale:
            a, b = log_mass_ratio, log(scale)
            log_denominator = max(a, b) + log1p(exp(-abs(a-b)))
        else:
            log_denominator = log_mass_ratio
        for n in range(1, derivative_order + 1):
            log_value = log_b + log(factorial(n-1)) - n*log_denominator
            scale_derivatives[str(n)] = {'bound': _from_log(log_value), 'log_bound': log_value}
    else:
        scale_derivatives = {str(n): {'bound': 0.0, 'log_bound': None}
                             for n in range(1, derivative_order+1)}
    result['common_scale_derivative_bounds'] = scale_derivatives
    if maximum and scale:
        first_log = log(2) + log_b + (log(scale)+log(maximum))/2-log(m2)
        result['cartesian_first_bound'] = _from_log(first_log)
        momentum_log = log(8)+log(g.loops)+log_b+(log(scale)+log(maximum)-log(m2))/2
        result['original_momentum_L1_bound'] = _from_log(momentum_log)
        result['original_momentum_L1_log_bound'] = momentum_log
    else:
        result['cartesian_first_bound'] = 0.0
        result['original_momentum_L1_bound'] = 0.0
        result['original_momentum_L1_log_bound'] = None
    # Add the positive terms without overflowing a physical input ratio.
    base = log(2)+log_b-log(m2)
    if maximum and scale:
        extra = log(4)+log_b+log(scale)+log(maximum)-2*log(m2)
        second_log = max(base, extra)+log1p(exp(-abs(base-extra)))
    else:
        second_log = base
    result['cartesian_mixed_second_bound'] = _from_log(second_log)
    result['cartesian_mixed_second_log_bound'] = second_log
    return result
