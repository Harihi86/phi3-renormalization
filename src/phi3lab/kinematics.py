"""Shared binary64 Euclidean three-point input domain and Gram vectors.

All momenta are incoming: p1+p2+p3=0, si=|pi|^2.  Validation is
permutation invariant, including near a zero leg and at extreme scales.
Only a relative Cauchy rounding tolerance is allowed at a Gram boundary;
``kinematics_status`` reports when reconstruction projects to that boundary.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from math import isfinite, nextafter


_RELATIVE_TOLERANCE = Fraction(1, 10**14)


def _nonnegative_binary64(value):
    if isinstance(value, bool):
        raise ValueError('Three-point invariants must be real numbers, not booleans')
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Three-point invariants must be finite real numbers') from exc
    if not isfinite(result) or result < 0:
        raise ValueError('Three-point invariants must be finite and nonnegative')
    return result


@lru_cache(maxsize=256)
def _exact_gram(invariants):
    a, b, c = (Fraction.from_float(value) for value in invariants)
    # Exact binary64 arithmetic here is stronger than evaluating the dot
    # difference with fsum: neither cancellation nor overflow loses a tiny leg.
    difference = c-a-b
    gram = 4*a*b-difference*difference
    scales = (4*a*b+(c-a-b)**2,
              4*a*c+(b-a-c)**2,
              4*b*c+(a-b-c)**2)
    scale = min(scales)
    if gram < 0 and -gram > _RELATIVE_TOLERANCE*scale:
        raise ValueError('The invariants must describe real Euclidean momenta summing to zero; relative Gram condition failed')
    return (a, b, c), gram, scale


def _coerce_invariants(invariants):
    if isinstance(invariants, (str, bytes)):
        raise ValueError('Supply a sequence of three Euclidean invariants')
    try:
        result = tuple(_nonnegative_binary64(value) for value in invariants)
    except TypeError as exc:
        raise ValueError('Supply three Euclidean invariants') from exc
    if len(result) != 3:
        raise ValueError('Exactly three Euclidean invariants are required')
    return result


@lru_cache(maxsize=256)
def _canonical_invariants(invariants):
    _, gram, _ = _exact_gram(invariants)
    if gram >= 0:
        return invariants, 0
    # An invalid nonnegative triple has its largest length exceeding the sum
    # of the other two. Move that unique largest invariant inward, retaining
    # exact binary64 values. This rule commutes with every leg permutation.
    values = list(invariants)
    axis = max(range(3), key=values.__getitem__)
    for steps in range(1, 1025):
        values[axis] = nextafter(values[axis], 0.0)
        effective = tuple(values)
        if _exact_gram(effective)[1] >= 0:
            return effective, steps
    raise ValueError('Gram boundary correction exceeded the supported rounding range')


def validate_three_point(invariants):
    """Return a legal binary64 triple, applying a marked boundary correction.

    With G=4ab-(c-a-b)^2, an exact negative requested Gram value is accepted
    only if -G <= 1e-14*min_pairs[4ab+(c-a-b)^2]. Testing all three Cauchy
    pairs makes the domain independent of external-leg ordering. A zero leg
    therefore requires the other two binary64 invariants to be equal.
    For an accepted negative G, the largest invariant is moved toward zero
    by nextafter until the effective triple has an exact nonnegative Gram.
    Every kernel and bound must use this returned effective triple; call
    ``kinematics_status`` on the requested input to record the correction.
    """
    return _canonical_invariants(_coerce_invariants(invariants))[0]


def kinematics_status(invariants):
    """Describe exact-domain validity versus a tolerated boundary projection."""
    requested = _coerce_invariants(invariants)
    effective, steps = _canonical_invariants(requested)
    _, gram, scale = _exact_gram(requested)
    effective_gram = _exact_gram(effective)[1]
    projected = bool(steps)
    return {'exact_euclidean_gram': not projected,
            'gram_sign': (gram > 0)-(gram < 0),
            'effective_gram_sign': (effective_gram > 0)-(effective_gram < 0),
            'boundary_projection': projected,
            'changed': projected,
            'requested_invariants': list(requested),
            'effective_invariants': list(effective),
            'projection_steps': steps,
            'relative_cauchy_defect': float(-gram/scale) if projected else 0.0,
            'relative_gram_tolerance': float(_RELATIVE_TOLERANCE),
            'tolerance_scale': 'Minimum of the three exact binary64 Cauchy product scales',
            'reconstruction': ('Largest invariant moved inward by binary64 nextafter until the effective Gram matrix is exactly nonnegative'
                               if projected else 'Exact nonnegative Gram matrix reconstructed up to binary64 rounding')}


def _decimal(value):
    return Decimal(value.numerator)/Decimal(value.denominator)


@lru_cache(maxsize=256)
def _normalized_vectors(invariants):
    (a, b, c), gram, _ = _exact_gram(invariants)
    scale = max(a, b, c)
    if not scale:
        return 0.0, 0.0, 0.0, 0.0
    axis = max(a, b)
    # axis/scale >= 1/4 for a physical triple, so the axis never becomes a
    # tiny divisor.  Square roots are taken before converting to binary64:
    # sqrt(1e-616) is representable even though 1e-616 itself is not.
    with localcontext() as context:
        context.prec = 80
        axis_length = _decimal(axis/scale).sqrt()
        projection = _decimal((c-a-b)/(2*scale))/axis_length
        perpendicular = _decimal(gram/(4*axis*scale)).sqrt()
        axis_length, projection, perpendicular = (float(value) for value in
                                                  (axis_length, projection, perpendicular))
    if a >= b:
        return axis_length, 0.0, projection, perpendicular
    return projection, perpendicular, axis_length, 0.0


def normalized_gram_vectors(invariants):
    """Return (p1x,p1y,p2x,p2y), scaled by sqrt(max(si)).

    p3=-p1-p2.  The larger of p1 and p2 is the reconstruction axis.
    All-zero input returns four zeros. A tolerated negative requested Gram
    is first replaced by the same effective physical triple returned by
    ``validate_three_point``; call ``kinematics_status`` to record it.
    No normalization of the triple precedes exact Gram algebra.
    """
    return _normalized_vectors(validate_three_point(invariants))
