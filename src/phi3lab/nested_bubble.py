"""Recursive massive six-dimensional BPH integrals with a bubble insertion.

The inserted edge momentum k is independent of the physical external momenta.
Integrating its (1-t_k^2)-subtracted bubble first gives B_R(k^2), including the
finite -k^2/6 Taylor term. Overall subtraction then acts only on the remaining
outer propagators. This is the original recursive BPH scheme in that routing.

L2E2G2 is a one-dimensional integral; L2E3G4/G6/G7 are two-dimensional.
Both loop measures are d^6k/pi^3, with no graph weight, coupling or mu factor.
The angularly reduced densities have fixed signs and explicit amplitude bounds.
These bounds are not quadrature errors or twelve-dimensional momentum L1 bounds.
The separate absolute-convergence proof is in docs/nested-bubble-analysis.html.
"""
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from functools import lru_cache
from math import asinh, fsum, hypot, inf, isfinite, log, nextafter, sqrt

from .finite import _gauss_legendre
from .kinematics import kinematics_status, validate_three_point


def _number(value, name, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a real number, not a boolean')
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real number') from exc
    if not isfinite(value) or value < 0 or (positive and value == 0):
        qualifier = 'positive' if positive else 'nonnegative'
        raise ValueError(f'{name} must be finite and {qualifier}')
    return value


def _three_kinematics(a, b, c):
    return validate_three_point((a, b, c))


def _scaled_bubble(t, w):
    """w*B_R(t/w) at unit mass, for t=1-w and 0<w<1.

    This is an exact elementary formula away from r=0; the convergent Beta
    moment series near r=0 retains all finite Taylor subtraction terms.
    """
    r = t/w
    if r < 1:
        coefficient = 1/60
        terms = [coefficient]
        for n in range(2, 80):
            coefficient *= -r*(n-1)/(2*(2*n+3))
            terms.append(coefficient)
            if abs(coefficient) < 1e-18:
                break
        return t*t/w*fsum(terms)
    return ((t+4*w)**1.5/(3*sqrt(t))*asinh(sqrt(t/(4*w)))
            - 4*w/3 - 4*t/9)


@lru_cache(maxsize=128)
def _two_coefficient(sigma, order):
    """Dimensionless J_2/sigma^2, evaluated without squaring sigma."""
    terms = []
    for y, weight in _gauss_legendre(order):
        w = y**4
        t = 1-w
        q = sigma*w
        if q <= 1:
            radical = hypot(1-q, 2*sqrt(w*q))
            ratio = 2*w/(radical+1-q)  # d/q, continuous at q=0
        else:
            inverse_q = 1/q
            radical_over_q = hypot(1-inverse_q, 2*sqrt(w*inverse_q))
            ratio = (1-inverse_q+radical_over_q)/2
        d = q*ratio
        inverse = 1/(1+d)
        fraction = d*inverse
        numerator = ((1+2*w)*fraction*fraction
                     + 6*w*fraction*inverse + 3*w*inverse*inverse)
        remainder_coefficient = ratio*ratio*numerator/(w+d)
        terms.append(weight*(2/3)*y**3*t*t*_scaled_bubble(t, w)
                     * remainder_coefficient)
    return fsum(terms)


def nested_bubble_two_point(s, m2=1, *, order=64):
    """Finite L2E2G2 BPH amplitude at s=p_E^2, with mass dimension two.

    The overall subtraction is (1-t_p^2), at fixed inserted edge momentum k.
    order is an integer in 2..256. Compare orders for the chosen s/m2; a small
    refinement difference is a diagnostic, not a rigorous error estimate.
    """
    s = _number(s, 's')
    m2 = _number(m2, 'm2', positive=True)
    _gauss_legendre(order)
    if s and log(s)-log(m2) > log(100):
        from .bubble_extreme import two_point_scaled
        value = s*two_point_scaled(log(s)-log(m2), order)
        if not isfinite(value):
            raise ValueError('The finite nested bubble value exceeds the floating-point range; use scaled amplitude')
        return value
    sigma = s/m2
    if not isfinite(sigma):
        raise ValueError('s/m2 exceeds the supported floating-point range')
    if s == 0:
        return 0.0
    value = s*(sigma*_two_coefficient(sigma, order))
    if not isfinite(value):
        raise ValueError('The finite nested bubble value exceeds the floating-point range')
    return value


@lru_cache(maxsize=128)
def _three_normalized(a, b, c, order):
    # Adjacent-leg symmetry is implemented before caching as well as satisfied
    # by the Feynman parameter integral itself.
    nodes = _gauss_legendre(order)
    terms = []
    for y, weight_y in nodes:
        w = y**4
        t = 1-w
        inner = []
        for u, weight_u in nodes:
            external = u*a+(1-u)*b
            difference = u*(1-u)*c
            q = w*external
            radical = hypot(1-q, 2*w*sqrt(external), 2*sqrt(t*w*difference))
            if q <= 1:
                e = 4*(q+t*difference)/(radical+1-q)
                inverse_l = 1/(2+w*e)
                e_over_l = e*inverse_l
            else:
                inverse_l = (1/q)/(1+1/q+radical/q)
                e_over_l = ((q-1)*inverse_l+radical*inverse_l)/w
            inner.append(weight_u*e_over_l*(1+2*inverse_l))
        terms.append(-2*weight_y*y**3*t*t*_scaled_bubble(t, w)*fsum(inner))
    return fsum(terms)


def nested_bubble_three_point(a, b, c, m2=1, *, order=64):
    """Finite bubble-inserted triangle, with dimensionless BPH amplitude.

    a,b are squared physical momenta at the two ends of the inserted edge;
    c is the squared momentum at the opposite vertex. Thus G4=(s1,s2,s3),
    G6=(s1,s3,s2), G7=(s2,s3,s1). The overall subtraction has order zero.
    The invariants must have nonnegative Euclidean Gram determinant. order
    is an integer in 2..256; compare refinements for the selected invariants.
    """
    a, b, c = _three_kinematics(a, b, c)
    m2 = _number(m2, 'm2', positive=True)
    _gauss_legendre(order)
    scale=max(a,b,c)
    if scale and log(scale)-log(m2) > log(100):
        from .bubble_extreme import three_point_value
        x,y=sorted((a/scale,b/scale))
        return three_point_value(log(scale)-log(m2),(x,y,c/scale),order)
    normalized = a/m2, b/m2, c/m2
    if not all(isfinite(value) for value in normalized):
        raise ValueError('An invariant/m2 exceeds the supported floating-point range')
    if not any(normalized):
        return 0.0
    x, y = sorted(normalized[:2])
    value = _three_normalized(x, y, normalized[2], order)
    if not isfinite(value):
        raise ValueError('The finite nested triangle value exceeds the floating-point range')
    return value


def _bound_data(inputs, numerator, denominator, formula, sign, scope,
                global_bound, global_formula):
    # Decimal.from_float keeps the exact converted input numbers. Every
    # positive operation rounds upward, including the final float conversion.
    converted = {key: Decimal.from_float(value) for key, value in inputs.items()}
    precision = max(100, *(len(value.as_tuple().digits)+8 for value in converted.values()))
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        polynomial = numerator(converted)/denominator(converted)
        logarithmic = global_bound(converted)
        bound = min(polynomial, logarithmic)
        decimals = {'absolute_integral_bound': bound, 'polynomial_bound': polynomial,
                    'global_log_bound': logarithmic}
        floats = {}
        for key, decimal in decimals.items():
            upper = float(decimal)
            if not isfinite(upper):
                upper = None
            elif Decimal.from_float(upper) < decimal:
                upper = nextafter(upper, inf)
            floats[key] = upper
    return {
        **inputs, **floats,
        'decimal_upper_bounds': {key: str(value) for key, value in decimals.items()},
        'formula': 'min(polynomial_bound, global_log_bound)',
        'formulas': {'polynomial_bound': formula, 'global_log_bound': global_formula},
        'sign': sign, 'scope': scope,
        'normalization': 'Unweighted J; d^6k/pi^3 d^6q/pi^3; no 1/S or coupling',
        'not_claimed': ['quadrature error bound', 'the same constant for momentum-space L1',
                        'arbitrary-loop BPH convergence theorem'],
        'numeric_rounding': 'Upward Decimal arithmetic and outward binary64 conversion; None means float range exceeded.',
        'input_interpretation': 'Enclosures apply to the exact binary64-converted input numbers.',
    }


def _two_log_bound(d):
    s, mass = d['s'], d['m2']
    if not s:
        return Decimal(0)
    logarithm = (1+s/mass).ln().next_plus()
    # A lower enclosure is necessary for every non-exact denominator.
    with localcontext() as context:
        context.rounding = ROUND_FLOOR
        denominator = 12*(mass+s)
    return s*logarithm*logarithm/24+s*s*(logarithm+1)/denominator


def _three_log_bound(d):
    scale, mass = max(d['a'], d['b']), d['m2']
    if not scale:
        return Decimal(0)
    logarithm = (1+2*scale/mass).ln().next_plus()
    with localcontext() as context:
        context.rounding = ROUND_FLOOR
        denominator = 6*(mass+2*scale)
    return logarithm*logarithm/24+scale*(logarithm+1)/denominator


def nested_bubble_two_point_bounds(s, m2=1):
    """Proven positive-amplitude bounds, quadratic and O(s*log(s)^2)."""
    s = _number(s, 's')
    m2 = _number(m2, 'm2', positive=True)
    return _bound_data({'s': s, 'm2': m2}, lambda d: 25*d['s']*d['s'],
                       lambda d: 576*d['m2'], '25*s^2/(576*m2)', 'nonnegative',
                       'L2E2G2; absolute angularly reduced radial integral and BPH amplitude',
                       _two_log_bound, 's*L^2/24+s^2*(L+1)/(12*(m2+s)), L=log(1+s/m2)')


def nested_bubble_three_point_bounds(a, b, c, m2=1):
    """Proven negative-amplitude bounds, linear and O(log(scale)^2)."""
    status = kinematics_status((a, b, c))
    a, b, c = _three_kinematics(a, b, c)
    m2 = _number(m2, 'm2', positive=True)
    result = _bound_data({'a': a, 'b': b, 'c': c, 'm2': m2},
                       lambda d: 25*(d['a']+d['b']), lambda d: 576*d['m2'],
                       '25*(a+b)/(576*m2)', 'nonpositive',
                       'Bubble-inserted triangle; absolute angularly reduced radial/Feynman integral and BPH amplitude',
                       _three_log_bound, 'L^2/24+A*(L+1)/(6*(m2+2*A)), A=max(a,b), L=log(1+2*A/m2)')
    result['kinematics_status'] = status
    return result


def nested_bubble_two_point_scaled(s, m2=1, *, order=64):
    """J_2/s, allowing extreme positive scale ratios without ratio overflow."""
    s=_number(s,'s');m2=_number(m2,'m2',positive=True);_gauss_legendre(order)
    if not s:return 0.
    length=log(s)-log(m2)
    if length>log(100):
        from .bubble_extreme import two_point_scaled
        return two_point_scaled(length,order)
    return nested_bubble_two_point(s,m2,order=order)/s
