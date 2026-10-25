"""Analytic bounds for the particular massive two-loop overlapping diamond.

J is unweighted, with two Euclidean loop measures d^6k/pi^3. M=m2>0,
s=p_E^2>=0. The BPH operation is (1-t_G^2)(1-t_gamma1^0-t_gamma2^0).
The proof controls the absolute integral of its COMBINED four-dimensional
parameter representation and the common-Schwinger-regulator limit. It does
not use this parameter L1 bound as a bound on the twelve-dimensional momentum
L1 norm. A separate momentum-space absolute-convergence proof is documented;
no general arbitrary-loop BPH theorem is proved here.

See docs/overlap-analysis.html for the connected derivation and exact constants.
"""
from decimal import Decimal, ROUND_CEILING, localcontext
from math import inf, isfinite, nextafter


CONSTANT_EXACT = '11*log(2)/24 - 29/384'


def _number(value, name, positive=False):
    if isinstance(value, bool):
        raise ValueError(f'{name} must be a real number, not a boolean')
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real number') from exc
    if not isfinite(number) or number < 0 or (positive and number == 0):
        qualifier = 'positive' if positive else 'nonnegative'
        raise ValueError(f'{name} must be finite and {qualifier}')
    return number


def _upper_float(value):
    """Convert a nonnegative Decimal enclosure to an upper binary64 value."""
    result = float(value)
    if not isfinite(result):
        return None
    if Decimal.from_float(result) < value:
        result = nextafter(result, inf)
    return result


def overlap_bound(s, m2=1):
    """JSON-ready proven bounds, including their exact mathematical formulas.

    C=11*log(2)/24-29/384, L=log(1+s/(4*m2)). The parameter L1 norm is
    at most min(C*s^2/m2, s*(7*L/4+L^2/4)). This also bounds |J|.
    Decimal computations round upward; correctly rounded ln is advanced by
    one Decimal ulp before subsequent upward-rounded arithmetic. Conversion
    to binary64 also rounds outward. Decimal strings remain available when a
    bound is outside the binary64 range (the corresponding float is None).
    This is a bound on amplitude magnitude, NOT a quadrature error estimate.
    """
    s = _number(s, 's')
    m2 = _number(m2, 'm2', positive=True)
    sd, md = Decimal.from_float(s), Decimal.from_float(m2)
    # Multiplication by four must keep the denominator exact; all binary64
    # values have a finite decimal expansion. Other positive arithmetic may
    # round upward without invalidating an upper bound.
    precision = max(100, len(md.as_tuple().digits)+4, len(sd.as_tuple().digits)+4)
    with localcontext() as context:
        context.prec = precision
        context.rounding = ROUND_CEILING
        log_two_upper = Decimal(2).ln().next_plus()
        # A lower enclosure for the subtracted positive rational is needed.
        rational_lower = (Decimal(29)/384).next_minus()
        constant = 11*log_two_upper/24-rational_lower
        if not s:
            scale = quadratic = logarithmic = first = Decimal(0)
            second = 2*constant/md
        else:
            scale = (1+sd/(4*md)).ln().next_plus()
            quadratic = constant*sd*sd/md
            logarithmic = sd*(7*scale/4+scale*scale/4)
            first = min(2*constant*sd/md, 1+7*scale/4+scale*scale/4)
            second = min(2*constant/md, (Decimal(5)/4+scale/2)/sd)
        combined = min(quadratic, logarithmic)
        values = {'quadratic_constant':constant, 'log_scale':scale,
                  'quadratic_bound':quadratic, 'global_log_bound':logarithmic,
                  'absolute_integral_bound':combined,
                  'first_s_derivative_bound':first,
                  'second_s_derivative_bound':second}
        result = {key:_upper_float(value) for key,value in values.items()}
        result['decimal_upper_bounds'] = {key:str(value) for key,value in values.items()}
    return {
        's':s, 'm2':m2, **result,
        'constant_exact':CONSTANT_EXACT,
        'formulas':{
            'quadratic':'C*s^2/m2',
            'global':'s*(7*L/4+L^2/4), L=log(1+s/(4*m2))',
            'first_s_derivative':'min(2*C*s/m2, 1+7*L/4+L^2/4)',
            'second_s_derivative':'min(2*C/m2, (5/4+L/2)/s) for s>0; 2*C/m2 at s=0',
            'compact_uniform':'For |p|<=P and m2>=m0_squared: replace s by P^2 and m2 by m0_squared in the amplitude bounds.'},
        'normalization':'Unweighted J; measure d^6k/pi^3 d^6q/pi^3; no 1/S or coupling',
        'scope':'This two-loop overlapping diamond; absolute combined parameter integral and common-Schwinger-regulator BPH limit',
        'not_claimed':['general all-loop BPH proof','termwise integration of unregulated divergent contributions',
                       'using the parameter bound C as the momentum-space L1 bound','quadrature error bound'],
        'numeric_rounding':'Upward Decimal enclosures and outward binary64 conversion; None means float range exceeded, with Decimal enclosure retained.',
        'input_interpretation':'Inputs are converted to binary64; the outward enclosures apply to those exact converted numbers.',
    }
