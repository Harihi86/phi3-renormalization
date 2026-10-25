"""Finite one-loop massive BPH benchmarks in six Euclidean dimensions.

The bare dimensionally regularized parameter formulas in the measure
    J = integral d^d k / pi^(d/2) product_e 1/(q_e^2+m2), d=6-2*eps,
are Gamma(-1+eps)*integral_0^1 F^(1-eps) for the bubble and
Gamma(eps)*integral_simplex F^(-eps) for the triangle. Applying respectively
(1-t_p^2) and (1-t_p^0) removes their local pole and local finite terms. Thus

 R J_bubble = integral_0^1 [ (m2+delta)*log1p(delta/m2)-delta ] dx,
 delta=s*x*(1-x), s=p_E^2;
 R J_triangle = -integral_simplex log1p((s1*x*y+s2*x*z+s3*y*z)/m2).

Here m2 denotes mass squared. No coupling, diagram weight, mu factor, or
(2*pi)^(-d) measure is included. These are actual finite one-loop integrals,
not finite two-loop results or general BPH convergence proofs. Deterministic
quadrature differences are convergence diagnostics, not rigorous error bounds.
"""
from functools import lru_cache
from math import cos, fsum, isfinite, log, log1p, pi
from .kinematics import validate_three_point


def _nonnegative(value, name):
    value = float(value)
    if not isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')
    return value


def _mass(value):
    value = _nonnegative(value, 'm2')
    if value == 0:
        raise ValueError('m2 must be positive')
    return value


@lru_cache(maxsize=16)
def _gauss_legendre(order):
    """Gauss-Legendre nodes and weights on [0,1], computed in binary64."""
    if type(order) is not int or not 2 <= order <= 256:
        raise ValueError('Quadrature order must be an integer between 2 and 256')
    nodes = [0.0] * order
    weights = [0.0] * order

    def polynomial(z):
        previous, current = 1.0, z
        for n in range(2, order + 1):
            previous, current = current, ((2*n-1)*z*current-(n-1)*previous)/n
        derivative = order * (z*current-previous)/(z*z-1)
        return current, derivative

    for i in range((order + 1)//2):
        z = cos(pi*(i + 0.75)/(order + 0.5))
        for _ in range(64):
            value, derivative = polynomial(z)
            next_z = z - value/derivative
            if abs(next_z-z) < 4e-16:
                z = next_z
                break
            z = next_z
        else:
            raise ArithmeticError('Legendre root iteration failed to converge')
        _, derivative = polynomial(z)
        weight = 1/((1-z*z)*derivative*derivative)
        nodes[i], nodes[order-1-i] = (1-z)/2, (1+z)/2
        weights[i] = weights[order-1-i] = weight
    return tuple(zip(nodes, weights))


def _bubble_series_coefficient(z):
    """Stable [(1+z)*log1p(z)-z]/z^2 for 0<=z<1/8."""
    term = 0.5
    terms = [term]
    for n in range(2, 80):
        term *= -z*(n-1)/(n+1)
        terms.append(term)
        if abs(term) <= 5e-18:
            break
    return fsum(terms)


def _log1p_ratio(delta, m2):
    ratio = delta/m2
    # Positive finite inputs can still overflow their ratio. In that regime
    # log(1+delta/m2) differs from log(delta)-log(m2) below binary64 precision.
    return log1p(ratio) if isfinite(ratio) else log(delta)-log(m2)


def bubble_bph(s, m2=1, *, order=64):
    """Finite (1-t_p^2)J_bubble for s=p_E^2>=0 and common m2>0.

    The result has dimension mass squared and begins at s^2/(60*m2). The
    optional order controls deterministic Gauss-Legendre quadrature.
    """
    s = _nonnegative(s, 's')
    m2 = _mass(m2)
    nodes = _gauss_legendre(order)
    if s == 0:
        return 0.0
    terms = []
    for x, weight in nodes:
        delta = s*x*(1-x)
        ratio = delta/m2
        if ratio < 0.125:
            # Compute delta^2/m2 as delta*ratio before the small-z series;
            # squaring ratio first can underflow even when the answer does not.
            value = delta*ratio*_bubble_series_coefficient(ratio)
        elif ratio <= 1:
            value = m2*((1+ratio)*log1p(ratio)-ratio)
        else:
            value = (m2+delta)*_log1p_ratio(delta, m2)-delta
        terms.append(weight*value)
    result = fsum(terms)
    if not isfinite(result):
        raise ValueError('The finite bubble value exceeds the floating-point range')
    return result


def triangle_bph(invariants=(1, 1, 1), m2=1, *, order=64):
    """Finite (1-t^0)J_triangle with F=m2+s1*x*y+s2*x*z+s3*y*z.

    Invariants are squared momenta of three real Euclidean vectors summing to
    zero. Shared exact-binary64 Gram validation uses a relative Cauchy
    rounding tolerance and canonical effective inputs at a boundary. A Duffy map from
    [0,1]^2 uses x=u,y=(1-u)*v,z=(1-u)*(1-v), Jacobian 1-u. The result is
    dimensionless; equal small invariants s give -s/(8*m2)+O(s^2).
    """
    invariants = validate_three_point(invariants)
    m2 = _mass(m2)
    nodes = _gauss_legendre(order)
    if not any(invariants):
        return 0.0
    s1, s2, s3 = invariants
    terms = []
    for u, weight_u in nodes:
        x = u
        for v, weight_v in nodes:
            y, z = (1-u)*v, (1-u)*(1-v)
            delta = fsum((s1*x*y, s2*x*z, s3*y*z))
            terms.append(-weight_u*weight_v*(1-u)*_log1p_ratio(delta, m2))
    result = fsum(terms)
    if not isfinite(result):
        raise ValueError('The finite triangle value exceeds the floating-point range')
    return result


def finite_benchmark():
    """Two reproducible one-loop results with three quadrature refinements."""
    records = []
    cases = [('L1E2G1', 'bubble', lambda n: bubble_bph(1, order=n),
              {'s': 1, 'm2': 1}, '(1-t_p^2) J', 'mass squared'),
             ('L1E3G1', 'triangle', lambda n: triangle_bph(order=n),
              {'invariants': [1, 1, 1], 'm2': 1}, '(1-t^0) J', 'dimensionless')]
    for graph_id, name, evaluator, kinematics, subtraction, dimension in cases:
        values = [{'order': order, 'value': evaluator(order)} for order in (16, 32, 64)]
        difference = abs(values[-1]['value']-values[-2]['value'])
        zero = bubble_bph(0) if name == 'bubble' else triangle_bph((0, 0, 0))
        records.append({
            'id': graph_id, 'name': name, 'dimension': 6,
            'kinematics': kinematics, 'subtraction': subtraction,
            'value': values[-1]['value'], 'mass_dimension': dimension,
            'quadrature': values, 'last_refinement_difference': difference,
            'difference_interpretation': ('Absolute difference between orders 32 and 64; '
                                          'a convergence diagnostic, not a rigorous error bound'),
            'checks': {'zero_subtraction_point': zero == 0.0,
                       'quadrature_refinement_stable': difference < 1e-12},
            'normalization': ('Euclidean integral d^6 k/pi^3; unit propagator powers; '
                              'no 1/S, coupling, or mu factor; positive common mass squared'),
            'method': ('Deterministic binary64 Gauss-Legendre quadrature; '
                       'the triangle uses a Duffy simplex map'),
            'limits': ('Actual finite one-loop massive BPH amplitude only. '
                       'No finite two-loop amplitude, external numerical integrator run, '
                       'or general convergence proof is claimed.'),
        })
    return records
