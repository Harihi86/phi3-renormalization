"""Finite two-loop overlapping self energy by the recursive BPH R operation.

I = D(k)D(p-k)D(q)D(p-q)D(k-q), D(l)=1/(l*l+m2).
R I = (1-t_p^2)[I-D(k)^3 D(q)D(p-q)-D(q)^3 D(k)D(p-k)].
The two proper triangle subgraphs overlap; their counterterms never multiply.

Common Schwinger parameters allow the three terms to be combined before
integration. Both Gaussian loop integrals and the overall Schwinger scale are
integrated analytically. Projective pair sums, k<->q symmetry and two blow-up
charts leave bounded four-dimensional unit-cube densities. Every cancellation
at a logarithmically divergent subgraph face is performed algebraically.

J has Euclidean measure d^6k/pi^3 d^6q/pi^3, no 1/S, coupling or mu factors.
Results are numerical finite amplitudes. Refinement differences are diagnostics,
not rigorous quadrature-error bounds. The derivation and analytic domination
bound are documented separately; sampling does not prove convergence.
"""
from functools import lru_cache
from math import fsum, isfinite, log, log1p

from .finite import _gauss_legendre, _mass, _nonnegative, _bubble_series_coefficient


def _h(z):
    """[(1+z)log(1+z)-z]/z^2, continued by h(0)=1/2.

    h(z)=integral_0^1(1-v)/(1+z*v)dv lies between zero and 1/2.
    The small-z series avoids cancellation; the large-z form avoids z*z.
    """
    if z < 0.125:
        return _bubble_series_coefficient(z)
    if z > 1:
        return ((1+1/z)*log1p(z)-1)/z
    return ((1+z)*log1p(z)-z)/(z*z)


def _log1p_divided(z):
    return log1p(z)/z if z else 1.0


def _parameters(chart, t, y):
    if chart == 0:
        x, w, inverse_mass = t, t*y, y/(1+y)
    else:
        x, w, inverse_mass = t*y, t, 1/(1+y)
    r = 1-w
    c = 1+x*r
    d = 1+y-t*y
    return x, w, r, c, d, inverse_mass


def _scaled_density(chart, t, y, u, v, sigma):
    """Density divided by sigma squared, continued smoothly at sigma=0."""
    x, w, r, c, d, inverse_mass = _parameters(chart, t, y)
    a, b = u*(1-u), v*(1-v)
    tangent = a+r*(u-v)**2/c
    sb = sigma*b
    inverse_shift = inverse_mass/(1+sb*inverse_mass)
    divided_difference = (tangent*b*inverse_mass*_log1p_divided(sb*inverse_mass)
        + x*tangent*tangent*inverse_shift
          *_h(sigma*x*tangent*inverse_shift))
    denominator_difference = (r*(3+3*x*r+x*x*r*r)*b*b*inverse_mass
                              *_h(sb*inverse_mass))
    raw_minus_first = 2*t*w/c**3*(divided_difference-denominator_difference)
    second = (2*y*inverse_mass/d**3*a*a
              *_h(sigma*x*a*inverse_mass))
    return raw_minus_first-second


def overlap_sector_density(chart, t, y, u, v, sigma=1):
    """Continuous dimensionless density on either resolved unit-cube chart.

    Including endpoints is useful for inspecting the explicit UV-face limits.
    The full physical answer is m2 times the sum of both chart integrals at
    sigma=s/m2. This function is not the original momentum-space integrand.
    """
    if type(chart) is not int or chart not in (0, 1):
        raise ValueError('chart must be 0 or 1')
    coordinates = tuple(float(value) for value in (t, y, u, v))
    if any(not isfinite(value) or not 0 <= value <= 1 for value in coordinates):
        raise ValueError('Each chart coordinate must lie in [0,1]')
    sigma = _nonnegative(sigma, 'sigma')
    value = sigma*sigma*_scaled_density(chart, *coordinates, sigma)
    if not isfinite(value):
        raise ValueError('The chart density exceeds the floating-point range')
    return value


@lru_cache(maxsize=64)
def _normalized_overlap_coefficient(sigma, order):
    nodes = _gauss_legendre(order)
    uv = tuple((u, v, u*(1-u), v*(1-v), (u-v)**2, wu*wv)
               for u, wu in nodes for v, wv in nodes)
    outer = []
    for chart in (0, 1):
        for t, wt in nodes:
            for y, wy in nodes:
                x, w, r, c, d, inverse_mass = _parameters(chart, t, y)
                prefactor = 2*t*w/c**3
                subtract_prefactor = 2*y*inverse_mass/d**3
                denominator_factor = r*(3+3*x*r+x*x*r*r)
                inner = []
                for u, v, a, b, distance, weight in uv:
                    tangent = a+r*distance/c
                    sb = sigma*b
                    inverse_shift = inverse_mass/(1+sb*inverse_mass)
                    divided_difference = (tangent*b*inverse_mass*_log1p_divided(sb*inverse_mass)
                        + x*tangent*tangent*inverse_shift
                          *_h(sigma*x*tangent*inverse_shift))
                    denominator_difference = (denominator_factor*b*b*inverse_mass
                                              *_h(sb*inverse_mass))
                    value = prefactor*(divided_difference-denominator_difference)
                    value -= subtract_prefactor*a*a*_h(sigma*x*a*inverse_mass)
                    inner.append(weight*value)
                outer.append(wt*wy*fsum(inner))
    return fsum(outer)


def overlap_bph(s, m2=1, *, order=16):
    """Return the finite six-dimensional overlap R J at s=p_E^2>=0.

    order is the tensor Gauss-Legendre order per coordinate; total work scales
    as twice order**4. The standard-library backend intentionally caps order
    at32. Values at different orders should be compared for the chosen s/m2.
    Overall mass scaling is exact up to binary64 arithmetic: J=m2*f(s/m2).
    """
    s = _nonnegative(s, 's')
    m2 = _mass(m2)
    if type(order) is not int or not 2 <= order <= 32:
        raise ValueError('Quadrature order must be an integer between 2 and32')
    if s and log(s)-log(m2)>log(100):
        from .overlap_extreme import overlap_scaled
        value=s*overlap_scaled(log(s)-log(m2),order)
        if not isfinite(value):
            raise ValueError('The finite overlap value exceeds floating-point range; use overlap_bph_scaled')
        return value
    sigma = s/m2
    if not isfinite(sigma):
        raise ValueError('s/m2 exceeds the supported floating-point range')
    if sigma == 0:
        return 0.0
    # Keep the leading s^2/m2 factor as s*(s/m2). Squaring sigma first can
    # underflow while multiplication by m2 would have restored a finite value.
    value = s*(sigma*_normalized_overlap_coefficient(sigma, order))
    if not isfinite(value):
        raise ValueError('The finite overlap value exceeds the floating-point range')
    return value


def overlap_benchmark():
    """Reproducible full two-loop BPH result, not an integrand sample."""
    values = [{'order': order, 'value': overlap_bph(1, order=order)}
              for order in (8, 12, 16)]
    difference = abs(values[-1]['value']-values[-2]['value'])
    return {
        'id': 'L2E2G1', 'name': 'overlapping_two_point', 'dimension': 6,
        'kinematics': {'s': 1, 'm2': 1}, 'value': values[-1]['value'],
        'method': 'Recursive BPH; analytic Gaussian/radial integration; two resolved four-cube charts; tensor Gauss-Legendre',
        'subtraction': '(1-t_G^2)[I-t_gamma1^0 I-t_gamma2^0 I]',
        'forbidden_product': 'No t_gamma1*t_gamma2: the proper triangle subgraphs overlap',
        'recursive_counterterms': {
            'gamma1_edges':[0,1,4], 'gamma2_edges':[2,3,4],
            'gamma_orders':[0,0], 'overall_order':2,
            'proper_counterterm':'C_gamma,eta = - integral d^6l/pi^3 D_eta(l)^3',
            'overall_counterterm':'C_G,eta(p) = -RbarJ_eta(0) - s*d(RbarJ_eta)/ds(0)',
            'overall_local_basis':['constant (mass)','s=p_E^2 (kinetic)'],
            'renormalization_conditions':['J_R(0)=0','dJ_R/ds(0)=0'],
            'qualification':'Regulated counterterms are defined before the common cutoff is removed; divergent pieces are never separately evaluated without a regulator.'},
        'quadrature': values, 'last_refinement_difference': difference,
        'difference_interpretation': 'An order-refinement diagnostic, not a rigorous quadrature-error bound',
        'normalization': 'Euclidean d^6k/pi^3 d^6q/pi^3; no1/S, coupling, mu factors; diagram symmetry weight1/2 is external',
        'checks': {'zero_subtraction_point': overlap_bph(0) == 0,
                   'quadrature_refinement_stable': difference < 1e-10},
        'integrated': True,
        'limits': 'Finite two-loop result for this massive Euclidean graph and Taylor scheme. No arbitrary-loop theorem or general diagram integrator is claimed.',
    }


def overlap_bph_scaled(s,m2=1,*,order=16):
    """J_R/s for a positive external scale, without forming s/m2."""
    s=_nonnegative(s,'s');m2=_mass(m2)
    if type(order)is not int or not 2<=order<=32:
        raise ValueError('Quadrature order must be an integer between2and32')
    if not s:return 0.
    length=log(s)-log(m2)
    if length>log(100):
        from .overlap_extreme import overlap_scaled
        return overlap_scaled(length,order)
    return overlap_bph(s,m2,order=order)/s
