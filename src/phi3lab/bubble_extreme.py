"""Log-radius integration of the same finite recursive BPH bubble kernels.

The variable is x=log(1+k^2/M); positive panels cover the physical layer at
x=log(S/M) without ever constructing S/M. Both endpoints are resolved even
when that ratio overflows binary64. The discarded radial tail has a separate
analytic bound; it is not a bound on the quadrature error.
"""
from functools import lru_cache
from math import ceil, exp, expm1, fsum, hypot, log, pi, sin, sqrt
from .finite import _gauss_legendre

TAIL = 48.0


def scaled_bubble_log(x):
    """B_R(r)/(1+r), r=exp(x)-1, in the original full Taylor scheme."""
    w, t = exp(-x), -expm1(-x)
    if t < .5:
        r = t/w
        coefficient = 1/60
        terms = [coefficient]
        for n in range(2, 80):
            coefficient *= -r*(n-1)/(2*(2*n+3))
            terms.append(coefficient)
            if abs(coefficient) < 1e-20:
                break
        return t*t/w*fsum(terms)
    angle = x/2 + log((sqrt(t)+sqrt(t+4*w))/2)
    return (t+4*w)**1.5/(3*sqrt(t))*angle - 4*w/3 - 4*t/9


def log_radius_nodes(length, order):
    """Finite positive panels split at both the mass and external scales."""
    edges = sorted(set([0., length]
                       + [min(length, 4.*j) for j in range(1, ceil(length/4))]
                       + [length+4.*j for j in range(1, 13)]))
    nodes = _gauss_legendre(order)
    for left, right in zip(edges, edges[1:]):
        for u, weight in nodes:
            yield left+(right-left)*u, (right-left)*weight


@lru_cache(maxsize=64)
def two_point_scaled(length, order):
    """J_2/s, including the finite kinetic counterterm; length=log(s/M)."""
    terms = []
    for x, weight in log_radius_nodes(length, order):
        w, t = exp(-x), -expm1(-x)
        logq = length-x
        if logq <= 0:
            q = exp(logq)
            denominator = hypot(1-q, 2*sqrt(w*q))+1-q
            ratio = 2*w/denominator if denominator else 0.
            kappa = 2*q/denominator if denominator else 0.
            d_fraction = kappa/(1+kappa)
            d = q*ratio
            inverse = 1/(1+d)
            fraction = d*inverse
        else:
            iq = exp(-logq)
            ratio = (1-iq+hypot(1-iq, 2*sqrt(w*iq)))/2
            inverse = iq/(ratio+iq)
            fraction = ratio/(ratio+iq)
            d_fraction = ratio/(ratio+w*iq)
        numerator = ((1+2*w)*fraction*fraction
                     + 6*w*fraction*inverse + 3*w*inverse*inverse)
        terms.append(weight*t*t*scaled_bubble_log(x)*ratio*d_fraction*numerator/6)
    return fsum(terms)


@lru_cache(maxsize=64)
def three_point_value(length, shape, order):
    """Dimensionless J_3, with max(shape)=1 and length=log(S/M)."""
    a, b, c = shape
    nodes = tuple((sin(pi*u/2)**2, weight*pi*sin(pi*u)/2)
                  for u, weight in _gauss_legendre(order))
    terms = []
    for x, weight in log_radius_nodes(length, order):
        w, t = exp(-x), -expm1(-x)
        inner = []
        for u, wu in nodes:
            external = u*a+(1-u)*b
            if not external:
                inner.append(0.)
                continue
            difference_ratio = u*(1-u)*c/external
            logq = length-x+log(external)
            if logq <= 0:
                q = exp(logq)
                difference = q*difference_ratio
                radical = hypot(1-q, 2*sqrt(w*q), 2*sqrt(t*difference))
                denominator = radical+1-q
                increment = 4*(w*q+t*difference)/denominator if denominator else 0.
                inverse = 1/(2+increment)
                factor = increment*inverse*(1+2*inverse)
            else:
                iq = exp(-logq)
                radical_over_q = hypot(1-iq, 2*sqrt(w*iq),
                                       2*sqrt(t*difference_ratio*iq))
                inverse = iq/(1+iq+radical_over_q)
                factor = 1-4*inverse*inverse
            inner.append(wu*factor)
        terms.append(-weight*t*t*scaled_bubble_log(x)*fsum(inner)/2)
    return fsum(terms)


def radial_tail_bound(length, *, two_point=False):
    """Upper bound for |discarded tail| of J_3 or J_2/s only."""
    return (length+TAIL+1)*exp(-TAIL)/(12 if two_point else 6)
