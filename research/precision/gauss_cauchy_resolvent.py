"""Exact Gauss/Cauchy bound using positive even moment defects.

Let A be the normalized integral on [-1,1], G_n its normalized Gauss-
Legendre rule, and f analytic on |z-c| <= R with |f| <= M. Write h for
the physical half-width and q=h/R<1. Cauchy's estimate bounds the kth
coefficient of f(c+h*x) by M*q**k. Odd defects vanish by symmetry and
even defects d_m=(A-G_n)x**(2*m) are nonnegative: Hermite interpolation
at the n Gauss nodes has remainder f^(2*n)(xi)*pi_n(x)**2/(2*n)!, and
this derivative is nonnegative for every even monomial of degree >=2*n.
Consequently

 |(A-G_n)f| <= M * (A-G_n)[1/(1-q*q*x*x)].

For the two resolvents 1/(1 +/- q*x), their exact Hermite remainder is
pi_n(x)**2 / ((1 +/- q*x)*pi_n(1/q)**2). Symmetry of the roots makes
their averaged remainder

 pi_n(x)**2 / ((1-q*q*x*x)*pi_n(1/q)**2).

Since A[P_n**2]=1/(2*n+1), cancellation of the monic normalization gives

 |(A-G_n)f| <= M / ((2*n+1)*(1-q*q)*P_n(1/q)**2).

The physical tensor integral therefore has error <= volume*sum_j M_j
times this factor. The tensor operators have normalized norm one because
their weights are positive and sum to one. This bound uses no sampled
integral, order difference, or floating point. Geometry and actual Gauss
rules must be certified separately. All intermediate rationals retain the
existing 1,000,000-bit guard.
"""
from fractions import Fraction as Q
import theta_two_radial_majorant_audit_v5 as A


def legendre(n, x):
    assert type(n) is int and 0 <= n <= 64
    x = A.bounded(Q(x))
    if n == 0:
        return Q(1)
    before, current = Q(1), x
    for j in range(1, n):
        after = A.bounded(((2*j+1)*x*current-j*before)/(j+1))
        before, current = current, after
    return current


def factor(n, q):
    assert type(n) is int and 2 <= n <= 64
    q = A.bounded(Q(q))
    assert 0 < q < 1
    p = legendre(n, A.bounded(1/q))
    assert p > 0
    return A.bounded(1/A.bounded((2*n+1)*(1-q*q)*p*p))


def error(volume, moduli, ratios, orders):
    volume = A.bounded(Q(volume))
    assert volume > 0 and len(moduli) == len(ratios) == len(orders)
    assert len(orders) > 0
    parts = []
    for m, q, n in zip(moduli, ratios, orders):
        m = A.bounded(Q(m))
        assert m >= 0
        parts.append(A.bounded(volume*m*factor(n, q)))
    return A.bounded(sum(parts)), parts
