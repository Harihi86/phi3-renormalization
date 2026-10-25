"""Exact angular moments and whole-axis control for two decorated core paths.

No graph certificate is claimed here.  The radial integral remains to be
computed.  Fractions are used throughout; no numerical fit fixes coefficients.
"""
from fractions import Fraction as Q
from math import factorial as fac, comb
from functools import lru_cache


def rising(x, n):
    answer = Q(1)
    for k in range(n):
        answer *= x + k
    return answer


def spectral_frame_coefficient(a, b):
    """Coefficient of x^a y^b in E exp(-2 sqrt(lambda) Tr(H^T R))."""
    n = a + b
    return rising(Q(5, 2), n) / (
        rising(Q(3), n) * fac(a) * fac(b)
        * rising(Q(5, 2), a) * rising(Q(5, 2), b))


@lru_cache(None)
def invariant_moment(n):
    """Coefficients of T^(n-2k) D^k, T=tr(H^T H), D=det(H^T H)."""
    assert isinstance(n, int) and 0 <= n <= 24
    residual = [spectral_frame_coefficient(n-b, b) for b in range(n+1)]
    coefficients = []
    for k in range(n//2+1):
        c = residual[k]
        coefficients.append(c)
        for b in range(k, n-k+1):
            residual[b] -= c * comb(n-2*k, b-k)
    assert not any(residual)
    return tuple(coefficients)


def moment(n, trace, determinant):
    return sum((c*trace**(n-2*k)*determinant**k
                for k, c in enumerate(invariant_moment(n))), Q())


def two_path_shape(tree):
    assert tree['kind'] == 'theta_K2_3' and len(tree['children']) == 2
    assert all(c['loops'] == 1 and not c['children'] for c in tree['children'])
    decorated = []
    plain = []
    for i, path in enumerate(tree['paths']):
        lengths = [len(a['edges']) for a in path['arcs']]
        if sorted(lengths) == [1, 2]:
            decorated.append((i, lengths.index(2)))
        else:
            assert lengths == [1, 1]
            plain.append(i)
    assert len(decorated) == 2 and len(plain) == 1
    (i, bi), (j, bj) = decorated
    return {'path_order': (i, j, plain[0]), 'bubble_sides': (bi, bj),
            'orientation': 'aligned' if bi == bj else 'opposite',
            'extra_multiplicity': 1}


def shifts(bi, bj):
    """External vectors in the (p_i,p_j) basis after shifting both bubbles."""
    assert bi in (0, 1) and bj in (0, 1)
    return ((1-2*bi, 0), (0, 1-2*bj), (-bi, -bj), (1-bi, 1-bj))


def vacuum_response_bound():
    return Q(41, 9) + Q(71, 3628800)


def integrated_remainder(n):
    assert isinstance(n, int) and n >= 4
    coefficients = [Q(comb(4, j), 1)/rising(Q(3), j) for j in range(5)]
    factor = Q(1, (n+1)*4**(n+1))
    for j in range(1, 5):
        factor += Q(5, 2)*coefficients[j]*fac(j-1)*comb(n, j-1)/4**(n-j+1)
    return vacuum_response_bound()*factor


@lru_cache(None)
def high_angular_numerator(k):
    """C_k(a,b)=P_k(a,b^2)/(a^2-b^2)^(k-5/2), k>=5.

    Dict keys are powers (a,b^2).  The derivative recurrence has positive
    coefficients and avoids subtracting nearly equal a^2 and b^2 terms.
    """
    assert 5 <= k <= 60
    if k == 5:
        return {(0, 0): Q(1)}
    prior = high_angular_numerator(k-1)
    h = 2*(k-1)-5
    out = {}
    for (r, s), c in prior.items():
        out[r+1, s] = out.get((r+1, s), Q()) + Q(h-r, k-1)*c
        if r:
            out[r-1, s+1] = out.get((r-1, s+1), Q()) + Q(r, k-1)*c
    assert all(c > 0 for c in out.values())
    return out
