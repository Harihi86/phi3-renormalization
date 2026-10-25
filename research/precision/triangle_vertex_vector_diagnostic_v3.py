"""Small floating diagnostic of the proved L3E3G6 rank-one representation.

No acceptance or interval bounds are produced by this module. The independent
oracle evaluates both original positive U/F polynomials, without rank-one data.
"""
from fractions import Fraction as Q
import numpy as np
import triangle_vertex_rank_one_v1 as F


def mapped(points, permutation, transport):
    assert transport in ('uniform', 'endpoint')
    assert tuple(permutation) in F.PERMUTATIONS
    points = np.asarray(points, dtype=np.longdouble)
    assert points.ndim == 2 and points.shape[1] == 8
    assert np.isfinite(points).all() and ((points > 0) & (points < 1)).all()
    s, v, f0, f1, f2, r0, t, z0 = points.T
    den = 1+s+s*v
    weights = [None]*3
    for i, j in enumerate(permutation):
        weights[j] = (1/den, s/den, s*v/den)[i]
    a0, a1, a2 = weights
    psi = a0*a1+a0*a2+a1*a2
    fs = (f0, f1, f2)
    p0, p1 = f0-f2, f1-f2
    norm = lambda x, y: x*x+y*y-x*y
    q0 = sum(w*f*(1-f) for w, f in zip(weights, fs)) + a0*a1*a2*norm(p0, p1)/psi
    mu0, mu1 = a1*a2*p0/psi-f0, a1*a2*p1/psi
    if transport == 'endpoint':
        rd = r0*r0+(1-r0)**2
        r = r0*r0/rd
        complement_r = (1-r0)**2/rd
        dr = 2*r0*(1-r0)/rd**2
        z, dy = z0*z0, 2*z0
        y = (1-z0)*(1+z0)
    else:
        r, dr, z, dy = r0, np.ones_like(r0), z0, np.ones_like(z0)
        complement_r = 1-r0
        y = 1-z0
    a = r*complement_r
    d = complement_r*t*(1-t)
    n = a*norm(mu0+t, mu1)
    u = a*(a1+a2)/psi
    cc = 1+u
    sigma = y/(cc*z)
    d1, d2, d3 = u+z, 1+u*z, 1+z*(u+q0*cc)
    x0 = q0*cc*z/d2
    x1 = y*(n*cc*z/d1+d)/d3
    ca = u*(d1*d1+d1*cc*z+(cc*z)**2)/d1**3*q0*cc/d2
    cb = cc**3*z*z/d1**3*(n*cc*z/d1+d)/d3
    def l(x):
        out = np.ones_like(x)
        np.divide(np.log1p(x), x, out=out, where=x != 0)
        return out
    common = complement_r*dr*dy
    measure = v/(1+v+s*v)**3*common
    positive, negative = measure*ca*l(x0), measure*cb*l(x1)
    alphas = [weights[i]*(1-fs[i]) for i in range(3)] + [weights[i]*fs[i] for i in range(3)]
    alphas += [sigma*r, sigma*complement_r*t, sigma*complement_r*(1-t)]
    jac = a0*a1*a2*s/den**3*common*sigma*sigma/(cc*z*z)
    assert np.isfinite(positive).all() and np.isfinite(negative).all()
    return {'value': positive-negative, 'positive': positive, 'negative': negative,
            'alphas': np.asarray(alphas), 'jacobian': jac, 'mass_sum': 1+sigma}


def polynomial(rows, alphas, kind):
    total = np.zeros(alphas.shape[1], dtype=np.longdouble)
    for row in rows:
        mask, slot, channel = row[:3]
        assert slot == 0
        assert channel is None if kind == 'U' else type(channel) is int and channel in (0, 1, 2)
        coefficient = Q(row[3]) if len(row) == 4 else Q(1)
        value = np.full_like(total, np.longdouble(coefficient.numerator)/coefficient.denominator)
        for i, alpha in enumerate(alphas):
            if mask & (1 << i):
                value *= alpha
        total += value
    return total


def original(raw, chart):
    assert tuple(map(Q, raw['positive_invariants'])) == (Q(1), Q(1), Q(1))
    states = raw['positive_full_descriptor']['states']
    terms = []
    for term in raw['gaussian_full_descriptor']['terms']:
        state = next(s for s in states if s['word'] == term['subgraphs_inner_first'])
        psi = polynomial(state['U'], chart['alphas'], 'U')
        phi = polynomial(state['F'], chart['alphas'], 'F')
        assert (psi > 0).all() and (phi >= 0).all()
        terms.append(-term['sign']*np.log1p(phi/psi/chart['mass_sum'])/psi**3*chart['jacobian'])
    assert len(terms) == 2
    return {'value': sum(terms), 'absolute_terms': sum(np.abs(t) for t in terms)}


def condition(numerator, denominator):
    result = np.full_like(denominator, np.inf)
    np.divide(numerator, np.abs(denominator), out=result, where=denominator != 0)
    return result
