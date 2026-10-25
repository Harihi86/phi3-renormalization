"""Exact first/second invariant derivatives by simplex boundary identities.

No finite differences, Gram clipping, integration nodes or added counterterms.
The independent validation differentiates the separately verified scalar.
"""
from math import comb
from triangle_complete_general_v1 import triangle_log, inverse_triangle, edge_log

WEIGHTS = ({(1, 0): 1, (2, 0): -1, (1, 1): -1},
           {(0, 1): 1, (1, 1): -1, (0, 2): -1}, {(1, 1): 1})


def edge_inverse(ctx, t):
    assert t >= 0
    if t == 0:
        return [ctx.mpf(1)/(k+1) for k in range(4)]
    j = 4*ctx.asinh(ctx.sqrt(t)/2)/ctx.sqrt(t*(t+4))
    return [j, j/2, j/2+(j-1)/t, j/2+ctx.mpf('1.5')*(j-1)/t]


def moments(ctx, a, b, c):
    d = (a+b-c)/2
    determinant = a*b-d*d
    assert a > 0 and b > 0 and c > 0 and determinant > 0
    inv = ctx.matrix([[b, -d], [-d, a]])/determinant
    centre = inv*ctx.matrix([a, b])/2
    height = 1+(a*centre[0]+b*centre[1])/2
    value = triangle_log(ctx, a, b, c)
    J = inverse_triangle(ctx, a, b, c)
    Ea, Eb, Ec = [edge_log(ctx, v) for v in (a, b, c)]
    first = centre*J-inv*ctx.matrix([Ec-Eb, Ec-Ea])/2
    boundary = ctx.matrix([[Ec/2, (Ec-Ea)/2], [(Ec-Eb)/2, Ec/2]])
    second = (first*ctx.matrix([[a, b]])+value*ctx.eye(2)-boundary)*inv/2
    M = {(0, 0): J, (1, 0): first[0], (0, 1): first[1],
         (2, 0): second[0, 0], (1, 1): second[0, 1], (0, 2): second[1, 1]}
    edges = [edge_inverse(ctx, v) for v in (a, b, c)]
    N = {(0, 0): (centre[0]*edges[1][0]+centre[1]*edges[0][0]
                  +(1-centre[0]-centre[1])*edges[2][0])/(2*height)}
    redundancies = [second[1, 0]-second[0, 1]]
    for degree in range(4):
        for m in range(degree+1):
            n = degree-m
            Bx = sum((-1)**k*comb(m, k)*edges[2][n+k] for k in range(m+1))
            By = sum((-1)**k*comb(n, k)*edges[2][m+k] for k in range(n+1))
            if m == 0: Bx -= edges[1][n]
            if n == 0: By -= edges[0][m]
            derivative = ctx.matrix([m*M[(m-1, n)] if m else 0,
                                     n*M[(m, n-1)] if n else 0])
            nexts = centre*N[m, n]+inv*(ctx.matrix([Bx, By])-derivative)/2
            for key, v in zip(((m+1, n), (m, n+1)), nexts):
                if key in N:
                    redundancies.append(N[key]-v)
                else:
                    N[key] = v
    return value, M, N, redundancies


def invariant_jets(ctx, a, b, c):
    value, M, N, redundant = moments(ctx, a, b, c)
    gradient = [sum(coef*M[key] for key, coef in poly.items()) for poly in WEIGHTS]
    hessian = []
    for left in WEIGHTS:
        row = []
        for right in WEIGHTS:
            row.append(-sum(c1*c2*N[m1+m2, n1+n2]
                            for (m1, n1), c1 in left.items()
                            for (m2, n2), c2 in right.items()))
        hessian.append(row)
    return value, gradient, hessian, redundant


def physical_jet(ctx, invariant, linear, quadratic):
    """Coefficient of lambda^k in -L(a+lambda*l+lambda^2*c), k=0..2."""
    value, g, h, unused = invariant_jets(ctx, *invariant)
    first = -sum(g[i]*linear[i] for i in range(3))
    second = -sum(g[i]*quadratic[i] for i in range(3))
    second -= sum(h[i][j]*linear[i]*linear[j] for i in range(3) for j in range(3))/2
    return [-value, first, second]
