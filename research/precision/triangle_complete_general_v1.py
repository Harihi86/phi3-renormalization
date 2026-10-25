"""Exact massive triangle scalar for arbitrary positive Euclidean invariants.

The Gram determinant must be strictly positive. Degenerate c=0,a=b is
handled separately and exactly, never by clipping a nearly degenerate input.
"""
from nested_triangle_complete_v1 import edge_log, log_linear_pair


def inverse_triangle(ctx, a, b, c):
    alpha = (a-b)**2 - 2*c*(a+b) + c*c
    beta = 2*c*(a+b-c)
    assert a > 0 and b > 0 and c > 0 and alpha < 0
    root = ctx.sqrt(c*(c+4))
    kappa = ctx.sqrt(-alpha)
    z0 = beta/(2*root)
    q1 = ctx.sqrt((a-b)**2+4*c)
    delta = (4*root*root*alpha-beta*beta)/(2*root*(2*root*q1+2*root*root+beta))
    assert delta < 0
    kernel = delta/(z0-ctx.j*kappa)
    result = 2*ctx.log((root+c)/(root-c))*ctx.log(1+kernel)
    for linear, square in ((a-b-c, b), (-a+b-c, a)):
        for sign in (1, -1):
            leading = sign*c-root
            middle = beta-2*root*sign*linear
            discriminant = 4*c*ctx.sqrt(square*(square+4))
            for side in (1, -1):
                zero = (-middle+side*discriminant)/(2*leading)
                factor = delta/(z0-zero)
                assert 1+factor > 0
                result += sign*log_linear_pair(ctx, factor, kernel)
    return -2*ctx.im(result)/kappa


def triangle_log(ctx, a, b, c):
    assert a >= 0 and b >= 0 and c >= 0
    if c == 0:
        assert a == b
        return edge_log(ctx, a)/2
    cross = (a+b-c)/2
    determinant = a*b-cross*cross
    assert a > 0 and b > 0 and determinant > 0
    cx = b*(a-cross)/(2*determinant)
    cy = a*(b-cross)/(2*determinant)
    height = 1+(a*cx+b*cy)/2
    boundary = cx*edge_log(ctx, b)+cy*edge_log(ctx, a)+(1-cx-cy)*edge_log(ctx, c)
    return boundary/2-ctx.mpf('.5')+height*inverse_triangle(ctx, a, b, c)
