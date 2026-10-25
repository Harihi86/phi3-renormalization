"""Exact semicircle marginals for the two cross-host propagators of G371."""


def moments(ctx, a, b):
    if not a > abs(b):
        raise ValueError('strict massive angular domain required')
    r = ctx.sqrt((a - b) * (a + b))
    return (2 / (a + r), 2 / (r * (a + r)), 1 / r**3, a / r**5)


def pair(ctx, a, c, b):
    if not min(a, c) > abs(b):
        raise ValueError('both massive angular denominators required')
    ra = ctx.sqrt((a - b) * (a + b))
    rc = ctx.sqrt((c - b) * (c + b))
    # Continuous at a == c; no numerical divided difference.
    return 2 * (1 + (a + c) / (ra + rc)) / ((a + ra) * (c + rc))


def angular_integral(ctx, function, a, b, c=None):
    """Independent real-angle quadrature, including near-endpoint peaks."""
    c = a if c is None else c
    gap = min(a - abs(b), c - abs(b)) / max(a, c)
    epsilon = min(ctx.mpf('.25'), ctx.sqrt(gap))
    nodes = sorted(set([ctx.mpf(0), epsilon, ctx.pi / 2,
                        ctx.pi - epsilon, ctx.pi]))
    return 2 / ctx.pi * ctx.quad(
        lambda theta: ctx.sin(theta)**2 * function(ctx.cos(theta)), nodes)
