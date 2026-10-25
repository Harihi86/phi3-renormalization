"""High-precision shifted-weight reference, optionally integrating the complete triangle directly."""
import mpmath as mp
import dynamic_weighted_triangle_v1 as M


def evaluate(points,binding=None,digits=640,chart='dynamic-triangle-four',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert chart in ('dynamic-triangle-four','dynamic-triangle-eight') and binding['id']=='L4E2G17'
    return ctx,M.evaluate(ctx,list(map(number,points)),chart.rsplit('-',1)[1],quadrature)
