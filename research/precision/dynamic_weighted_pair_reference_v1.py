"""Independent high-precision closed pair and full Feynman-parameter reference."""
import mpmath as mp
import dynamic_weighted_pair_v1 as M
def evaluate(points,binding=None,digits=640,chart='dynamic-pair-cubic',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert binding['id']=='L4E2G8' and chart in ('dynamic-pair-cubic','dynamic-pair-diagonal')
    return ctx,M.evaluate(ctx,list(map(number,points)),chart.rsplit('-',1)[1],quadrature)
