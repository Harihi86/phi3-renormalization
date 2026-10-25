"""Independent closed scalar reference, with a separate 1024-term positive-series check."""
import mpmath as mp
import double_weighted_angular_v1 as M

def evaluate(points,binding=None,digits=640,chart='angular-two-linear',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert chart in ('angular-two-linear','angular-two-cubic')
    return ctx,M.evaluate(ctx,list(map(number,points)),chart.rsplit('-',1)[1],series=quadrature,terms=1024)
