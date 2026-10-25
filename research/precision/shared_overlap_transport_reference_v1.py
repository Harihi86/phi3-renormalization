"""Independent scalar reference and optional complete triangle quadrature."""
import mpmath as mp
import shared_overlap_transport_v1 as M

def evaluate(points,binding=None,digits=640,chart='transport-overlap-quadratic',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert chart in ('transport-overlap-quadratic','transport-overlap-cubic')
    return ctx,M.evaluate(ctx,list(map(number,points)),quadrature,chart.rsplit('-',1)[1])
