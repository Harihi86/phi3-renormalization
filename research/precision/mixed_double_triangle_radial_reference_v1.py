"""High-precision full original double-triangle remainder with independent parameter integrals."""
import mpmath as mp
import mixed_double_triangle_radial_v1 as M
def evaluate(points,binding=None,digits=640,chart='mixed-double-triangle-four',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert binding['id']=='L4E2G63' and chart in ('mixed-double-triangle-four','mixed-double-triangle-eight')
    return ctx,M.evaluate(ctx,list(map(number,points)),chart.rsplit('-',1)[1],quadrature)
