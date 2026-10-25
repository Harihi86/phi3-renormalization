"""Independent full Feynman pair and triangle parameter integrals for G46."""
import mpmath as mp
import mixed_weighted_triangle_pair_v1 as M
def evaluate(points,binding=None,digits=640,chart='mixed-triangle-pair-cubic',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert binding['id']=='L4E2G46' and chart in ('mixed-triangle-pair-cubic','mixed-triangle-pair-diagonal')
    return ctx,M.evaluate(ctx,list(map(number,points)),chart.rsplit('-',1)[1],quadrature)
