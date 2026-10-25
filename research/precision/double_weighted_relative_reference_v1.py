"""Independent original-routed transverse integral, retaining each forest's Taylor jet."""
import mpmath as mp
import double_weighted_relative_v1 as M

def evaluate(points,binding=None,digits=640,chart='relative-four',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert chart in ('relative-four','relative-eight')
    assert binding is not None and binding['id']=='L4E2G20'
    return ctx,M.evaluate(ctx,list(map(number,points)),binding['id'],chart.rsplit('-',1)[1],binding['paired_states'] if quadrature else None)
