"""Original box parameter integral or independent high-precision triangle primitives."""
import mpmath as mp
import dynamic_shared_box_v1 as M


def evaluate(points,binding=None,digits=640,chart='dynamic-box-four',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert chart in ('dynamic-box-four','dynamic-box-eight') and binding['id']=='L4E2G1'
    x,u,v=map(number,points);count=4 if chart.endswith('four') else 8
    s=x/(1-x)**3;ds=(1+2*x)/(1-x)**4;total=ctx.mpf(0)
    for i in range(count):
        fraction=u/2+ctx.mpf(i)/count
        if fraction==1:continue
        co=ctx.cospi(fraction);weight=ctx.mpf(8)/3*ctx.sinpi(fraction)**4
        total+=weight*(M.remainder(ctx,s,co,v,quadrature)+M.remainder(ctx,s,co,1-v,quadrature))/2
    return ctx,ds*s*s*M.T.bubble(ctx,s)**2/(2*(1+s)**3)*total/count
