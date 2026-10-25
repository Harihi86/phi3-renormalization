"""Exact finite angular averaging of the full G48 density, with the true mean retained."""
import shared_overlap_control_v1 as M


def evaluate(ctx,point,direct=False,chart='four'):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    count={'four':4,'eight':8}[chart]
    s=x/(1-x)**3;ds=(1+2*x)/(1-x)**4;total=ctx.mpf(0)
    for index in range(count):
        fraction=u/2+ctx.mpf(index)/count
        # At this angle the exact angular measure vanishes. The massive integrand is finite.
        if fraction==1:continue
        angle=ctx.pi*fraction;c=ctx.cos(angle);weight=ctx.mpf(8)/3*ctx.sin(angle)**4
        residual=sum(M.M.remainder(ctx,s,c,w,direct)-M.control(ctx,s,c,w) for w in (v,1-v))/2
        total+=weight*residual
    return ds*s*s*M.M.bubble(ctx,s)**2/(2*(1+s)**3)*total/count
