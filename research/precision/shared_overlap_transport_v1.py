"""Exact radial tail transport and v-reflection of the full controlled G48 integral."""
import shared_overlap_control_v1 as M


def evaluate(ctx,point,direct=False,chart='quadratic'):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    assert chart in ('quadratic','cubic');power={'quadratic':2,'cubic':3}[chart]
    s=x/(1-x)**power;ds=(1+(power-1)*x)/(1-x)**(power+1)
    theta=ctx.pi*u/2;c=ctx.cos(theta);weight=ctx.mpf(8)/3*ctx.sin(theta)**4
    residual=sum(M.M.remainder(ctx,s,c,w,direct)-M.control(ctx,s,c,w) for w in (v,1-v))/2
    return ds*s*s*M.M.bubble(ctx,s)**2/(2*(1+s)**3)*weight*residual
