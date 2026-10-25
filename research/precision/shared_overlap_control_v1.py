"""Complete G48 remainder with exactly joint-zero-mean tail controls.

The radial variable and the actual external Taylor subtraction are unchanged.
The nonzero joint tail mean is retained, never discarded.
"""
import shared_overlap_radial_v1 as M


def coefficients(ctx,c,v):
    sigma=v*(1-v);ell=ctx.log1p(sigma);c2=c*c;c4=c2*c2
    a=-(1-12*c2+16*c4)/2-sigma*((v*v+1)*(4*c2-1)+4*v*c2)-sigma*sigma/2
    b=(c4*(2*v**4+ctx.mpf(8)*v**3/3+4*v*v+8*v+16)
       +c2*((-4*v**4+4*v*v+8*v+4)*ell-8*v**4+v*v-4*v-12)
       +(ctx.mpf(3)*v**4/2-2*v**3-v*v/2-ctx.mpf(1)/2)*ell
       +v**4/2+v**3/2+v*v/4-v/2+1)
    mean=-ctx.mpf(16)/27+5*ctx.sqrt(5)*ctx.asinh(ctx.mpf(1)/2)/9
    return a,b,mean


def control(ctx,s,c,v):
    a,b,mean=coefficients(ctx,c,v)
    return s*s/(1+s)**5*(a*ctx.log1p(s)+b-mean)


def evaluate(ctx,point,direct=False,chart='linear'):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    assert chart in ('linear','square')
    if chart=='square':jac=2*x;x=x*x
    else:jac=ctx.mpf(1)
    s=x/(1-x);ds=jac/(1-x)**2;theta=ctx.pi*u/2;c=ctx.cos(theta)
    weight=ctx.mpf(8)/3*ctx.sin(theta)**4
    residual=M.remainder(ctx,s,c,v,direct)-control(ctx,s,c,v)
    return ds*s*s*M.bubble(ctx,s)**2/(2*(1+s)**3)*weight*residual
