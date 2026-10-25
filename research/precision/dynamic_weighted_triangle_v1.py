"""G17: original shifted weighted edge, complete unweighted triangle, and its external Taylor jets."""
import shared_overlap_radial_v1 as T


def weighted(ctx,s):return T.bubble(ctx,s)/(1+s)**2


def weighted_derivatives(ctx,s):
    d=1+s;root=ctx.sqrt(s);h=ctx.asinh(root/2);b=T.bubble(ctx,s)
    j0=4*h/ctx.sqrt(s*(s+4))
    bp=ctx.sqrt(s+4)*(s-2)*h/(3*s*root)+(s+4)/(6*s)-ctx.mpf(4)/9
    bpp=1/(6*s)-(1-j0)/(s*s)
    return b/d**2,bp/d**2-2*b/d**3,bpp/d**2-4*bp/d**3+6*b/d**4


def mass_kinetic(ctx,s,h,v):
    w,wp,wpp=weighted_derivatives(ctx,s);asinh=ctx.asinh(ctx.sqrt(s)/2)
    l0=ctx.sqrt((s+4)/s)*asinh-1;j0=4*asinh/ctx.sqrt(s*(s+4));j00=2*(1+j0)/(s+4)
    i20=j0/2-(1-j0)/s;i30=j0/2-3*(1-j0)/(2*s)
    l2=v*i20-v*v*i30-h*h*v*v/(s*s)*(1-2*j0+j00)
    return -w*l0,-w*l2-2*h*h*v*(1-j0)/s*wp-l0*(wp+2*h*h*wpp)


def remainder(ctx,s,cosine,v,direct=False):
    h=ctx.sqrt(s)*cosine
    paired=(weighted(ctx,s+1-2*h)*T.triangle_slice(ctx,s,h,v,direct=direct)+weighted(ctx,s+1+2*h)*T.triangle_slice(ctx,s,-h,v,direct=direct))/2
    mass,kinetic=mass_kinetic(ctx,s,h,v)
    return paired-mass-kinetic-weighted(ctx,s)/(1+s)*T.bubble_parameter(ctx,v)


def evaluate(ctx,point,chart='eight',direct=False):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    count={'four':4,'eight':8}[chart];s=x/(1-x)**3;ds=(1+2*x)/(1-x)**4;total=ctx.mpf(0)
    for i in range(count):
        fraction=u/2+ctx.mpf(i)/count
        if fraction==1:continue
        co=ctx.cospi(fraction);weight=ctx.mpf(8)/3*ctx.sinpi(fraction)**4
        total+=weight*(remainder(ctx,s,co,v,direct)+remainder(ctx,s,co,1-v,direct))/2
    return ds*s*s*weighted(ctx,s)*total/(2*count)
