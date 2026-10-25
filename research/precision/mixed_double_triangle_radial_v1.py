"""G63: two complete finite triangle slices with original parent bubbles and T2."""
import shared_overlap_control_v1 as A
T=A.M


def triangle_jet(ctx,s,h,v):
    z=ctx.asinh(ctx.sqrt(s)/2);l0=ctx.sqrt((s+4)/s)*z-1
    j=4*z/ctx.sqrt(s*(s+4));j2=2*(1+j)/(s+4)
    i20=j/2-(1-j)/s;i30=j/2-3*(1-j)/(2*s)
    l2=v*i20-v*v*i30-h*h*v*v/(s*s)*(1-2*j+j2)
    return -l0,h*v*(1-j)/s,-l2


def mass_kinetic(ctx,s,h,v,w):
    a=triangle_jet(ctx,s,h,v);b=triangle_jet(ctx,s,h,w);d=1+s
    c0=a[0]*b[0];c1=a[0]*b[1]+a[1]*b[0];c2=a[0]*b[2]+a[1]*b[1]+a[2]*b[0]
    return c0/d,c2/d+2*h*c1/d**2+c0*(4*h*h/d**3-1/d**2)


def remainder(ctx,s,c,v,w,direct=False,controlled=True):
    h=ctx.sqrt(s)*c;d=1+s;t0=triangle_jet(ctx,s,h,v)[0]
    positive=T.triangle_slice(ctx,s,h,v,direct=direct)*T.triangle_slice(ctx,s,h,w,direct=direct)/(d+1-2*h)
    negative=T.triangle_slice(ctx,s,-h,v,direct=direct)*T.triangle_slice(ctx,s,-h,w,direct=direct)/(d+1+2*h)
    mass,kinetic=mass_kinetic(ctx,s,h,v,w)
    value=(positive+negative)/2-mass-kinetic-t0/d**2*(T.bubble_parameter(ctx,v)+T.bubble_parameter(ctx,w))
    if controlled:value-=t0*(A.control(ctx,s,c,v)+A.control(ctx,s,c,w))
    return value


def evaluate(ctx,point,chart='eight',direct=False):
    x,u,v,w=map(ctx.mpf,point);assert all(0<z<1 for z in (x,u,v,w))
    count={'four':4,'eight':8}[chart];s=x/(1-x)**3;ds=(1+2*x)/(1-x)**4;value=ctx.mpf(0)
    for i in range(count):
        fraction=u/2+ctx.mpf(i)/count
        if fraction==1:continue
        c=ctx.cospi(fraction);weight=ctx.mpf(8)/3*ctx.sinpi(fraction)**4
        local=sum(remainder(ctx,s,c,a,b,direct) for a in (v,1-v) for b in (w,1-w))/4
        value+=weight*local/count
    return ds*s*s*T.bubble(ctx,s)/(2*(1+s)**2)*value
