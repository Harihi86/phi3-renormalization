"""G1 shared weighted edge: complete unweighted parallelogram box and its actual Taylor coefficients."""
import shared_overlap_radial_v1 as T


def moments(ctx,s):
    j1=4*ctx.asinh(ctx.sqrt(s)/2)/ctx.sqrt(s*(s+4))
    j2=2*(1+j1)/(s+4);j3=(1+3*j2)/(s+4)
    return (1-j1)/s,(j1-j2)/s,(1-2*j1+j2)/s**2,(1-3*j1+3*j2-j3)/s**3


def mass_kinetic(ctx,s,h,v):
    m,i1,i2,i3=moments(ctx,s)
    return 2*v*m,-2*v*v*(i1-2*i2)+2*v**3*(i1-3*i2)+ctx.mpf(8)/3*h*h*v**3*i3


def box_slice(ctx,s,h,v,direct=False):
    if direct:
        def integrand(z):
            lam=1-z;sigma=z*lam;a=1+s*sigma+lam*v*(1-lam*v);b=2*v*h*sigma
            if not h:return 4*v*lam*sigma/a
            return lam*ctx.log1p(2*b/(a-b))/h
        return ctx.quad(integrand,[0,ctx.mpf('.5'),1])
    if h:return (T.triangle_slice(ctx,s,h,v)-T.triangle_slice(ctx,s,-h,v))/h
    a=s+v*v;b=s+v*(2*v-1);c=1+v*(1-v);delta=ctx.sqrt(b*b+4*a*c)
    q0=2/delta*(ctx.atanh((2*a-b)/delta)+ctx.atanh(b/delta))
    q1=(b*q0-ctx.log((c+b-a)/c))/(2*a)
    q2=(c*q0+b*q1-1)/a;q3=(c*q1+b*q2-ctx.mpf(1)/2)/a
    return 4*v*(q1-2*q2+q3)


def remainder(ctx,s,cosine,v,direct=False):
    h=ctx.sqrt(s)*cosine;mass,kinetic=mass_kinetic(ctx,s,h,v)
    return box_slice(ctx,s,h,v,direct)-mass-kinetic-2*T.bubble_parameter(ctx,v)/(1+s)**2


def evaluate(ctx,point,chart='eight',direct=False):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    count={'four':4,'eight':8}[chart];s=x/(1-x)**3;ds=(1+2*x)/(1-x)**4;total=ctx.mpf(0)
    for i in range(count):
        fraction=u/2+ctx.mpf(i)/count
        if fraction==1:continue
        angle=ctx.pi*fraction;co=ctx.cos(angle);weight=ctx.mpf(8)/3*ctx.sin(angle)**4
        total+=weight*(remainder(ctx,s,co,v,direct)+remainder(ctx,s,co,1-v,direct))/2
    return ds*s*s*T.bubble(ctx,s)**2/(2*(1+s)**3)*total/count
