"""G46: full weighted spherical pair times a finite triangle, with both parents."""
import spherical_pair_closed_v1 as P
import shared_overlap_radial_v1 as T


def triangle_jet(ctx,s,h,v):
    z=ctx.asinh(ctx.sqrt(s)/2);l0=ctx.sqrt((s+4)/s)*z-1
    j=4*z/ctx.sqrt(s*(s+4));j2=2*(1+j)/(s+4)
    i20=j/2-(1-j)/s;i30=j/2-3*(1-j)/(2*s)
    l2=v*i20-v*v*i30-h*h*v*v/(s*s)*(1-2*j+j2)
    return -l0,h*v*(1-j)/s,-l2


def pair_jet(ctx,s,t,h):
    k,r=P.G.shape(ctx,s,t);delta=ctx.sqrt((1+s+t)**2-4*s*t)
    u=1+(1+s-t)/delta
    first=-k**3*u;second=ctx.mpf(3)/2*k**4*u*u-4*t*k**3/delta**3
    j3=k*k/delta
    return k*k-1/(1+t)**2,-h*first,(first+2*h*h*second-j3)/3


def finite_pair(ctx,s,t,h,epsilon=1,direct=False):
    e=ctx.mpf(epsilon);u=s+e*e-2*e*h;dot=s-e*h
    if direct:value=P.feynman_reference(ctx,s,u,dot,e=t)
    else:
        ks,rs=P.G.shape(ctx,s,t);ku,ru=P.G.shape(ctx,u,t)
        pref=ks*ku;q=pref*t*dot;den=1-q
        # Exact geometric Gram formula, including collinear points; no subtraction
        # of two nearly equal O(s^2) Gram terms.
        gram=e*e*(s-h*h);assert gram>=0
        value=pref*(1-(rs*rs+ru*ru-2*q)/den*P.h_function(ctx,(pref*t)**2*gram/(den*den)))
    return value-1/(1+t)**2


def bubble_radial(ctx,s,direct=False):
    d=1+s
    if direct:
        angular=ctx.quad(lambda a:ctx.mpf(8)/(3*ctx.pi)*ctx.sin(a)**4/(d+1-2*ctx.sqrt(s)*ctx.cos(a)),[0,ctx.pi/2,ctx.pi])
    else:
        k,r=P.G.shape(ctx,s,ctx.mpf(1));angular=k*(1-r*r/3)
    kinetic=-1/d**2+2*s/(3*d**3)
    return (angular-1/d-kinetic)/d


def mass_kinetic(ctx,s,t,h,v):
    a=pair_jet(ctx,s,t,h);b=triangle_jet(ctx,s,h,v);d=1+s
    c0=a[0]*b[0];c1=a[0]*b[1]+a[1]*b[0];c2=a[0]*b[2]+a[1]*b[1]+a[2]*b[0]
    return c0/d**2,c2/d**2+2*h*c1/d**3+c0*(4*h*h/d**4-1/d**3)


def remainder(ctx,s,t,c,v,direct=False):
    h=ctx.sqrt(s)*c;d=1+s
    def product(sign):
        return finite_pair(ctx,s,t,sign*h,direct=direct)*T.triangle_slice(ctx,s,sign*h,v,direct=direct)/(d*(d+1-2*sign*h))
    mass,kinetic=mass_kinetic(ctx,s,t,h,v)
    parent3=pair_jet(ctx,s,t,h)[0]/d**3*T.bubble_parameter(ctx,v)
    parent4=triangle_jet(ctx,t,ctx.mpf(0),v)[0]/(1+t)**2*bubble_radial(ctx,s,direct)
    return (product(1)+product(-1))/2-mass-kinetic-parent3-parent4


def orbit(ctx,s,t,u,v,direct=False):
    value=ctx.mpf(0)
    for i in range(8):
        fraction=u/2+ctx.mpf(i)/8
        if fraction==1:continue
        c=ctx.cospi(fraction);weight=ctx.mpf(8)/3*ctx.sinpi(fraction)**4
        value+=weight*(remainder(ctx,s,t,c,v,direct)+remainder(ctx,s,t,c,1-v,direct))/16
    return value


def evaluate(ctx,point,chart='diagonal',direct=False):
    x,y,u,v=map(ctx.mpf,point);assert all(0<z<1 for z in (x,y,u,v))
    convert=lambda a:(a/(1-a)**3,(1+2*a)/(1-a)**4)
    a,da=convert(x);b,db=convert(y)
    weighted=lambda t:T.bubble(ctx,t)/(1+t)**2
    if chart=='cubic':return da*db*a*a*b*b*weighted(b)*orbit(ctx,a,b,u,v,direct)/4
    assert chart=='diagonal'
    s,t=(a+b)**2,a*a;jac=4*a*(a+b)*da*db
    return jac*s*s*t*t*(weighted(t)*orbit(ctx,s,t,u,v,direct)+weighted(s)*orbit(ctx,t,s,u,v,direct))/4
