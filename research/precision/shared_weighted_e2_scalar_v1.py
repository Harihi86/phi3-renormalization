"""Shared B(s)^2/(1+s)^3 edge and its complete overall Taylor remainder."""
import weighted_bubble_e2_integral_reference_v1 as B

def hyper_bubble(ctx,s):
    s=ctx.mpf(s)
    assert s>=0
    return s*s/60*ctx.hyp2f1(1,1,ctx.mpf(7)/2,-s/4)

def direct_bubble(ctx,s):
    s=ctx.mpf(s)
    def f(x):
        a=x*(1-x)*s
        return (1+a)*ctx.log1p(a)-a
    return ctx.quad(f,[0,ctx.mpf('.5'),1])

def density(ctx,mass,edge,z,v,bubble=None):
    m,d,z,v=map(ctx.mpf,(mass,edge,z,v))
    assert m>d>0 and z>=0 and 0<=v<=1
    b=B.bubble if bubble is None else bubble
    a=m+v*z
    return d*d*z*z*(1-v)*b(ctx,(a-d)/d)**2/a**3

def partial_coefficients(ctx,s,t):
    """Five exact partial fractions, for distinct nonzero spectral parameters."""
    s,t=map(ctx.mpf,(s,t));assert 0<s<=ctx.mpf('.25') and 0<t<=ctx.mpf('.25') and s!=t
    # Obtain the repeated-pole coefficients from derivatives at x=-1.
    f=lambda x:x**4/((1+s*x)*(1+t*x))
    c3=f(-1);c2=ctx.diff(f,-1);c1=ctx.diff(f,-1,2)/2
    cs=1/((s-1)**3*(s-t));ct=1/((t-1)**3*(t-s))
    return c1,c2,c3,cs,ct

def spectral_primitive(ctx,mass,edge,z,s,t):
    m,d,z,s,t=map(ctx.mpf,(mass,edge,z,s,t));c1,c2,c3,cs,ct=partial_coefficients(ctx,s,t)
    def r(a):return (a+z)*ctx.log1p(z/a)-z
    rm=ctx.log1p(z/m)-z/m
    rmm=z*z/(m*m*(m+z))
    return c1*r(m)-c2*d*rm+c3*d*d*rmm/2+cs/s*r(m+d*(1/s-1))+ct/t*r(m+d*(1/t-1))

def spectral_remainder(ctx,mass,edge,z,s,t):
    m,d,z,s,t=map(ctx.mpf,(mass,edge,z,s,t))
    def f(v):
        x=(m+v*z-d)/d
        return (1-v)*x**4/((1+x)**3*(1+s*x)*(1+t*x))/d
    return z*z*ctx.quad(f,[0,ctx.mpf('.5'),1])
