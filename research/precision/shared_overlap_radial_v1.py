"""G48: analytic triangle parameter, paired external directions and local counterterm."""

def bubble(ctx,s):
    if not s:return ctx.mpf(0)
    return (s+4)*ctx.sqrt((s+4)/s)*ctx.asinh(ctx.sqrt(s)/2)/3-4*s/9-ctx.mpf(4)/3

def bubble_parameter(ctx,v):
    a=v*(1-v)
    return (1+a)*ctx.log1p(a)-a

def weighted_log(ctx,a,b,c):
    """Integral_0^1 (1-z) log(c+b*z-a*z^2) dz, for positive endpoint factors."""
    assert a>0 and c>0 and c+b-a>0
    disc=ctx.sqrt(b*b+4*a*c)
    # Product of roots = -c/a. Avoid subtracting nearly equal positive numbers.
    if b>=0:
        plus=(b+disc)/(2*a);minus=-c/(a*plus)
    else:
        minus=(b-disc)/(2*a);plus=-c/(a*minus)
    assert plus>1 and minus<0
    def f0(t):return t*ctx.log(t)-t
    def f1(t):return t*t*ctx.log(t)/2-t*t/4
    return (ctx.log(a)/2+(1-plus)*(f0(plus)-f0(plus-1))+f1(plus)-f1(plus-1)
        +(1-minus)*(f0(1-minus)-f0(-minus))-f1(1-minus)+f1(-minus))

def triangle_slice(ctx,s,h,v,epsilon=1,direct=False):
    """Minus the complete triangle z integral, retaining the bubble-pair parameter v."""
    e=ctx.mpf(epsilon);a=s+e*e*v*v-2*e*v*h
    b=s+e*e*v*(2*v-1)-2*e*v*h;c=1+e*e*v*(1-v)
    if direct:return -ctx.quad(lambda z:(1-z)*ctx.log(c+b*z-a*z*z),[0,ctx.mpf('.5'),1])
    return -weighted_log(ctx,a,b,c)

def taylor_coefficients(ctx,s,h,v):
    """The actual epsilon^0 and epsilon^2 coefficients, before angular averaging."""
    D=1+s;H=ctx.asinh(ctx.sqrt(s)/2)
    l0=ctx.sqrt((s+4)/s)*H-1
    j0=4*H/ctx.sqrt(s*(s+4));j00=2*(1+j0)/(s+4)
    i20=j0/2-(1-j0)/s;i30=j0/2-3*(1-j0)/(2*s)
    l2=v*i20-v*v*i30-h*h*v*v/(s*s)*(1-2*j0+j00)
    f0=-l0/D
    f2=-l2/D+2*h*h*v*(1-j0)/(s*D*D)-l0*(4*h*h/(D**3)-1/(D*D))
    return f0,f2

def remainder(ctx,s,cosine,v,direct=False):
    h=ctx.sqrt(s)*cosine;D=1+s
    positive=triangle_slice(ctx,s,h,v,direct=direct)/(D+1-2*h)
    negative=triangle_slice(ctx,s,-h,v,direct=direct)/(D+1+2*h)
    f0,f2=taylor_coefficients(ctx,s,h,v)
    return (positive+negative)/2-f0-f2-bubble_parameter(ctx,v)/(D*D)

def evaluate(ctx,point,direct=False,chart='linear'):
    x,u,v=map(ctx.mpf,point);assert all(0<t<1 for t in (x,u,v))
    assert chart in ('linear','square')
    if chart=='square':jac=2*x;x=x*x
    else:jac=ctx.mpf(1)
    s=x/(1-x);ds=jac/(1-x)**2;theta=ctx.pi*u/2;c=ctx.cos(theta)
    weight=ctx.mpf(8)/3*ctx.sin(theta)**4
    return ds*s*s*bubble(ctx,s)**2/(2*(1+s)**3)*weight*remainder(ctx,s,c,v,direct)
