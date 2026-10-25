"""Complete equal-mass triangle log integral for one unit external invariant.

All quantities are evaluated in the caller's independent mpmath context.
The real quadratic roots and off-cut dilogarithms implement an exact scalar
identity; no integral sampling occurs in triangle_log_complete.
"""


def edge_log(ctx,t):
    if not t:return ctx.mpf(0)
    return 2*ctx.sqrt((t+4)/t)*ctx.asinh(ctx.sqrt(t)/2)-2


def log_linear_pair(ctx,a,b):
    # Integral from 0 to 1 of b*log(1+a*v)/(1+b*v).
    # a is real, 1+a>0, Im(b) != 0. The dilogarithm segment stays off its cut.
    assert a.imag==0 and 1+a.real>0 and b.imag!=0
    if not a:return ctx.mpc(0)
    z=a/(a-b)
    return ctx.log(1-a/b)*ctx.log(1+b)-ctx.polylog(2,z*(1+b))+ctx.polylog(2,z)


def inverse_triangle(ctx,a,b):
    alpha=(a-b)**2-2*(a+b)+1;beta=2*(a+b-1);r=ctx.sqrt(5)
    assert a>0 and b>0 and alpha<0
    kappa=ctx.sqrt(-alpha);z0=beta/(2*r);q1=ctx.sqrt((a-b)**2+4)
    delta=(4*r*r*alpha-beta*beta)/(2*r*(2*r*q1+2*r*r+beta))
    assert delta<0
    kernel=delta/(z0-ctx.j*kappa)
    result=2*ctx.log((r+1)/(r-1))*ctx.log(1+kernel)
    for linear,square in ((a-b-1,b),(-a+b-1,a)):
        for sign in (1,-1):
            leading=sign-r;middle=beta-2*r*sign*linear
            discriminant=4*ctx.sqrt(square*(square+4))
            for root in ((-middle+discriminant)/(2*leading),(-middle-discriminant)/(2*leading)):
                factor=delta/(z0-root)
                assert 1+factor>0
                result+=sign*log_linear_pair(ctx,factor,kernel)
    return -2*ctx.im(result)/kappa


def triangle_log_complete(ctx,a,b):
    # F=1+a*x+b*y-a*x^2-b*y^2+(1-a-b)*x*y, triangle x,y>=0,x+y<=1.
    cross=(a+b-1)/2;det=a*b-cross*cross
    assert det>0
    cx=b*(a-cross)/(2*det);cy=a*(b-cross)/(2*det)
    height=1+(a*cx+b*cy)/2
    boundary=cx*edge_log(ctx,b)+cy*edge_log(ctx,a)+(1-cx-cy)*edge_log(ctx,ctx.mpf(1))
    return boundary/2-ctx.mpf('.5')+height*inverse_triangle(ctx,a,b)


def integrated_control_moments(ctx,t):
    j0=4*ctx.asinh(ctx.sqrt(t)/2)/ctx.sqrt(t*(t+4))
    j2=2*(1+j0)/(t+4)
    return {'T0':-edge_log(ctx,t)/2,'U1':-(1-j0)/(2*t),
            'Cpair':-(1-j0)/(4*t),'Cfixed':-(j0-3*(1-j0)/t)/12,
            'I2':(1-2*j0+j2)/(t*t)}
