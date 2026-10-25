"""Closed scalar marginals; original graphs and remaining integrations are separate."""
import mpmath as mp


def atan_log_integral(ctx,q):
    """Integral_0^q atan(sqrt(3)*s)/(1-s^2) ds, 0<=q<1."""
    assert 0<=q<1
    a=ctx.j*ctx.sqrt(3)
    def piece(b):
        z=-a/(b-a);u=1+b*q
        return (ctx.log(1-a/b)*ctx.log(u)-ctx.polylog(2,z*u)+ctx.polylog(2,z))/b
    return ctx.im(piece(-1)+piece(1))/2


def bubble_primitive(ctx,t):
    """A'(t)=B(t)/(1+t)^2, A(0)=0. No numerical inner quadrature."""
    assert t>=0
    if not t:return ctx.mpf(0)
    q=ctx.sqrt(t/(t+4));h=ctx.asinh(ctx.sqrt(t)/2);angle=ctx.atan(ctx.sqrt(3)*q)
    return (h*h/3+2*ctx.sqrt(3)*(h*angle-atan_log_integral(ctx,q))
            +4*q*h/(1+3*q*q)-ctx.mpf(17)/18*ctx.log1p(t)-ctx.mpf(32)/9*q*q/(1+3*q*q))


def bubble_K(ctx,x,V):
    assert 0<x<1 and V>=0
    return -(bubble_primitive(ctx,(1-x+V)/x)-bubble_primitive(ctx,(1-x)/x))


def linear_log_mean(ctx,z):
    assert z>-1
    return ((1+z)*ctx.log1p(z)-z)/z if z else ctx.mpf(0)


def triangle_log_marginal(ctx,u,s01,s02,s12):
    """Integral_0^1 log(1+V(u,v)) dv; simplex measure (1-u) is external."""
    assert 0<=u<=1 and min(s01,s02,s12)>=0
    m=1+u*(1-u)*s02;end=1+u*(1-u)*s01
    C=(1-u)**2*s12;B=end-m+C
    if not C:return ctx.log(m)+linear_log_mean(ctx,B/m)
    root=ctx.sqrt(B*B+4*C*m)
    alpha=2*m/(root+B) if B>=0 else (root-B)/(2*C)
    tau=2*end/(root-B+2*C)
    # q(v)=C*(v+alpha)*(1+tau-v), both roots outside [0,1].
    return ctx.log(m)+linear_log_mean(ctx,1/alpha)+tau*ctx.log1p(1/tau)-1
