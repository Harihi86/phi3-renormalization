"""Complete normalized six-dimensional spherical average of two propagators.

No expansion in external momentum and no truncation of the angular domain.
The H function is elementary; its small-argument series has alternating terms.
"""
import double_weighted_angular_v1 as G


def h_function(ctx, z):
    z=ctx.mpf(z)
    assert z>=0
    if not z:return ctx.mpf(1)/3
    if z<=ctx.mpf(1)/4:
        total=ctx.mpf(0);power=ctx.mpf(1)
        for n in range(10000):
            term=power/((2*n+1)*(2*n+3));total+=term
            if abs(term)<ctx.eps*max(abs(total),ctx.mpf(1)):return total
            power*=-z
        raise ArithmeticError('H series did not converge')
    root=ctx.sqrt(z)
    return ((1+z)*ctx.atan(root)/root-1)/(2*z)


def pair(ctx,s,t,dot,e=1):
    s,t,dot,e=map(ctx.mpf,(s,t,dot,e))
    assert s>=0 and t>=0 and e>=0 and dot*dot<=s*t
    ks,rs=G.shape(ctx,s,e);kt,rt=G.shape(ctx,t,e)
    q=ks*kt*e*dot
    eta2=(ks*kt*e)**2*(s*t-dot*dot)
    u=1-q;v=rs*rs+rt*rt-2*q
    return ks*kt*(1-v/u*h_function(ctx,eta2/(u*u)))


def feynman_reference(ctx,s,t,dot,e=1):
    s,t,dot,e=map(ctx.mpf,(s,t,dot,e))
    def integrand(x):
        a=1+e+t+x*(s-t)
        norm2=t+2*x*(dot-t)+x*x*(s+t-2*dot)
        return (2/(a+ctx.sqrt(a*a-4*e*norm2)))**2
    return ctx.quad(integrand,[0,ctx.mpf(1)/4,ctx.mpf(1)/2,ctx.mpf(3)/4,1])
