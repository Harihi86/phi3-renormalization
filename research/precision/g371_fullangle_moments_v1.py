"""Normalized six-dimensional relative-angle inverse moments, without cancellation."""
def moments(ctx,a,b):
    if not a>abs(b):raise ValueError('strict massive relative-angle domain required')
    r=ctx.sqrt((a-b)*(a+b));s=a+r
    return (4*(a+2*r)/(3*s*s),4/(s*s),4/(r*s*s),4*(a+2*r)/(3*r**3*s*s))


def integral(ctx,fun,a,b):
    epsilon=min(ctx.mpf('.25'),ctx.sqrt((a-abs(b))/a))
    nodes=sorted(set([ctx.mpf(0),epsilon,ctx.pi/2,ctx.pi-epsilon,ctx.pi]))
    return 8/(3*ctx.pi)*ctx.quad(lambda theta:ctx.sin(theta)**4*fun(ctx.cos(theta)),nodes)
