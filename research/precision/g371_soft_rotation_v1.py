"""SO(6) rotation moments retaining a complete shifted soft propagator."""
def moments(ctx,w):
    if not w>=0:raise ValueError('nonnegative massive soft radius required')
    a=w+2;b=2*ctx.sqrt(w);r=ctx.sqrt(w*w+4);s=a+r
    m=4*(a+2*r)/(3*s*s)
    mu1=-b*(a+3*r)/(3*s**3)
    mu2=a*(a+3*r)/(3*s**3)
    return m,mu1,mu2


def weighted_projection_moments(ctx,t,v,h):
    """Return E[1/Y], E[projection/Y], E[projection*projection^T/Y]."""
    w=t+v+2*h
    if not t>0 or not v>0 or not t*v>h*h or not w>0:
        raise ValueError('strict two-vector Gram domain required')
    m,mu1,mu2=moments(ctx,w)
    gram=[[ctx.mpf(1),-ctx.mpf('.5')],[-ctx.mpf('.5'),ctx.mpf(1)]]
    projected=[[mu2/4+(gram[a][b]-ctx.mpf('.25'))*(m-mu2)/5 for b in range(2)] for a in range(2)]
    alpha=[(t+h)/w,(v+h)/w];rho2=(t*v-h*h)/w;sgn=[1,-1]
    first=[alpha[i]*ctx.sqrt(w)*mu1/2 for i in range(2) for a in range(2)]
    second=[[alpha[i]*alpha[j]*w*projected[a][b]+sgn[i]*sgn[j]*rho2/5*(gram[a][b]*m-projected[a][b])
             for j in range(2) for b in range(2)] for i in range(2) for a in range(2)]
    return m,first,second
