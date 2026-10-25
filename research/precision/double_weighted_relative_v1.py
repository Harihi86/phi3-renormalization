"""G20/G9: exact transverse-angle integral and actual external Taylor coefficients."""
import shared_overlap_radial_v1 as B


def transverse(ctx,centre,width):
    width=abs(width);assert centre>width
    if not width:return 1/centre,-1/centre**2,2/centre**3
    h=ctx.atanh(width/centre)
    value=ctx.mpf(3)/(2*width**2)*(centre-(centre**2-width**2)*h/width)
    first=ctx.mpf(3)/width**2*(1-centre*h/width)
    second=ctx.mpf(3)/width**2*(centre/(centre**2-width**2)-h/width)
    return value,first,second


def remainder(ctx,s,t,a,b,gid='L4E2G20'):
    assert gid in ('L4E2G20','L4E2G9')
    ds,dt=1+s,1+t;rs,rt=ctx.sqrt(s),ctx.sqrt(t)
    la,lb=-2*rs*a,2*rt*b;width=2*rs*rt*ctx.sqrt((1-a*a)*(1-b*b))
    if gid=='L4E2G20':centre=1+s+t-2*rs*rt*a*b;linear=la+lb;quadratic=ctx.mpf(1)
    else:centre=1+s+t+2*rs*rt*a*b;linear=ctx.mpf(0);quadratic=ctx.mpf(0)
    f0,f1,f2=transverse(ctx,centre,width)
    g0=1/(ds*dt);g1=-(la/ds+lb/dt)/(ds*dt)
    g2=(la*la/ds**2-1/ds+lb*lb/dt**2-1/dt+la*lb/(ds*dt))/(ds*dt)
    mass=g0*f0;kinetic=g2*f0+g1*linear*f1+g0*(quadratic*f1+linear*linear*f2/2)
    full=sum(transverse(ctx,centre+sign*linear+quadratic,width)[0]/((ds+sign*la+1)*(dt+sign*lb+1)) for sign in (-1,1))/2
    ps=(sum(1/(ds+sign*la+1) for sign in (-1,1))/2-1/ds-(la*la/ds**3-1/ds**2))/dt**2
    pt=(sum(1/(dt+sign*lb+1) for sign in (-1,1))/2-1/dt-(lb*lb/dt**3-1/dt**2))/ds**2
    return full-mass-kinetic-ps-pt


def routed_reference(ctx,s,t,a,b,paired_states):
    """Integrate each exact routed rational forest after its own epsilon Taylor subtraction."""
    rs,rt=ctx.sqrt(s),ctx.sqrt(t);perp=ctx.sqrt((1-a*a)*(1-b*b))
    def integrand(u):
        qk=rs*rt*(a*b+perp*u);result=ctx.mpf(0)
        for state in paired_states:
            rows=[r for i,r in enumerate(state['routes']) if i not in state['weighted_edges']]
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)];plus=minus=ctx.mpf(1)
            for alpha,beta,gamma in rows:
                d0=1+alpha*alpha*s+beta*beta*t+2*alpha*beta*qk
                d1=2*gamma*(alpha*rs*a+beta*rt*b);d2=gamma*gamma
                inverse=[1/d0,-d1/d0**2,d1*d1/d0**3-d2/d0**2]
                jet=[sum(jet[j]*inverse[k-j] for j in range(k+1)) for k in range(3)]
                plus/=d0+d1+d2;minus/=d0-d1+d2
            result+=state['sign']*((plus+minus)/2-jet[0]-jet[2])
        return ctx.mpf(3)/4*(1-u*u)*result
    return ctx.quad(integrand,[-1,0,1])


def evaluate(ctx,point,gid='L4E2G20',chart='four',paired_states=None):
    x,y,u,v=map(ctx.mpf,point);assert all(0<w<1 for w in (x,y,u,v))
    count={'four':4,'eight':8}[chart];s=x/(1-x)**3;t=y/(1-y)**3
    ds=(1+2*x)/(1-x)**4;dt=(1+2*y)/(1-y)**4;total=ctx.mpf(0)
    for ia in range(count):
        fa=u+ctx.mpf(2*ia)/count
        if fa in (1,2):continue
        aa=ctx.pi*fa;a=ctx.cos(aa);wa=ctx.mpf(8)/3*ctx.sin(aa)**4
        for ib in range(count):
            fb=v+ctx.mpf(2*ib)/count
            if fb in (1,2):continue
            bb=ctx.pi*fb;b=ctx.cos(bb);wb=ctx.mpf(8)/3*ctx.sin(bb)**4
            value=remainder(ctx,s,t,a,b,gid) if paired_states is None else routed_reference(ctx,s,t,a,b,paired_states)
            total+=wa*wb*value
    return ds*dt*s*s*t*t*B.bubble(ctx,s)*B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)*total/(count*count)
