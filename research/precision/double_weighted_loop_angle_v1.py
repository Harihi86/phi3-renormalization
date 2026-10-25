"""G20: average the q direction first, then cancel low angular modes by an exact orbit."""
import double_weighted_full_angle_v1 as F
G=F.G


def mass_kinetic(ctx,s,t):
    ds,dt=1+s,1+t;k,r=G.shape(ctx,s,t)
    j1=k*(1-r*r/3);j2=k*k;j3=k**3/(1-r*r)
    a=ctx.mpf(2)/3*(-1/ds-1/dt+s/ds**2+t/dt**2)-(1+s+t)/(3*ds*dt)
    b=-ctx.mpf(1)/3+(s-t-1)/(3*ds)+(t-s-1)/(3*dt)
    mass=j1/(ds*dt);kinetic=(a*j1+b*j2-ctx.mpf(2)*j3/3+1/(3*ds*dt))/(ds*dt)
    return mass,kinetic


def parents(ctx,s,t):
    js=-1/(1+s)**2+ctx.mpf(2)*s/(3*(1+s)**3)
    jt=-1/(1+t)**2+ctx.mpf(2)*t/(3*(1+t)**3)
    return (G.coefficient(ctx,s,ctx.mpf(1),0)-1/(1+s)-js)/(1+t)**2+(G.coefficient(ctx,t,ctx.mpf(1),0)-1/(1+t)-jt)/(1+s)**2


def bare(ctx,s,t,b,terms=1024,chebyshev=False):
    dot=1+ctx.sqrt(t)*b;u=t+1+2*ctx.sqrt(t)*b
    assert u>0
    return F.pair(ctx,ctx.mpf(1),u,dot,e=s,terms=terms,chebyshev=chebyshev)/(1+u)


def angular_orbit(ctx,s,t,u,count=8,terms=1024,chebyshev=False):
    mass,kinetic=mass_kinetic(ctx,s,t);sub=mass+kinetic+parents(ctx,s,t);value=ctx.mpf(0)
    for i in range(count):
        fraction=u/2+ctx.mpf(i)/count
        if fraction==1:continue
        b=ctx.cospi(fraction);weight=ctx.mpf(8)/3*ctx.sinpi(fraction)**4
        value+=weight*((bare(ctx,s,t,b,terms,chebyshev)+bare(ctx,s,t,-b,terms,chebyshev))/2-sub)
    return value/count


def evaluate(ctx,point,chart='eight',terms=1024,chebyshev=False):
    x,y,u=map(ctx.mpf,point);assert all(0<v<1 for v in (x,y,u))
    s=x/(1-x)**3;t=y/(1-y)**3;ds=(1+2*x)/(1-x)**4;dt=(1+2*y)/(1-y)**4
    count={'four':4,'eight':8}[chart]
    return ds*dt*s*s*t*t*G.B.bubble(ctx,s)*G.B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)*angular_orbit(ctx,s,t,u,count,terms,chebyshev)
