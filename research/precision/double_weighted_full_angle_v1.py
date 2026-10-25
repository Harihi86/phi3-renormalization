"""G20: full external-direction average via exact two-propagator spherical convolutions."""
import double_weighted_angular_v1 as G


def pair(ctx,s,t,dot,e=1,terms=1024,derivative=False,chebyshev=False):
    e=ctx.mpf(e);assert s>0 and t>0 and e>0 and dot*dot<=s*t
    ks,rs=G.shape(ctx,s,e);kt,rt=G.shape(ctx,t,e)
    q=ks*kt*e*dot;z2=(ks*kt*e)**2*s*t;aa,bb=rs*rs,rt*rt
    def logderivative(radius):
        a=1+radius+e;delta=ctx.sqrt(a*a-4*radius*e)
        return -(1+(a-2*radius)/delta)/(a+delta)
    gs,gt=logderivative(s),logderivative(t);g=gs+gt
    old2=ctx.mpf(1);old1=4*q;uold2=ctx.mpf(1);uold1=2*q;total=ctx.mpf(0)
    for n in range(terms):
        if n==0:poly=old2
        elif n==1:poly=old1
        elif chebyshev:
            unew=2*q*uold1-z2*uold2;poly=(n+1)*unew+z2*old2
            uold2,uold1=uold1,unew;old2,old1=old1,poly
        else:
            poly=(ctx.mpf(2)*(n+1)*q*old1-(n+2)*z2*old2)/n;old2,old1=old1,poly
        da=ctx.mpf(1)/(n+1)-aa/(n+3);db=ctx.mpf(1)/(n+1)-bb/(n+3)
        if derivative:
            dap=-2*aa*(gs+ctx.mpf(1)/(2*e))/(n+3);dbp=-2*bb*(gt+ctx.mpf(1)/(2*e))/(n+3)
            factor=((n+1)*g+n/e)*da*db+dap*db+da*dbp
        else:factor=da*db
        total+=ctx.mpf(2)/(n+2)*poly*factor
    return ks*kt*total


def triangle(ctx,s,t,h,e=1,terms=1024,chebyshev=False):
    e=ctx.mpf(e);third=s+t-2*h;den=1+e+2*h
    assert third>0 and h*h<s*t
    derivative=not den
    args={'e':e,'terms':terms,'derivative':derivative,'chebyshev':chebyshev}
    numerator=pair(ctx,s,third,s-h,**args)+pair(ctx,t,third,t-h,**args)-pair(ctx,s,t,-h,**args)
    return numerator if derivative else numerator/den


def mass_kinetic(ctx,s,t,h):
    third=s+t-2*h;ds,dt,dl=1+s,1+t,1+third
    mass=1/(ds*dt*dl)
    kinetic=mass*(-1/ds-1/dt-1/dl+ctx.mpf(2)/3*(s/ds**2+t/dt**2+third/dl**2-h/(ds*dt)+(s-h)/(ds*dl)+(t-h)/(dt*dl)))
    return mass,kinetic


def remainder(ctx,s,t,h,terms=1024,chebyshev=False):
    mass,kinetic=mass_kinetic(ctx,s,t,h)
    js=-1/(1+s)**2+ctx.mpf(2)*s/(3*(1+s)**3)
    jt=-1/(1+t)**2+ctx.mpf(2)*t/(3*(1+t)**3)
    ps=(G.coefficient(ctx,s,ctx.mpf(1),0)-1/(1+s)-js)/(1+t)**2
    pt=(G.coefficient(ctx,t,ctx.mpf(1),0)-1/(1+t)-jt)/(1+s)**2
    return triangle(ctx,s,t,h,terms=terms,chebyshev=chebyshev)-mass-kinetic-ps-pt


def evaluate(ctx,point,chart='cubic',terms=1024,chebyshev=False):
    x,y,u=map(ctx.mpf,point);assert all(0<v<1 for v in (x,y,u))
    power={'linear':1,'cubic':3}[chart]
    convert=lambda v:(v/(1-v)**power,(1+(power-1)*v)/(1-v)**(power+1))
    s,ds=convert(x);t,dt=convert(y);angle=ctx.pi*u;h=ctx.sqrt(s*t)*ctx.cos(angle)
    weight=ctx.mpf(8)/3*ctx.sin(angle)**4
    return ds*dt*s*s*t*t*G.B.bubble(ctx,s)*G.B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)*weight*remainder(ctx,s,t,h,terms,chebyshev)
