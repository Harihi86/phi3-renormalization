"""G8: complete spherical pair, original shifted weight and full T2 forest."""
import spherical_pair_closed_v1 as P
import dynamic_weighted_triangle_v1 as W


def mass_kinetic(ctx,s,t):
    ds,dt=1+s,1+t;ws=W.weighted(ctx,s);wt,wp,wpp=W.weighted_derivatives(ctx,t)
    k,r=P.G.shape(ctx,s,t);j2=k*k;j3=k**3/(1-r*r);j4=k**4*(3-r*r)/(3*(1-r*r)**3)
    numerator=ws*wt-wt*wt
    mass=(numerator*j2-ws*wt/ds**2)/dt
    kinetic=(ws*(wp+t*wpp/3)*(j2-1/ds**2)+numerator*(-j3/3-2*j4/3)
             +ws*wp/3*((s-t+1)*j3-j2))/dt
    return mass,kinetic


def forest(ctx,s,t,b,epsilon=1,direct=False):
    ep=ctx.mpf(epsilon);kp=ctx.sqrt(t)*b;u=t+ep*ep-2*ep*kp;dot=t-ep*kp
    ws,wt,wu=W.weighted(ctx,s),W.weighted(ctx,t),W.weighted(ctx,u)
    pair=(P.feynman_reference if direct else P.pair)(ctx,t,u,dot,e=s)
    return (pair*(ws*wu-wt*wt)-ws*wu/(1+s)**2)/(1+t)


def orbit(ctx,s,t,u,direct=False):
    total=ctx.mpf(0);count=8
    for i in range(count):
        f=u/2+ctx.mpf(i)/count
        if f==1:continue
        b=ctx.cospi(f);weight=ctx.mpf(8)/3*ctx.sinpi(f)**4
        total+=weight*(forest(ctx,s,t,b,direct=direct)+forest(ctx,s,t,-b,direct=direct))/2
    mass,kinetic=mass_kinetic(ctx,s,t)
    return total/count-mass-kinetic


def evaluate(ctx,point,chart='cubic',direct=False):
    x,y,u=map(ctx.mpf,point);assert all(0<z<1 for z in (x,y,u))
    convert=lambda v:(v/(1-v)**3,(1+2*v)/(1-v)**4)
    a,da=convert(x);b,db=convert(y)
    if chart=='cubic':return da*db*a*a*b*b*orbit(ctx,a,b,u,direct)/4
    assert chart=='diagonal'
    s,t=(a+b)**2,a*a;jac=4*a*(a+b)*da*db
    return jac*s*s*t*t*(orbit(ctx,s,t,u,direct)+orbit(ctx,t,s,u,direct))/4
