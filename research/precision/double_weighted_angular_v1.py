"""G9: both angular integrations closed, with the full original T2 forest retained."""
import shared_overlap_radial_v1 as B


def shape(ctx,s,t):
    a=1+s+t;disc=ctx.sqrt(a*a-4*s*t);k=2/(a+disc);r=k*ctx.sqrt(s*t)
    return k,r


def coefficient(ctx,s,t,n):
    k,r=shape(ctx,s,t)
    return k*r**n*(ctx.mpf(1)/(n+1)-r*r/(n+3))


def angular(ctx,s,t,e=1,series=False,terms=256):
    e=ctx.mpf(e)
    if not e:return coefficient(ctx,s,t,0)/((1+s)*(1+t))
    ks,rs=shape(ctx,s,e);kt,rt=shape(ctx,t,e);kh,rh=shape(ctx,s,t)
    z=rs*rt*rh;a,b,c=rs*rs,rt*rt,rh*rh
    if series:
        total=ctx.mpf(0);zn=ctx.mpf(1)
        for n in range(terms):
            factor=ctx.mpf(2)*(n+1)*(n+3)/(3*(n+2))
            total+=zn*factor*(ctx.mpf(1)/(n+1)-a/(n+3))*(ctx.mpf(1)/(n+1)-b/(n+3))*(ctx.mpf(1)/(n+1)-c/(n+3))
            zn*=z
        return ks*kt*kh*total
    e1=a+b+c;e2=a*b+a*c+b*c;e3=a*b*c;log=-ctx.log1p(-z);dilog=ctx.polylog(2,z)
    one1=log/z;two1=dilog/z;one2=(log-z)/(z*z)
    one3=(log-z-z*z/2)/(z**3);two3=(dilog-z-z*z/4)/(z**3)
    return ks*kt*kh*ctx.mpf(2)/3*(2*two1-(1+e1)*one1+(1+e1+e2+e3)*one2-(e2+e3)*one3-2*e3*two3)


def mass_kinetic(ctx,s,t):
    ds,dt=1+s,1+t
    j=lambda w:-1/(1+w)**2+ctx.mpf(2)*w/(3*(1+w)**3)
    h0=coefficient(ctx,s,t,0)/(ds*dt)
    h1=coefficient(ctx,s,t,0)*(j(s)/dt+j(t)/ds)+ctx.mpf(4)*ctx.sqrt(s*t)*coefficient(ctx,s,t,1)/(9*ds*ds*dt*dt)
    return h0,h1,j(s),j(t)


def remainder(ctx,s,t,series=False,terms=256):
    h0,h1,js,jt=mass_kinetic(ctx,s,t)
    parent_s=(coefficient(ctx,s,ctx.mpf(1),0)-1/(1+s)-js)/(1+t)**2
    parent_t=(coefficient(ctx,t,ctx.mpf(1),0)-1/(1+t)-jt)/(1+s)**2
    return angular(ctx,s,t,series=series,terms=terms)-h0-h1-parent_s-parent_t


def evaluate(ctx,point,chart='linear',series=False,terms=256):
    x,y=map(ctx.mpf,point);assert all(0<v<1 for v in (x,y))
    power={'linear':1,'cubic':3}[chart]
    convert=lambda v:(v/(1-v)**power,(1+(power-1)*v)/(1-v)**(power+1))
    s,ds=convert(x);t,dt=convert(y)
    return ds*dt*s*s*t*t*B.bubble(ctx,s)*B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)*remainder(ctx,s,t,series,terms)
