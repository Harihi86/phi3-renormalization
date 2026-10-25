"""G20 complete remaining angular average with an explicit Cauchy strip error."""
import spherical_pair_closed_v1 as P
import double_weighted_loop_angle_v2 as A


def bare(ctx,s,t,b,direct=False):
    u=t+1+2*ctx.sqrt(t)*b;dot=1+ctx.sqrt(t)*b
    if direct:return P.feynman_reference(ctx,ctx.mpf(1),u,dot,e=s)/(1+u)
    ks,rs=P.G.shape(ctx,ctx.mpf(1),s);kt,rt=P.G.shape(ctx,u,s)
    pref=ks*kt;q=pref*s*dot;eta2=(pref*s)**2*t*(1-b*b);den=1-q
    return pref*(1-(rs*rs+rt*rt-2*q)/den*P.h_function(ctx,eta2/(den*den)))/(1+u)


def cauchy_bound(ctx,s,t,count):
    ds=2+s-ctx.mpf(5)*ctx.sqrt(s)/2;dt=2+t-ctx.mpf(5)*ctx.sqrt(t)/2
    return ctx.mpf(1250)/(42*ds*dt*(ctx.mpf(4)**count-1))


def trapezoid(ctx,s,t,count,direct=False):
    value=ctx.mpf(0)
    for j in range(count):
        f=(ctx.mpf(j)+ctx.mpf(1)/2)/count;b=ctx.cospi(f);weight=ctx.mpf(8)/3*ctx.sinpi(f)**4
        value+=weight*bare(ctx,s,t,b,direct)
    return value/count


def angular(ctx,s,t,relative_digits=100):
    mass,kinetic=A.mass_kinetic(ctx,s,t);sub=mass+kinetic+A.parents(ctx,s,t)
    for count in (32,64,128,256,512,1024,2048):
        value=trapezoid(ctx,s,t,count)-sub;bound=cauchy_bound(ctx,s,t,count)
        if bound<ctx.mpf(10)**(-relative_digits)*abs(value):return value,bound,count
    raise ArithmeticError('Complete angular Cauchy enclosure did not meet the fixed reference tolerance')


def evaluate(ctx,point,chart='cubic',relative_digits=100):
    x,y=map(ctx.mpf,point);assert all(0<v<1 for v in (x,y))
    convert=lambda v:(v/(1-v)**3,(1+2*v)/(1-v)**4)
    a,da=convert(x);b,db=convert(y)
    if chart=='cubic':
        value,bound,count=angular(ctx,a,b,relative_digits)
        pref=da*db*a*a*b*b*A.G.B.bubble(ctx,a)*A.G.B.bubble(ctx,b)/(4*(1+a)**2*(1+b)**2)
        return pref*value,abs(pref)*bound,count
    assert chart=='diagonal'
    s,t=(a+b)**2,a*a;v1,e1,n1=angular(ctx,s,t,relative_digits);v2,e2,n2=angular(ctx,t,s,relative_digits)
    pref=4*a*(a+b)*da*db*s*s*t*t*A.G.B.bubble(ctx,s)*A.G.B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)
    return pref*(v1+v2),abs(pref)*(e1+e2),max(n1,n2)
