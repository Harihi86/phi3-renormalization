"""Eight-dimensional complete forest, built from exact original-route bindings."""
from fractions import Fraction as Q
import mpmath as mp


def components(point, binding, ctx, sign=1):
    p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    s,z,theta,x0,y,u,tu,tv=p
    t=(s/(1-s))**2;A=1+t;den=x0*x0+(1-x0)**2
    x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    co,si=ctx.cospi(theta),ctx.sinpi(theta);rad=ctx.sqrt(t*z)
    q1=sign*rad*co;q2=sign*rad*(-co+ctx.sqrt(3)*si)/2
    def norm(row):
        a,b,c=row
        return a*a*t+2*a*(b*q1+c*q2)+b*b+c*c-b*c
    def external(row):
        _,b,c=row;return b*b+c*c-b*c
    a,b=binding['parent_boundary_shifts'];ab=[i-j for i,j in zip(a,b)]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    V0=sum(w*r[0]**2*t for w,r in zip(weights,(a,b,ab)))
    DV=sum(w*(2*r[0]*(r[1]*q1+r[2]*q2)+external(r)) for w,r in zip(weights,(a,b,ab)))
    v0=1-ctx.sqrt(1-tu);v1=ctx.sqrt(1-tu)*tv;v2=ctx.sqrt(1-tu)*(1-tv)
    triweights=[v0*v1,v0*v2,v1*v2]
    VT0=sum(w*r[0]**2*t for w,r in zip(triweights,binding['triangle_pair_differences']))
    DVT=sum(w*(2*r[0]*(r[1]*q1+r[2]*q2)+external(r)) for w,r in zip(triweights,binding['triangle_pair_differences']))
    host_changes=[(2*r[0]*(r[1]*q1+r[2]*q2)+external(r))/A for r in binding['host_routes']]
    logratio=-sum(ctx.log1p(v) for v in host_changes)
    common=2/(s*(1-s))*(t/A)**3*(1-z)*(1-x)*jac/2
    return x,u,V0,DV,VT0,DVT,logratio,common


def bubble(ctx,t):
    return (t+4)**ctx.mpf('1.5')/(3*ctx.sqrt(t))*ctx.asinh(ctx.sqrt(t)/2)-4*t/9-ctx.mpf(4)/3


def local(ctx,x,u,V0,DV,mode):
    if mode=='primitive':
        log0=ctx.log1p(V0);t0=(1-x+ctx.expm1(u*log0))/x
        K0=-log0*bubble(ctx,t0)/(1+t0)
        shift=ctx.log1p(DV/(1+V0));t1=((1+V0)*ctx.exp(u*shift)-x)/x
        DK=-shift*bubble(ctx,t1)/(1+t1)
    else:
        assert mode=='spectral'
        b=1-x+V0;ell=ctx.log1p(b/(4*x));a=x/b*ctx.expm1(u*ell);d=1-a;base=x+a*(1-x)
        weight=(1-4*a)**ctx.mpf('1.5')/6*x/b*ell*ctx.exp(u*ell)
        K0=weight*((1+d)/d**2*ctx.log1p(V0)-x/d*V0/(1+V0)-ctx.log1p(a*V0/base)/(a*d*d))
        DK=weight*((1+d)/d**2*ctx.log1p(DV/(1+V0))-x/d*DV/((1+V0)*(1+V0+DV))
                   -ctx.log1p(a*DV/(base+a*V0))/(a*d*d))
    return K0,DK


def one(point,binding,ctx,mode,stable=True):
    result=ctx.mpf(0)
    for sign in (-1,1):
        x,u,V0,DV,VT0,DVT,logratio,common=components(point,binding,ctx,sign)
        K0,DK=local(ctx,x,u,V0,DV,mode)
        T0=-ctx.log1p(VT0);DT=-ctx.log1p(DVT/(1+VT0));ratio=ctx.exp(logratio)
        if stable:value=ctx.expm1(logratio)*K0*T0+ratio*(DK*T0+K0*DT+DK*DT)
        else:value=ratio*(K0+DK)*(T0+DT)-K0*T0
        result+=common*value/2
    return result


def evaluate(point,binding,mode='primitive',orbit=1,digits=300,stable=True):
    assert mode in ('primitive','spectral') and orbit in (1,4,8)
    ctx=mp.mp.clone();ctx.dps=digits;result=ctx.mpf(0)
    for j in range(orbit):
        pp=list(point);pp[2]=(pp[2]+Q(j,orbit))%1
        result+=one(pp,binding,ctx,mode,stable)
        pp[5]=1-pp[5];result+=one(pp,binding,ctx,mode,stable)
    return ctx,result/(2*orbit)
