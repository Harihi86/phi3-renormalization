"""Eight-dimensional complete forest, built from exact original-route bindings."""
from fractions import Fraction as Q
import mpmath as mp

def number(ctx,v):return ctx.mpf(v.numerator)/v.denominator if isinstance(v,(Q,int)) else ctx.mpf(v)


def components(point, binding, ctx, sign=1):
    p=[number(ctx,v) for v in point]
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


def one(point,binding,ctx,mode,stable=True,control=False):
    result=ctx.mpf(0)
    for sign in (-1,1):
        x,u,V0,DV,VT0,DVT,logratio,common=components(point,binding,ctx,sign)
        K0,DK=local(ctx,x,u,V0,DV,mode)
        T0=-ctx.log1p(VT0);DT=-ctx.log1p(DVT/(1+VT0));ratio=ctx.exp(logratio)
        if stable:value=ctx.expm1(logratio)*K0*T0+ratio*(DK*T0+K0*DT+DK*DT)
        else:value=ratio*(K0+DK)*(T0+DT)-K0*T0
        result+=common*value/2
    if control:result+=control_term(point,binding,ctx,mode)
    return result


def evaluate(point,binding,mode='primitive',orbit=1,digits=300,stable=True,control=False):
    assert mode in ('primitive','spectral','primitive-triangle','spectral-triangle','primitive-joint','spectral-joint') and orbit in (1,4,8)
    ctx=mp.mp.clone();ctx.dps=digits;result=ctx.mpf(0)
    for j in range(orbit):
        pp=list(point);pp[2]=(pp[2]+Q(j,orbit))%1
        result+=mapped_one(pp,binding,ctx,mode,stable,control)
        pp[5]=1-pp[5];result+=mapped_one(pp,binding,ctx,mode,stable,control)
    return ctx,result/(2*orbit)


def control_term(point,binding,ctx,mode):
    p=[number(ctx,v) for v in point]
    s,z,theta,x0,y,u,tu,tv=p;t=(s/(1-s))**2;A=1+t;den=x0*x0+(1-x0)**2
    x=x0*x0/den;jac=2*x0*(1-x0)/den**2;r=ctx.sqrt(t*z)
    q=[r*ctx.cospi(theta),r*(-ctx.cospi(theta)+ctx.sqrt(3)*ctx.sinpi(theta))/2]
    a,b=binding['parent_boundary_shifts'];rows=[a,b,[i-j for i,j in zip(a,b)]]
    pw=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    rem=ctx.sqrt(1-tu);v=[tu/(1+rem),rem*tv,rem*(1-tv)];tw=[v[0]*v[1],v[0]*v[2],v[1]*v[2]]
    norm=lambda r:r[0]**2+r[1]**2-r[0]*r[1]
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]-(a[0]*b[1]+a[1]*b[0])/2
    scalar=lambda a:sum(i*j for i,j in zip(q,a))
    def vertex(weights,rs):
        V0=t*sum(w*r[0]**2 for w,r in zip(weights,rs))
        C=sum(w*norm(r[1:]) for w,r in zip(weights,rs))
        vector=[2*sum(w*r[0]*r[j] for w,r in zip(weights,rs)) for j in (1,2)]
        return V0,C,vector
    V0,Cp,pv=vertex(pw,rows);VT0,Ct,tv=vertex(tw,binding['triangle_pair_differences'])
    hs=[[2*r[0]*r[1],2*r[0]*r[2]] for r in binding['host_routes']]
    hv=[sum(r[j] for r in hs) for j in range(2)];rv=[-v/A for v in hv]
    hc=sum(norm(r[1:]) for r in binding['host_routes'])
    R1=scalar(rv);R2=(scalar(hv)**2+sum(scalar(r)**2 for r in hs))/(2*A*A)-hc/A
    meanR2=t*(norm(hv)+sum(norm(r) for r in hs))/(12*A*A)-hc/A
    K0=local(ctx,x,u,V0,ctx.mpf(0),mode)[0];arg=(1-x+V0)/x;B=bubble(ctx,arg)
    # Differentiate the independent asinh expression, not the native atanh formula.
    Bprime=ctx.diff(lambda v:bubble(ctx,v),arg)
    Kp=-x*B/(1+V0)**2;Kpp=2*x*B/(1+V0)**3-Bprime/(1+V0)**2
    T0=-ctx.log1p(VT0);Tp=-1/(1+VT0);Tpp=1/(1+VT0)**2;Lp=scalar(pv);Lt=scalar(tv)
    value=R2*K0*T0+R1*(Kp*Lp*T0+K0*Tp*Lt)+(Kp*Cp+Kpp*Lp**2/2)*T0+Kp*Tp*Lp*Lt+K0*(Tp*Ct+Tpp*Lt**2/2)
    mean=meanR2*K0*T0+t/6*(dot(rv,pv)*Kp*T0+dot(rv,tv)*K0*Tp)+(Kp*Cp+Kpp*t*norm(pv)/12)*T0+Kp*Tp*t*dot(pv,tv)/6+K0*(Tp*Ct+Tpp*t*norm(tv)/12)
    common=2/(s*(1-s))*(t/A)**3*(1-x)*jac/2
    return common*(mean/2-(1-z)*value)


def mapped_one(point,binding,ctx,mode,stable=True,control=False):
    if '-' not in mode:return one(point,binding,ctx,mode,stable,control)
    base,kind=mode.split('-');p=[number(ctx,v) for v in point]
    t=(p[0]/(1-p[0]))**2;ell=ctx.log1p(t/2)
    def triangle(pp):
        xi=ctx.expm1(pp[6]*ell)/t;jac=ell*ctx.exp(pp[6]*ell)/t;result=ctx.mpf(0)
        for side in (0,1):
            first=xi if side==0 else 1-xi;qq=list(pp);qq[6]=xi*(2-xi) if side==0 else 1-xi*xi
            result+=2*(1-first)*jac*one(qq,binding,ctx,base,stable,control)
        return result
    if kind=='triangle':return triangle(p)
    assert kind=='joint';q=[0]+[r[0] for r in binding['parent_boundary_shifts']]
    index=next(i for i,value in enumerate(q) if q.count(value)==1)
    xi=ctx.expm1(p[3]*ell)/t;jac=ell*ctx.exp(p[3]*ell)/t;result=ctx.mpf(0)
    for first in (xi,1-xi):
        bary=[None]*3;bary[index]=first;other=[j for j in range(3) if j!=index]
        bary[other[0]]=(1-first)*p[4];bary[other[1]]=(1-first)*(1-p[4])
        x=bary[0];old=ctx.sqrt(x)/(ctx.sqrt(x)+ctx.sqrt(1-x));den=old*old+(1-old)**2
        oldjac=2*old*(1-old)/den**2;pp=list(p);pp[3]=old;pp[4]=bary[1]/(1-x)
        result+=jac*(1-first)/((1-x)*oldjac)*triangle(pp)
    return result
