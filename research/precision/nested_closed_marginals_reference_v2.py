"""Independent formal-series control, exact simplex moments and sphere-axis averaging."""
import mpmath as mp
import nested_closed_marginals_reference_v1 as B
M=B.M


def correction(ctx,p,binding):
    s,z,theta,x0,y,u=p;t=(s/(1-s))**2;A=1+t;den=x0*x0+(1-x0)**2
    x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    parent=binding['parent_boundary_shifts'];parent=[*parent,[a-b for a,b in zip(*parent)]]
    tri=binding['triangle_pair_differences'];host=binding['host_routes']
    assert [r[0]**2 for r in tri]==[1,1,0]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    V0=sum(w*r[0]**2*t for w,r in zip(weights,parent));VT0=u*(1-u)*t;arg=(1-x+V0)/x
    def bubble(t):return (t+4)**ctx.mpf('1.5')/(3*ctx.sqrt(t))*ctx.asinh(ctx.sqrt(t)/2)-4*t/9-ctx.mpf(4)/3
    K0=M.bubble_K(ctx,x,V0);Kp=-x*bubble(arg)/(1+V0)**2
    Kpp=2*x*bubble(arg)/(1+V0)**3-ctx.diff(bubble,arg)/(1+V0)**2
    T0=-ctx.log1p(VT0);Tp=-1/(1+VT0);Tpp=1/(1+VT0)**2
    def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
        rh=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            L,C=linear(r)/A,constant(r)/A;rh=mul(rh,[ctx.mpf(1),-L,L*L-C])
        Lp=sum(w*linear(r) for w,r in zip(weights,parent));Cp=sum(w*constant(r) for w,r in zip(weights,parent))
        kp=[K0,Kp*Lp,Kp*Cp+Kpp*Lp*Lp/2];total=ctx.mpf(0)
        # Simpson is exact here: the physical-lambda coefficient is degree <=2 in v.
        for v,weight in ((ctx.mpf(0),ctx.mpf(1)/6),(ctx.mpf('.5'),ctx.mpf(2)/3),(ctx.mpf(1),ctx.mpf(1)/6)):
            tw=[u*(1-u)*v,u*(1-u)*(1-v),(1-u)**2*v*(1-v)]
            Lt=sum(w*linear(r) for w,r in zip(tw,tri));Ct=sum(w*constant(r) for w,r in zip(tw,tri))
            tp=[T0,Tp*Lt,Tp*Ct+Tpp*Lt*Lt/2];total+=weight*mul(mul(rh,kp),tp)[2]
        return total
    radial=ctx.sqrt(t*z);point=coefficient(radial*ctx.cospi(theta),radial*ctx.sinpi(theta));r=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((r,0),(-r,0),(0,r),(0,-r)))+8*coefficient(0,0))/12
    common=2/(s*(1-s))*(t/A)**3*(1-x)*jac*(1-u)
    return common*(mean/2-(1-z)*point)


def evaluate(point,binding,digits=640,mode='direct-control',quadrature=False):
    control='control' in mode;chart=mode.split('-')[0];assert chart in ('direct','joint')
    if not control:return B.evaluate(point,binding,digits,chart,quadrature)
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    def value(p):return B.direct(ctx,p,binding,quadrature)+correction(ctx,p,binding)
    if chart=='direct':return ctx,value(p)
    index=0 if binding['parent_boundary_shifts'][0][0]==binding['parent_boundary_shifts'][1][0] else 1
    t=(p[0]/(1-p[0]))**2;ell=ctx.log1p(t/2);xi=ctx.expm1(p[3]*ell)/t;dx=ell*ctx.exp(p[3]*ell)/t
    eta=ctx.expm1(p[5]*ell)/t;du=ell*ctx.exp(p[5]*ell)/t;total=ctx.mpf(0)
    for first in (xi,1-xi):
        bary=[None]*3;bary[index]=first;left=[j for j in range(3) if j!=index]
        bary[left[0]]=(1-first)*p[4];bary[left[1]]=(1-first)*(1-p[4]);x=bary[0]
        old=ctx.sqrt(x)/(ctx.sqrt(x)+ctx.sqrt(1-x));den=old*old+(1-old)**2;oldjac=2*old*(1-old)/den**2
        for last in (eta,1-eta):
            pp=p.copy();pp[3]=old;pp[4]=bary[1]/(1-x);pp[5]=last
            total+=dx*(1-first)/((1-x)*oldjac)*du*value(pp)
    return ctx,total
