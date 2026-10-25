"""Independent original routing, scalar quadratures and formal Taylor products."""
import mpmath as mp
import nested_analytic_marginals_v1 as M


def bubble(ctx,t):
    if not t:return ctx.mpf(0)
    return (t+4)**ctx.mpf('1.5')/(3*ctx.sqrt(t))*ctx.asinh(ctx.sqrt(t)/2)-4*t/9-ctx.mpf(4)/3


def density(ctx,p,binding,control=True,quadrature=False):
    s,z,theta,x0,y=p;t=(s/(1-s))**2;A=1+t
    den=x0*x0+(1-x0)**2;x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    parent=binding['parent_boundary_shifts'];parent=[*parent,[a-b for a,b in zip(*parent)]]
    host=binding['host_distinct_routes'];br=binding['disjoint_bubble_boundary']
    constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    V0=sum(w*r[0]**2*t for w,r in zip(weights,parent))
    def parent_value(V):
        if quadrature:return -ctx.quad(lambda u:V*x*bubble(ctx,(1-x+V*u)/x)/(1+V*u)**2,[0,ctx.mpf('.5'),1])
        return M.bubble_K(ctx,x,V)
    def J(value):
        if quadrature:
            numerator=ctx.quad(lambda u:(1+value*u*(1-u))*ctx.log1p(value*u*(1-u))-value*u*(1-u),[0,ctx.mpf('.5'),1])
        else:numerator=bubble(ctx,value)
        return numerator/(1+value)
    K0=parent_value(V0);J0=J(t);radial=ctx.sqrt(t*z);total=ctx.mpf(0)
    for sign in (-1,1):
        qx=sign*radial*ctx.cospi(theta);qy=sign*radial*ctx.sinpi(theta)
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        norm=lambda r:r[0]**2*t+2*r[0]*(r[1]*q1+r[2]*q2)+constant(r)
        V=sum(w*norm(r) for w,r in zip(weights,parent));K=parent_value(V)
        ratio=ctx.fprod(A/(1+norm(r)) for r in host)
        total+=(ratio*K*J(norm(br))-K0*J0)/2
    common=2/(s*(1-s))*(t/A)**3*(1-x)*jac
    result=common*(1-z)*total
    if not control:return result
    arg=(1-x+V0)/x;B=bubble(ctx,arg)
    Kp=-x*B/(1+V0)**2
    Kpp=2*x*B/(1+V0)**3-ctx.diff(lambda value:bubble(ctx,value),arg)/(1+V0)**2
    # Differentiate the full independent bubble quotient; do not reuse native B''.
    jfun=lambda value:bubble(ctx,value)/(1+value)
    Jp=ctx.diff(jfun,t);Jpp=ctx.diff(jfun,t,2)
    def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            l,c=linear(r)/A,constant(r)/A;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
        lp=sum(w*linear(r) for w,r in zip(weights,parent));cp=sum(w*constant(r) for w,r in zip(weights,parent))
        jet=mul(jet,[K0,Kp*lp,Kp*cp+Kpp*lp*lp/2])
        lb=linear(br)
        jet=mul(jet,[J0,Jp*lb,Jp*constant(br)+Jpp*lb*lb/2])
        return jet[2]
    point=coefficient(radial*ctx.cospi(theta),radial*ctx.sinpi(theta));r=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((r,0),(-r,0),(0,r),(0,-r)))+8*coefficient(0,0))/12
    return result+common*(mean/2-(1-z)*point)


def evaluate(point,binding,digits=640,mode='bubble-direct',quadrature=False):
    assert mode in ('bubble-plain','bubble-direct','bubble-radial') and len(point)==5
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    if mode!='bubble-radial':return ctx,density(ctx,p,binding,mode!='bubble-plain',quadrature)
    s=p[0];x0=p[3];y=p[4];x=x0*x0/(x0*x0+(1-x0)**2)
    a,b=binding['parent_boundary_shifts'];loops=[a[0],b[0],a[0]-b[0]]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    ell=sum(w*l*l for w,l in zip(weights,loops));h=1/ctx.sqrt(ell);den=1-s+h*s
    pp=p.copy();pp[0]=h*s/den
    return ctx,h/(den*den)*density(ctx,pp,binding,True,quadrature)
