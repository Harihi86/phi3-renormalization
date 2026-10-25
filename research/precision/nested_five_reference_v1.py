"""Independent complete original routes with formal integrated Taylor jets."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M
import nested_closed_marginals_reference_v2 as V


def density(ctx,p,binding,quadrature=False):
    s,z,theta,x0,y=p;t=(s/(1-s))**2;A=1+t;den=x0*x0+(1-x0)**2
    x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    parent=binding['parent_boundary_shifts'];parent=[*parent,[a-b for a,b in zip(*parent)]]
    tri=binding['triangle_pair_differences'];host=binding['host_routes']
    assert [r[0]**2 for r in tri]==[1,1,0]
    constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    assert constant(tri[2])==1
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    V0=sum(w*r[0]**2*t for w,r in zip(weights,parent));K0=M.bubble_K(ctx,x,V0);T0=-H.edge_log(ctx,t)/2
    total=ctx.mpf(0);radial=ctx.sqrt(t*z)
    for sign in (-1,1):
        qx=sign*radial*ctx.cospi(theta);qy=sign*radial*ctx.sinpi(theta);q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        norm=lambda r:r[0]**2*t+2*r[0]*(r[1]*q1+r[2]*q2)+constant(r)
        V=sum(w*norm(r) for w,r in zip(weights,parent));K=M.bubble_K(ctx,x,V)
        aa,bb=norm(tri[0]),norm(tri[1])
        if quadrature:
            T=-ctx.quad(lambda u:(1-u)*M.triangle_log_marginal(ctx,u,aa,bb,ctx.mpf(1)),[0,ctx.mpf('.5'),1])
        else:T=-H.triangle_log_complete(ctx,aa,bb)
        ratio=ctx.fprod(A/(1+norm(r)) for r in host)
        # Complete subtraction, independent of the native difference rearrangement.
        total+=(ratio*K*T-K0*T0)/2
    common=2/(s*(1-s))*(t/A)**3*(1-x)*jac
    base=common*(1-z)*total
    if quadrature:
        # Independently integrate the frozen six-dimensional formal-series control.
        correction=ctx.quad(lambda u:V.correction(ctx,[*p,u],binding),[0,ctx.mpf('.5'),1])
        return base+correction
    # Derivatives of the fully integrated scalar log give the control moments.
    first=ctx.diff(lambda v:H.edge_log(ctx,v),t)
    second=ctx.diff(lambda v:H.edge_log(ctx,v),t,2)
    U1=-first/2;Cpair=-first/4;Cfixed=-(1-(t+3)*first)/12;I2=-second
    Kp=ctx.diff(lambda v:M.bubble_K(ctx,x,v),V0)
    Kpp=ctx.diff(lambda v:M.bubble_K(ctx,x,v),V0,2)
    def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        rh=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            L,C=linear(r)/A,constant(r)/A;rh=mul(rh,[ctx.mpf(1),-L,L*L-C])
        Lp=sum(w*linear(r) for w,r in zip(weights,parent));Cp=sum(w*constant(r) for w,r in zip(weights,parent))
        kp=[K0,Kp*Lp,Kp*Cp+Kpp*Lp*Lp/2]
        qa,qb=linear(tri[0])/2,linear(tri[1])/2
        tp=[T0,U1*(qa+qb),Cpair*(constant(tri[0])+constant(tri[1]))+Cfixed*constant(tri[2])+I2*(qa*qa+qa*qb+qb*qb)/3]
        return mul(mul(rh,kp),tp)[2]
    point=coefficient(radial*ctx.cospi(theta),radial*ctx.sinpi(theta));r=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((r,0),(-r,0),(0,r),(0,-r)))+8*coefficient(0,0))/12
    return base+common*(mean/2-(1-z)*point)


def evaluate(point,binding,digits=640,mode='five-direct',quadrature=False):
    assert mode in ('five-direct','five-radial')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    if mode=='five-direct':return ctx,density(ctx,p,binding,quadrature)
    s=p[0];x0=p[3];y=p[4];x=x0*x0/(x0*x0+(1-x0)**2)
    a,b=binding['parent_boundary_shifts'];loops=[a[0],b[0],a[0]-b[0]]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    ell=sum(w*l*l for w,l in zip(weights,loops));h=1/ctx.sqrt(ell);den=1-s+h*s
    pp=p.copy();pp[0]=h*s/den
    return ctx,h/(den*den)*density(ctx,pp,binding,quadrature)
