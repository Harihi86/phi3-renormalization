"""Independent six-cube reference from original routed quadratic forms."""
import mpmath as mp
import nested_analytic_marginals_v1 as M


def direct(ctx,p,binding,quadrature=False):
    s,z,theta,x0,y,u=p;t=(s/(1-s))**2;A=1+t
    den=x0*x0+(1-x0)**2;x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    parent=binding['parent_boundary_shifts'];parent=[*parent,[a-b for a,b in zip(*parent)]]
    tri=binding['triangle_pair_differences'];weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    v0=sum(w*r[0]**2*t for w,r in zip(weights,parent));tri0=[r[0]**2*t for r in tri]
    def K(V):
        if not quadrature:return M.bubble_K(ctx,x,V)
        a=(1-x)/x;b=(1-x+V)/x
        def f(t):
            B=(t+4)**ctx.mpf('1.5')/(3*ctx.sqrt(t))*ctx.asinh(ctx.sqrt(t)/2)-4*t/9-ctx.mpf(4)/3
            return -B/(1+t)**2
        return ctx.quad(f,[a,(a+b)/2,b])
    def triangle(rows):
        if not quadrature:return -M.triangle_log_marginal(ctx,u,*rows)
        a,b,c=rows
        return -ctx.quad(lambda v:ctx.log1p(u*(1-u)*(v*a+(1-v)*b)+(1-u)**2*v*(1-v)*c),[0,ctx.mpf('.5'),1])
    K0,T0=K(v0),triangle(tri0);value=ctx.mpf(0)
    for sign in (-1,1):
        radial=sign*ctx.sqrt(t*z);q1=radial*ctx.cospi(theta);q2=radial*(-ctx.cospi(theta)+ctx.sqrt(3)*ctx.sinpi(theta))/2
        def norm(r):
            a,b,c=r
            return a*a*t+2*a*(b*q1+c*q2)+b*b+c*c-b*c
        V=sum(w*norm(r) for w,r in zip(weights,parent));tri1=[norm(r) for r in tri]
        ratio=ctx.fprod(A/(1+norm(r)) for r in binding['host_routes'])
        value+=(ratio*K(V)*triangle(tri1)-K0*T0)/2
    return 2/(s*(1-s))*(t/A)**3*(1-z)*(1-x)*jac*(1-u)*value


def evaluate(point,binding,digits=640,mode='direct',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    if mode=='direct':return ctx,direct(ctx,p,binding,quadrature)
    assert mode=='joint'
    index=0 if binding['parent_boundary_shifts'][0][0]==binding['parent_boundary_shifts'][1][0] else 1
    t=(p[0]/(1-p[0]))**2;ell=ctx.log1p(t/2);xi=ctx.expm1(p[3]*ell)/t;dx=ell*ctx.exp(p[3]*ell)/t
    eta=ctx.expm1(p[5]*ell)/t;du=ell*ctx.exp(p[5]*ell)/t;total=ctx.mpf(0)
    for first in (xi,1-xi):
        bary=[None]*3;bary[index]=first
        remaining=[j for j in range(3) if j!=index]
        bary[remaining[0]]=(1-first)*p[4];bary[remaining[1]]=(1-first)*(1-p[4])
        x=bary[0];old=ctx.sqrt(x)/(ctx.sqrt(x)+ctx.sqrt(1-x));den=old*old+(1-old)**2;oldjac=2*old*(1-old)/den**2
        for last in (eta,1-eta):
            pp=p.copy();pp[3]=old;pp[4]=bary[1]/(1-x);pp[5]=last
            total+=dx*(1-first)/((1-x)*oldjac)*du*direct(ctx,pp,binding,quadrature)
    return ctx,total
