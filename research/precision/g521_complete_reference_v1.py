"""Independent full routed G521 product and formal Taylor convolution."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='three-control',quadrature=False):
    assert mode in ('three-direct','three-control') and binding['id']=='L4E3G521'
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    s,z,theta=p;t=(s/(1-s))**2;A=1+t;host=binding['host_routes'];children=binding['child_pair_differences']
    T0=-H.edge_log(ctx,t)/2;rad=ctx.sqrt(t*z);total=ctx.mpf(0)
    constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    for sign in (-1,1):
        qx=sign*rad*ctx.cospi(theta);qy=sign*rad*ctx.sinpi(theta);q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        norm=lambda r:r[0]**2*t+2*r[0]*(r[1]*q1+r[2]*q2)+constant(r)
        factors=[]
        for rows in children:
            assert [r[0]**2 for r in rows]==[1,1,0] and constant(rows[2])==1
            a,b=norm(rows[0]),norm(rows[1])
            if quadrature:factors.append(-ctx.quad(lambda u:(1-u)*M.triangle_log_marginal(ctx,u,a,b,ctx.mpf(1)),[0,ctx.mpf('.5'),1]))
            else:factors.append(-H.triangle_log_complete(ctx,a,b))
        ratio=ctx.fprod(A/(1+norm(r)) for r in host)
        total+=(ratio*ctx.fprod(factors)-T0**3)/2
    common=2/(s*(1-s))*(t/A)**3;value=common*(1-z)*total
    if mode=='three-direct':return ctx,value
    if quadrature:
        br=[0,ctx.mpf('.5'),1];ell=lambda u:u*(1-u)
        U1=ctx.quad(lambda u:-(1-u)*ell(u)/(1+t*ell(u)),br)
        Cpair=ctx.quad(lambda u:-(1-u)*ell(u)/(2*(1+t*ell(u))),br)
        Cfixed=ctx.quad(lambda u:-(1-u)**3/(6*(1+t*ell(u))),br)
        I2=ctx.quad(lambda u:ell(u)**2/(1+t*ell(u))**2,br)
    else:
        d1=ctx.diff(lambda v:H.edge_log(ctx,v),t);d2=ctx.diff(lambda v:H.edge_log(ctx,v),t,2)
        U1=-d1/2;Cpair=-d1/4;Cfixed=-(1-(t+3)*d1)/12;I2=-d2
    def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2;linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            a,b=linear(r)/A,constant(r)/A;jet=mul(jet,[ctx.mpf(1),-a,a*a-b])
        for rows in children:
            qa,qb=linear(rows[0])/2,linear(rows[1])/2
            jt=[T0,U1*(qa+qb),Cpair*(constant(rows[0])+constant(rows[1]))+Cfixed*constant(rows[2])+I2*(qa*qa+qa*qb+qb*qb)/3]
            jet=mul(jet,jt)
        return jet[2]
    point=coefficient(rad*ctx.cospi(theta),rad*ctx.sinpi(theta));r=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((r,0),(-r,0),(0,r),(0,-r)))+8*coefficient(0,0))/12
    return ctx,value+common*(mean/2-(1-z)*point)
