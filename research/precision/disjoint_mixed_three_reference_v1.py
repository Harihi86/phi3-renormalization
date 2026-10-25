"""Original routes and two independent T2 bubbles, with formal Taylor products."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='mixed-control',quadrature=False):
    assert mode in ('mixed-direct','mixed-control') and binding['id'] in ('L4E3G204','L4E3G253','L4E3G411','L4E3G465')
    ctx=mp.mp.clone();ctx.dps=digits;s,z,theta=[ctx.mpf(v.numerator)/v.denominator for v in point]
    t=(s/(1-s))**2;A=1+t;host=binding['all_host_routes'];rows=binding['triangle_pair_differences'];bubbles=binding['bubble_routes']
    constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    def bubble(w):
        if quadrature:return ctx.quad(lambda x:(1+w*x*(1-x))*ctx.log1p(w*x*(1-x))-w*x*(1-x),[0,ctx.mpf('.5'),1])
        return (w+4)**ctx.mpf('1.5')/(3*ctx.sqrt(w))*ctx.asinh(ctx.sqrt(w)/2)-4*w/9-ctx.mpf(4)/3
    # Keep all five original host denominators here, unlike the native J factors.
    B=bubble(t);T0=-H.edge_log(ctx,t)/2;P0=B*B*T0;rad=ctx.sqrt(t*z);total=ctx.mpf(0)
    for sign in (-1,1):
        qx=sign*rad*ctx.cospi(theta);qy=sign*rad*ctx.sinpi(theta);q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        norm=lambda r:r[0]**2*t+2*r[0]*(r[1]*q1+r[2]*q2)+constant(r)
        a,b=norm(rows[0]),norm(rows[1])
        if quadrature:tri=-ctx.quad(lambda u:(1-u)*M.triangle_log_marginal(ctx,u,a,b,ctx.mpf(1)),[0,ctx.mpf('.5'),1])
        else:tri=-H.triangle_log_complete(ctx,a,b)
        ratio=ctx.fprod(A/(1+norm(r)) for r in host)
        total+=(ratio*ctx.fprod(bubble(norm(r)) for r in bubbles)*tri-P0)/2
    common=2/(s*(1-s))*(t/A)**3/A**2;answer=common*(1-z)*total
    if mode=='mixed-direct':return ctx,answer
    if quadrature:
        B1=ctx.quad(lambda x:x*(1-x)*ctx.log1p(t*x*(1-x)),[0,ctx.mpf('.5'),1])
        B2=ctx.quad(lambda x:(x*(1-x))**2/(1+t*x*(1-x)),[0,ctx.mpf('.5'),1])
    else:B1,B2=[ctx.diff(bubble,t,i) for i in (1,2)]
    E1,E2=[ctx.diff(lambda w:H.edge_log(ctx,w),t,i) for i in (1,2)]
    U1=-E1/2;Cpair=-E1/4;Cfixed=-(1-(t+3)*E1)/12;I2=-E2
    def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2;linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            l,c=linear(r)/A,constant(r)/A;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
        for r in bubbles:
            l=linear(r);jet=mul(jet,[B,B1*l,B1*constant(r)+B2*l*l/2])
        a,b=linear(rows[0])/2,linear(rows[1])/2
        jet=mul(jet,[T0,U1*(a+b),Cpair*(constant(rows[0])+constant(rows[1]))+Cfixed*constant(rows[2])+I2*(a*a+a*b+b*b)/3])
        return jet[2]
    point=coefficient(rad*ctx.cospi(theta),rad*ctx.sinpi(theta));root=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((root,0),(-root,0),(0,root),(0,-root)))+8*coefficient(0,0))/12
    return ctx,answer+common*(mean/2-(1-z)*point)
