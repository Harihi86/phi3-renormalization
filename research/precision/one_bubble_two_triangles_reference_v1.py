"""Four original host propagators,B and two independent complete triangle factors."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='one-bubble-control',quadrature=False):
    assert mode in ('one-bubble-direct','one-bubble-control') and binding['id'] in ('L4E3G249','L4E3G404')
    ctx=mp.mp.clone();ctx.dps=digits;s,z,theta=[ctx.mpf(v.numerator)/v.denominator for v in point]
    t=(s/(1-s))**2;A=1+t;host=binding['all_host_routes'];triangles=binding['triangle_pair_differences'];br=binding['bubble_route']
    constant=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    def bubble(w):
        if quadrature:return ctx.quad(lambda x:(1+w*x*(1-x))*ctx.log1p(w*x*(1-x))-w*x*(1-x),[0,ctx.mpf('.5'),1])
        return (w+4)**ctx.mpf('1.5')/(3*ctx.sqrt(w))*ctx.asinh(ctx.sqrt(w)/2)-4*w/9-ctx.mpf(4)/3
    def triangle(a,b):
        if quadrature:return -ctx.quad(lambda u:(1-u)*M.triangle_log_marginal(ctx,u,a,b,ctx.mpf(1)),[0,ctx.mpf('.5'),1])
        return -H.triangle_log_complete(ctx,a,b)
    B=bubble(t);T0=-H.edge_log(ctx,t)/2;P0=B*T0*T0;rad=ctx.sqrt(t*z);total=ctx.mpf(0)
    for sign in (-1,1):
        qx=sign*rad*ctx.cospi(theta);qy=sign*rad*ctx.sinpi(theta);q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        norm=lambda r:r[0]**2*t+2*r[0]*(r[1]*q1+r[2]*q2)+constant(r)
        ratio=ctx.fprod(A/(1+norm(r)) for r in host)
        total+=(ratio*bubble(norm(br))*ctx.fprod(triangle(norm(rows[0]),norm(rows[1])) for rows in triangles)-P0)/2
    common=2/(s*(1-s))*(t/A)**3/A;answer=common*(1-z)*total
    if mode=='one-bubble-direct':return ctx,answer
    if quadrature:
        B1=ctx.quad(lambda x:x*(1-x)*ctx.log1p(t*x*(1-x)),[0,ctx.mpf('.5'),1])
        B2=ctx.quad(lambda x:(x*(1-x))**2/(1+t*x*(1-x)),[0,ctx.mpf('.5'),1])
    else:B1,B2=[ctx.diff(bubble,t,i) for i in (1,2)]
    E1,E2=[ctx.diff(lambda w:H.edge_log(ctx,w),t,i) for i in (1,2)]
    U1=-E1/2;Cp=-E1/4;Cf=-(1-(t+3)*E1)/12;I2=-E2
    def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2;linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
        for r in host:
            l,c=linear(r)/A,constant(r)/A;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
        l=linear(br);jet=mul(jet,[B,B1*l,B1*constant(br)+B2*l*l/2])
        for rows in triangles:
            a,b=linear(rows[0])/2,linear(rows[1])/2
            jet=mul(jet,[T0,U1*(a+b),Cp*(constant(rows[0])+constant(rows[1]))+Cf*constant(rows[2])+I2*(a*a+a*b+b*b)/3])
        return jet[2]
    point=coefficient(rad*ctx.cospi(theta),rad*ctx.sinpi(theta));root=ctx.sqrt(t)
    mean=(sum(coefficient(a,b) for a,b in ((root,0),(-root,0),(0,root),(0,-root)))+8*coefficient(0,0))/12
    return ctx,answer+common*(mean/2-(1-z)*point)
