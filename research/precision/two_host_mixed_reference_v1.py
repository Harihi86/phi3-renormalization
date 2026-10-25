"""All seven original hosts, independent B derivatives, AD triangle and 24 frames."""
import mpmath as mp
import triangle_general_ad_reference_v1 as AD
import triangle_complete_general_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='mixed-pair',quadrature=False):
    assert mode in ('mixed-pair','mixed-scale')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    original=binding['original_binding'];inverse=binding['inverse_change']
    hosts=original['all_host_routes'];br=original['bubble_route'];rows=original['triangle_pair_differences'].copy()
    fixed=[i for i,r in enumerate(rows) if r[:2]==[0,0]]
    if fixed:assert len(fixed)==1;rows.append(rows.pop(fixed[0]))
    radial=lambda x:(x/(1-x))**2
    if mode=='mixed-scale':
        total=radial(p[0]);t,v=total*p[1],total*(1-p[1]);jac=total*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    radmeasure=t*t*v*v*jac/4;z,Z=p[2:4]
    ct,st=ctx.cospi(2*p[4]),ctx.sinpi(2*p[4]);cp,sp=ctx.cospi(2*p[5]),ctx.sinpi(2*p[5]);cpsi,spsi=ctx.cospi(p[6]),ctx.sinpi(p[6])
    e=ctx.matrix([ctx.sqrt(z)*ct,ctx.sqrt(z)*st,ctx.sqrt(1-z),0])
    n1=ctx.matrix([-st,ct,0,0]);n2=ctx.matrix([ctx.sqrt(1-z)*ct,ctx.sqrt(1-z)*st,-ctx.sqrt(z),0]);n3=ctx.matrix([0,0,0,1])
    n=ctx.sqrt(Z)*(cp*n1+sp*n2)+ctx.sqrt(1-Z)*n3
    u=ctx.sqrt(t)*e;vv=ctx.sqrt(v)*(cpsi*e+spsi*n)
    dot=lambda a,b:(a.T*b)[0]
    extnorm=lambda rr:ctx.mpf(rr[2]**2+rr[3]**2-rr[2]*rr[3])
    external=lambda rr:ctx.matrix([rr[2]-ctx.mpf('.5')*rr[3],ctx.sqrt(3)*rr[3]/2,0,0])
    project=lambda q,r:[q[0],(-q[0]+ctx.sqrt(3)*q[1])/2,r[0],(-r[0]+ctx.sqrt(3)*r[1])/2]
    linear=lambda rr,pr:2*(rr[0]*(rr[2]*pr[0]+rr[3]*pr[1])+rr[1]*(rr[2]*pr[2]+rr[3]*pr[3]))
    def triangle(args):
        if quadrature:return -ctx.quad(lambda x:(1-x)*M.triangle_log_marginal(ctx,x,*args),[0,ctx.mpf('.5'),1])
        return -H.triangle_log(ctx,*args)
    def bubble(w):
        if quadrature:return ctx.quad(lambda x:(1+w*x*(1-x))*ctx.log1p(w*x*(1-x))-w*x*(1-x),[0,ctx.mpf('.5'),1])
        return (w+4)**ctx.mpf('1.5')/(3*ctx.sqrt(w))*ctx.asinh(ctx.sqrt(w)/2)-4*w/9-ctx.mpf(4)/3
    def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
    total_residual=total_mean=ctx.mpf(0)
    for cq,cr,excluded in ((u,vv,2),(u,vv-u,1),(vv-u,u,0)):
        prop=[1+dot(cq,cq),1+dot(cr,cr),1+dot(cq+cr,cq+cr)]
        weight=prop[excluded]**3/sum(x**3 for x in prop)
        q=inverse[0][0]*cq+inverse[0][1]*cr;r=inverse[1][0]*cq+inverse[1][1]*cr
        route=lambda rr:rr[0]*q+rr[1]*r
        baseline=[1+dot(route(rr),route(rr)) for rr in hosts];denominator=ctx.fprod(baseline)
        bt=dot(route(br),route(br));B=bubble(bt)
        if quadrature:
            B1=ctx.quad(lambda x:x*(1-x)*ctx.log1p(bt*x*(1-x)),[0,ctx.mpf('.5'),1])
            B2=ctx.quad(lambda x:(x*(1-x))**2/(1+bt*x*(1-x)),[0,ctx.mpf('.5'),1])
        else:B1,B2=[ctx.diff(bubble,bt,i) for i in (1,2)]
        args=[dot(route(rr),route(rr)) for rr in rows];T0=-H.triangle_log(ctx,*args)
        if args[2]==0:
            assert args[0]==args[1]
            E1=ctx.diff(lambda x:H.edge_log(ctx,x),args[0]);E2=ctx.diff(lambda x:H.edge_log(ctx,x),args[0],2)
            ell,Cp,Cf,I2=-E1/2,-E1/4,-(1-(args[0]+3)*E1)/12,-E2
        else:unused,gradient,hessian=AD.invariant_jets(ctx,*args)
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
            for rr,d in zip(hosts,baseline):
                l,c=linear(rr,pr)/d,extnorm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
            bl=linear(br,pr);jet=mul(jet,[B,B1*bl,B1*extnorm(br)+B2*bl*bl/2])
            if args[2]==0:
                a,b=linear(rows[0],pr)/2,linear(rows[1],pr)/2
                factor=[T0,ell*(a+b),Cp*(extnorm(rows[0])+extnorm(rows[1]))+Cf*extnorm(rows[2])+I2*(a*a+a*b+b*b)/3]
            else:
                ll=[linear(rr,pr) for rr in rows];cc=[extnorm(rr) for rr in rows]
                factor=[T0,-sum(gradient[i]*ll[i] for i in range(3)),
                        -sum(gradient[i]*cc[i] for i in range(3))-sum(hessian[i][j]*ll[i]*ll[j] for i in range(3) for j in range(3))/2]
            return mul(jet,factor)
        point=coefficients(project(q,r))[2];bare=ctx.mpf(0)
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            factor=ctx.fprod(d/(1+dot(full(rr),full(rr))) for d,rr in zip(baseline,hosts))
            factor*=bubble(dot(full(br),full(br)))*triangle([dot(full(rr),full(rr)) for rr in rows])
            bare+=(factor-B*T0)/2
        tq,vr,hr=dot(q,q),dot(r,r),dot(q,r);root=ctx.sqrt(tq);transverse=ctx.sqrt(vr-hr*hr/tq)
        frames=[];linear_frames=[]
        for i in range(6):
            for sign in (-1,1):
                for side in (-1,1):
                    qq=[ctx.mpf(0)]*6;rr=[ctx.mpf(0)]*6;qq[i]=sign*root;rr[i]=sign*hr/root;rr[(i+1)%6]=side*transverse
                    c=coefficients(project(qq,rr));frames.append(c[2]);linear_frames.append(c[1])
        assert abs(sum(linear_frames))<ctx.power(10,-digits+30)*max(ctx.mpf(1),max(abs(x) for x in linear_frames))
        total_residual+=weight*(bare-point)/denominator;total_mean+=weight*sum(frames)/(24*denominator)
    return ctx,radmeasure*(ctx.mpf(8)/3*spsi**4)*(3*(1-z)*ctx.sqrt(1-Z)*total_residual+total_mean)
