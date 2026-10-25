"""Original routes, differentiated scalar, formal convolution and 24 rotation frames."""
import mpmath as mp
import triangle_complete_general_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='angular-pair',quadrature=False):
    assert binding['id'] in ('L4E3G244','L4E3G370') and mode in ('angular-pair','angular-scale')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    radial=lambda x:(x/(1-x))**2
    if mode=='angular-scale':
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
    def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
    total_residual=total_mean=ctx.mpf(0)
    for cq,cr,excluded in ((u,vv,2),(u,vv-u,1),(vv-u,u,0)):
        prop=[1+dot(cq,cq),1+dot(cr,cr),1+dot(cq+cr,cq+cr)]
        weight=prop[excluded]**3/sum(x**3 for x in prop)
        q,r=cq,(-cr if binding['id']=='L4E3G370' else cr)
        route=lambda rr:rr[0]*q+rr[1]*r
        baseline=[1+dot(route(rr),route(rr)) for rr in binding['host_routes']]
        denominator=ctx.fprod(baseline);scalars=[];t0=[]
        for rows in binding['child_pair_differences']:
            args=[dot(route(rr),route(rr)) for rr in rows];value=-H.triangle_log(ctx,*args);t0.append(value)
            if args[2]==0:
                assert args[0]==args[1]
                E1=ctx.diff(lambda x:H.edge_log(ctx,x),args[0]);E2=ctx.diff(lambda x:H.edge_log(ctx,x),args[0],2)
                scalars.append(('unit',(-E1/2,-E1/4,-(1-(args[0]+3)*E1)/12,-E2)))
            else:
                gradient=[ctx.diff(lambda a,b,c:H.triangle_log(ctx,a,b,c),tuple(args),tuple(int(k==i) for k in range(3))) for i in range(3)]
                hessian=[[None]*3 for _ in range(3)]
                for i in range(3):
                    for j in range(i,3):
                        hessian[i][j]=hessian[j][i]=ctx.diff(lambda a,b,c:H.triangle_log(ctx,a,b,c),tuple(args),tuple(int(k==i)+int(k==j) for k in range(3)))
                scalars.append(('general',(gradient,hessian)))
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
            for rr,d in zip(binding['host_routes'],baseline):
                l,c=linear(rr,pr)/d,extnorm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
            for k,rows in enumerate(binding['child_pair_differences']):
                kind,scalar=scalars[k]
                if kind=='unit':
                    ell,Cp,Cf,I2=scalar;a,b=linear(rows[0],pr)/2,linear(rows[1],pr)/2
                    factor=[t0[k],ell*(a+b),Cp*(extnorm(rows[0])+extnorm(rows[1]))+Cf*extnorm(rows[2])+I2*(a*a+a*b+b*b)/3]
                else:
                    gradient,hessian=scalar;ll=[linear(rr,pr) for rr in rows];cc=[extnorm(rr) for rr in rows]
                    factor=[t0[k],-sum(gradient[i]*ll[i] for i in range(3)),
                            -sum(gradient[i]*cc[i] for i in range(3))-sum(hessian[i][j]*ll[i]*ll[j] for i in range(3) for j in range(3))/2]
                jet=mul(jet,factor)
            return jet
        pr=project(q,r);point=coefficients(pr)[2];p0=ctx.fprod(t0);bare=ctx.mpf(0)
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            factor=ctx.fprod(d/(1+dot(full(rr),full(rr))) for d,rr in zip(baseline,binding['host_routes']))
            for rows in binding['child_pair_differences']:factor*=triangle([dot(full(rr),full(rr)) for rr in rows])
            bare+=(factor-p0)/2
        tq,vr,hr=dot(q,q),dot(r,r),dot(q,r);root=ctx.sqrt(tq);transverse=ctx.sqrt(vr-hr*hr/tq)
        frames=[];linear_frames=[]
        for i in range(6):
            for sign in (-1,1):
                for side in (-1,1):
                    qq=[ctx.mpf(0)]*6;rr=[ctx.mpf(0)]*6
                    qq[i]=sign*root;rr[i]=sign*hr/root;rr[(i+1)%6]=side*transverse
                    coeff=coefficients(project(qq,rr));frames.append(coeff[2]);linear_frames.append(coeff[1])
        assert len(frames)==24 and abs(sum(linear_frames))<ctx.power(10,-digits+30)*max(ctx.mpf(1),max(abs(x) for x in linear_frames))
        mean=sum(frames)/24
        total_residual+=weight*(bare-point)/denominator;total_mean+=weight*mean/denominator
    orientation=3*(1-z)*ctx.sqrt(1-Z);intrinsic=ctx.mpf(8)/3*spsi**4
    return ctx,radmeasure*intrinsic*(orientation*total_residual+total_mean)
