"""Original-route product and independent two-stage rotation-frame average."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M
import g371_fullangle_moments_v1 as A


def evaluate(point,binding,digits=640,mode='soft-pair',quadrature=False):
    assert binding['id']=='L4E3G371' and mode in ('soft-pair','soft-scale') and len(point)==7
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    radial=lambda x:(x/(1-x))**2
    if mode=='soft-scale':
        total=radial(p[0]);t,v=total*p[1],total*(1-p[1]);jac=total*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    radmeasure=t*t*v*v*jac/4;z,Z=p[2:4]
    ct,st=ctx.cospi(2*p[4]),ctx.sinpi(2*p[4]);cp,sp=ctx.cospi(2*p[5]),ctx.sinpi(2*p[5]);cpsi,spsi=ctx.cospi(p[6]),ctx.sinpi(p[6])
    e=ctx.matrix([ctx.sqrt(z)*ct,ctx.sqrt(z)*st,ctx.sqrt(1-z),0])
    n1=ctx.matrix([-st,ct,0,0]);n2=ctx.matrix([ctx.sqrt(1-z)*ct,ctx.sqrt(1-z)*st,-ctx.sqrt(z),0]);n3=ctx.matrix([0,0,0,1])
    n=ctx.sqrt(Z)*(cp*n1+sp*n2)+ctx.sqrt(1-Z)*n3
    u=ctx.sqrt(t)*e;vv=ctx.sqrt(v)*(cpsi*e+spsi*n)
    dot=lambda a,b:(a.T*b)[0]
    external=lambda rr:ctx.matrix([rr[2]-ctx.mpf('.5')*rr[3],ctx.sqrt(3)*rr[3]/2,0,0])
    extnorm=lambda rr:ctx.mpf(rr[2]**2+rr[3]**2-rr[2]*rr[3])
    project=lambda q,r:[q[0],(-q[0]+ctx.sqrt(3)*q[1])/2,r[0],(-r[0]+ctx.sqrt(3)*r[1])/2]
    other=[rr for rr in binding['host_routes'] if rr[0]*rr[1]==0];cross=[rr for rr in binding['host_routes'] if rr[0]*rr[1]!=0]
    assert len(other)==4 and cross==[[-1,-1,0,0],[1,1,1,1]]
    def triangle(a,b):
        if quadrature:return -ctx.quad(lambda x:(1-x)*M.triangle_log_marginal(ctx,x,a,b,ctx.mpf(1)),[0,ctx.mpf('.5'),1])
        return -H.triangle_log_complete(ctx,a,b)
    total_residual=total_mean=ctx.mpf(0)
    for q,r,excluded in ((u,vv,2),(u,vv-u,1),(vv-u,u,0)):
        tq,vr,hr=dot(q,q),dot(r,r),dot(q,r);sv=q+r;w=dot(sv,sv);X=1+w
        props=[1+tq,1+vr,X];weight=props[excluded]**3/sum(d**3 for d in props)
        route=lambda rr:rr[0]*q+rr[1]*r
        baseline=[1+dot(route(rr),route(rr)) for rr in other];d4=ctx.fprod(baseline)
        t0=[-H.edge_log(ctx,value)/2 for value in (tq,vr)];p0=ctx.fprod(t0);scalar=[]
        for value in (tq,vr):
            e1=ctx.diff(lambda x:H.edge_log(ctx,x),value);e2=ctx.diff(lambda x:H.edge_log(ctx,x),value,2)
            scalar.append((-e1/2,-e1/4,-(1-(value+3)*e1)/12,-e2))
        linear=lambda rr,pr:2*(rr[0]*(rr[2]*pr[0]+rr[3]*pr[1])+rr[1]*(rr[2]*pr[2]+rr[3]*pr[3]))
        def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
            for rr,d in zip(other,baseline):
                l,c=linear(rr,pr)/d,extnorm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
            for i,rows in enumerate(binding['child_pair_differences']):
                ell,Cp,Cf,I2=scalar[i];a,b=linear(rows[0],pr)/2,linear(rows[1],pr)/2
                jet=mul(jet,[t0[i],ell*(a+b),Cp*(extnorm(rows[0])+extnorm(rows[1]))+Cf*extnorm(rows[2])+I2*(a*a+a*b+b*b)/3])
            return jet
        pr=project(q,r);jet=coefficients(pr);residual=ctx.mpf(0)
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            factor=ctx.fprod(d/(1+dot(full(rr),full(rr))) for d,rr in zip(baseline,other))
            for rows in binding['child_pair_differences']:factor*=triangle(dot(full(rows[0]),full(rows[0])),dot(full(rows[1]),full(rows[1])))
            Y=1+dot(full(cross[1]),full(cross[1]))
            residual+=(factor-p0-sign*jet[1]-jet[2])/(2*X*Y)
        # Extract the low-degree polynomial independently from formal convolution.
        zero=[ctx.mpf(0)]*4
        approximate=lambda pr:sum(coefficients(pr))
        f0=approximate(zero);lin=[];diag=[];off={}
        for i in range(4):
            a=zero.copy();b=zero.copy();a[i]=1;b[i]=-1;fa,fb=approximate(a),approximate(b)
            lin.append((fa-fb)/2);diag.append((fa+fb)/2-f0)
        for i in range(4):
            for j in range(i+1,4):
                a=zero.copy();a[i]=a[j]=1;off[i,j]=(approximate(a)-f0-lin[i]-lin[j]-diag[i]-diag[j])/2
        def polynomial(pr):
            return f0+sum(l*x+d*x*x for l,d,x in zip(lin,diag,pr))+2*sum(c*pr[i]*pr[j] for (i,j),c in off.items())
        alpha,beta=(tq+hr)/w,(vr+hr)/w;rho=ctx.sqrt((tq*vr-hr*hr)/w);rootw=ctx.sqrt(w)
        def conditional(mu):
            # P is coordinate0 and (p1-p2)/sqrt(3) is coordinate1 here.
            # This is an independent global frame with the same external Gram.
            tangent=ctx.sqrt(1-mu*mu);values=[]
            for axis in range(1,6):
                for sign in (-1,1):
                    nn=[ctx.mpf(0)]*6;nn[0]=mu;nn[axis]=sign*tangent
                    b0=[ctx.mpf(0)]*6;b0[0]=-sign*tangent;b0[axis]=mu
                    bases=[b0]
                    for j in range(1,6):
                        if j!=axis:
                            bj=[ctx.mpf(0)]*6;bj[j]=1;bases.append(bj)
                    for tangent_axis in bases:
                        for side in (-1,1):
                            qq=[alpha*rootw*a+side*rho*b for a,b in zip(nn,tangent_axis)]
                            rr=[beta*rootw*a-side*rho*b for a,b in zip(nn,tangent_axis)]
                            pr=[(qq[0]+ctx.sqrt(3)*qq[1])/2,(qq[0]-ctx.sqrt(3)*qq[1])/2,(rr[0]+ctx.sqrt(3)*rr[1])/2,(rr[0]-ctx.sqrt(3)*rr[1])/2]
                            values.append(polynomial(pr))
            assert len(values)==100;return sum(values)/100
        c0=conditional(ctx.mpf(0));plus,minus=conditional(ctx.mpf('.5')),conditional(ctx.mpf('-.5'))
        c1=plus-minus;c2=2*(plus+minus-2*c0)
        for mu in (-1,1):
            actual=conditional(ctx.mpf(mu));assert abs(actual-(c0+c1*mu+c2*mu*mu))<=ctx.power(10,-digits+30)*max(1,abs(actual))
        aa=w+2;bb=2*ctx.sqrt(w)
        if quadrature:
            averaged=A.integral(ctx,lambda mu:(c0+c1*mu+c2*mu*mu)/(aa+bb*mu),aa,bb)
        else:
            m0=A.moments(ctx,aa,bb)[0];m1=(1-aa*m0)/bb;m2=-aa*m1/bb
            averaged=c0*m0+c1*m1+c2*m2
        total_residual+=weight*residual/d4;total_mean+=weight*(averaged/X-p0/(X*X))/d4
    orientation=3*(1-z)*ctx.sqrt(1-Z);intrinsic=ctx.mpf(8)/3*spsi**4
    return ctx,radmeasure*intrinsic*(orientation*total_residual+total_mean)
