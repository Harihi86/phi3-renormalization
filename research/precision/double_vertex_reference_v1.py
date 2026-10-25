"""All original host denominators, two independent children and the parent T0."""
import mpmath as mp
import triangle_general_ad_reference_v1 as AD
import triangle_complete_general_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='double-pair',quadrature=False):
    assert mode in ('double-pair','double-scale')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    original=binding['original_binding'];inner_sign=binding['inner_loop_sign']
    outer=original['outer_host_routes'];inner=original['all_inner_host_routes'];kinds=original['child_kinds']
    def ordered(rows):
        rows=[list(r) for r in rows];fixed=[i for i,r in enumerate(rows) if r[:2]==[0,0]]
        if fixed:
            assert len(fixed)==1;rows.append(rows.pop(fixed[0]))
            for i in range(2):
                if next(v for v in rows[i][:2] if v)<0:rows[i]=[-v for v in rows[i]]
            assert rows[0][:2]==rows[1][:2]
        return rows
    child_rows=[rows if kind=='bubble' else ordered(rows) for kind,rows in zip(kinds,original['child_routes'])]
    radial=lambda x:(x/(1-x))**2
    if mode=='double-scale':
        total=radial(p[0]);t,v=total*p[1],total*(1-p[1]);jac=total*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    radmeasure=t*t*v*v*jac/4;z,Z=p[2:4]
    ct,st=ctx.cospi(2*p[4]),ctx.sinpi(2*p[4]);cp,sp=ctx.cospi(2*p[5]),ctx.sinpi(2*p[5]);cpsi,spsi=ctx.cospi(p[6]),ctx.sinpi(p[6])
    e=ctx.matrix([ctx.sqrt(z)*ct,ctx.sqrt(z)*st,ctx.sqrt(1-z),0])
    n1=ctx.matrix([-st,ct,0,0]);n2=ctx.matrix([ctx.sqrt(1-z)*ct,ctx.sqrt(1-z)*st,-ctx.sqrt(z),0]);n3=ctx.matrix([0,0,0,1])
    n=ctx.sqrt(Z)*(cp*n1+sp*n2)+ctx.sqrt(1-Z)*n3;u=ctx.sqrt(t)*e;vv=ctx.sqrt(v)*(cpsi*e+spsi*n)
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
        prop=[1+dot(cq,cq),1+dot(cr,cr),1+dot(cq+cr,cq+cr)];weight=prop[excluded]**3/sum(x**3 for x in prop)
        q,r=cq,inner_sign*cr;route=lambda rr:rr[0]*q+rr[1]*r
        hosts=outer+inner;baselines=[1+dot(route(rr),route(rr)) for rr in hosts]
        outer_den=ctx.fprod(baselines[:len(outer)]);inner_den=ctx.fprod(baselines[len(outer):]);den=outer_den*inner_den
        vr=dot(r,r);parent0=ctx.fprod(bubble(vr)/(1+vr) if kind=='bubble' else (triangle([vr,vr,ctx.mpf(0)]) if quadrature else -H.edge_log(ctx,vr)/2) for kind in kinds)/(1+vr)**3
        def triangle_coefficients(rows):
            args=[dot(route(rr),route(rr)) for rr in rows];zero=-H.triangle_log(ctx,*args)
            if args[2]==0:
                assert args[0]==args[1]
                E1=ctx.diff(lambda x:H.edge_log(ctx,x),args[0]);E2=ctx.diff(lambda x:H.edge_log(ctx,x),args[0],2)
                ell,Cp,Cf,I2=-E1/2,-E1/4,-(1-(args[0]+3)*E1)/12,-E2
                def coeff(pr):
                    a,b=linear(rows[0],pr)/2,linear(rows[1],pr)/2
                    return [zero,ell*(a+b),Cp*(extnorm(rows[0])+extnorm(rows[1]))+Cf*extnorm(rows[2])+I2*(a*a+a*b+b*b)/3]
            else:
                unused,gradient,hessian=AD.invariant_jets(ctx,*args)
                def coeff(pr):
                    ll=[linear(rr,pr) for rr in rows];cc=[extnorm(rr) for rr in rows]
                    return [zero,-sum(gradient[i]*ll[i] for i in range(3)),-sum(gradient[i]*cc[i] for i in range(3))-sum(hessian[i][j]*ll[i]*ll[j] for i in range(3) for j in range(3))/2]
            return zero,coeff
        def child_coefficients(kind,rows):
            if kind=='triangle':return triangle_coefficients(rows)
            br=rows[0];bt=dot(route(br),route(br));zero=bubble(bt)
            if quadrature:
                b1=ctx.quad(lambda x:x*(1-x)*ctx.log1p(bt*x*(1-x)),[0,ctx.mpf('.5'),1]);b2=ctx.quad(lambda x:(x*(1-x))**2/(1+bt*x*(1-x)),[0,ctx.mpf('.5'),1])
            else:b1,b2=[ctx.diff(bubble,bt,i) for i in (1,2)]
            def coeff(pr):
                l=linear(br,pr);return [zero,b1*l,b1*extnorm(br)+b2*l*l/2]
            return zero,coeff
        primitive=[child_coefficients(kind,rows) for kind,rows in zip(kinds,child_rows)]
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)];outer_jet=None
            for index,(rr,d) in enumerate(zip(hosts,baselines)):
                l,c=linear(rr,pr)/d,extnorm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
                if index+1==len(outer):outer_jet=jet.copy()
            a=jet
            for zero,coeff in primitive:a=mul(a,coeff(pr))
            b=outer_jet
            return [a[i]/den-parent0*b[i]/outer_den for i in range(3)]
        point_coef=coefficients(project(q,r))[2];bare=ctx.mpf(0)
        baseline_value=(ctx.fprod(item[0] for item in primitive)/inner_den-parent0)/outer_den
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            shifted_den=[1+dot(full(rr),full(rr)) for rr in hosts]
            outer_product=ctx.fprod(shifted_den[:len(outer)]);inner_product=ctx.fprod(shifted_den[len(outer):])
            children=ctx.fprod(bubble(dot(full(rows[0]),full(rows[0]))) if kind=='bubble' else triangle([dot(full(rr),full(rr)) for rr in rows]) for kind,rows in zip(kinds,child_rows))
            value=(children/inner_product-parent0)/outer_product
            bare+=(value-baseline_value)/2
        tq,hr=dot(q,q),dot(q,r);root=ctx.sqrt(tq);transverse=ctx.sqrt(vr-hr*hr/tq);frames=[];linear_frames=[]
        for axis in range(6):
            for sign in (-1,1):
                for side in (-1,1):
                    qq=[ctx.mpf(0)]*6;rr=qq.copy();qq[axis]=sign*root;rr[axis]=sign*hr/root;rr[(axis+1)%6]=side*transverse
                    coeff=coefficients(project(qq,rr));frames.append(coeff[2]);linear_frames.append(coeff[1])
        assert abs(sum(linear_frames))<ctx.power(10,-digits+30)*max(ctx.mpf(1),max(abs(x) for x in linear_frames))
        chi=(1+((tq+vr)/64)**4)**(-3)
        total_residual+=weight*(bare-chi*point_coef);total_mean+=weight*chi*sum(frames)/24
    return ctx,radmeasure*(ctx.mpf(8)/3*spsi**4)*(3*(1-z)*ctx.sqrt(1-Z)*total_residual+total_mean)
