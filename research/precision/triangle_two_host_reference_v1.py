"""Original routes, scalar K autodifferentiation, direct primitive integrals and 24 frames."""
import mpmath as mp
import triangle_child_radial_reference_v2 as F
import triangle_general_ad_reference_v1 as AD
import triangle_complete_general_v1 as H
import nested_analytic_marginals_v1 as M
import nested_bubble_vertex_factor_v1 as N


def evaluate(point,binding,digits=640,mode='triangle-child-pair',quadrature=False):
    assert mode in ('triangle-child-pair','triangle-child-scale')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    original=binding['original_binding'];inner_sign=binding['inner_loop_sign']
    parent=original['parent_index'] is not None
    outer=original['outer_host_routes'] if parent else original['host_routes'];inner=original['middle_host_routes']
    inners=original['child_inner_shifts'];triangles=original['child_triangle_routes']
    zero_inners=original['parent_T0_child_inner_shifts'];zero_triangles=original['parent_T0_child_triangle_routes']
    fmap=lambda u:u*u/(u*u+(1-u)**2)
    x,y=fmap(p[7]),p[8];u,w=fmap(p[9]),p[10]
    a=[x,(1-x)*y,(1-x)*(1-y)];b=[u,(1-u)*w,(1-u)*(1-w)]
    weights=[b[0]*b[1],b[0]*b[2],b[1]*b[2]]
    auxiliary=(1-x)*ctx.diff(fmap,p[7])*(1-u)*ctx.diff(fmap,p[9])
    radial=lambda x:(x/(1-x))**2
    if mode=='triangle-child-scale':
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
    def child(args):return F.direct(ctx,*args) if quadrature else F.closed(ctx,*args)
    def geometry(inner,tri,pr):
        center=[sum(ai*ctx.mpf(row[j]) for ai,row in zip(a,inner)) for j in range(4)]
        ell=sum(w*row[0]**2 for w,row in zip(weights,tri))
        c=[sum(w*row[0]*row[j+1] for w,row in zip(weights,tri)) for j in range(4)]
        def normcoef(rr):return [dot(route(rr),route(rr)),linear(rr,pr),extnorm(rr)]
        output=[[ctx.mpf(0)]*3 for _ in range(3)]
        # Q uses the original host Schur complement, not pair differences.
        for ai,row in zip(a,inner):
            nn=normcoef(row)
            for j in range(3):output[0][j]+=ai*nn[j]
        nn=normcoef(center)
        for j in range(3):output[0][j]-=nn[j]
        for wi,row in zip(weights,tri):
            nn=normcoef([ctx.mpf(row[j+1])-row[0]*center[j] for j in range(4)])
            ee=normcoef(row[1:])
            for j in range(3):output[1][j]+=wi*nn[j];output[2][j]+=wi*ee[j]
        nn=normcoef(c)
        for j in range(3):output[2][j]-=nn[j]/ell
        return ell,output
    def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
    total_residual=total_mean=ctx.mpf(0)
    for cq,cr,excluded in ((u,vv,2),(u,vv-u,1),(vv-u,u,0)):
        prop=[1+dot(cq,cq),1+dot(cr,cr),1+dot(cq+cr,cq+cr)];weight=prop[excluded]**3/sum(x**3 for x in prop)
        q,r=cq,inner_sign*cr;route=lambda rr:rr[0]*q+rr[1]*r
        hosts=outer+inner;baselines=[1+dot(route(rr),route(rr)) for rr in hosts]
        outer_den=ctx.fprod(baselines[:len(outer)]);inner_den=ctx.fprod(baselines[len(outer):]);den=outer_den*inner_den
        vr=dot(r,r);ell,geo0=geometry(inners,triangles,project(q,r))
        args=[geo0[0][0],ell,geo0[1][0],geo0[2][0]]
        scalar,gradient,hessian=F.invariant_jets(ctx,*args);K0=child(args) if quadrature else scalar;parent0=ctx.mpf(0)
        if parent:
            ellz,geo=geometry(zero_inners,zero_triangles,[ctx.mpf(0)]*4)
            parent0=child([geo[0][0],ellz,geo[1][0],geo[2][0]])/(1+vr)**3
        def child_coeff(pr):
            ell,geo=geometry(inners,triangles,pr);ll=[g[1] for g in geo];cc=[g[2] for g in geo]
            return [K0,sum(gradient[i]*ll[i] for i in range(3)),
                    sum(gradient[i]*cc[i] for i in range(3))+sum(hessian[i][j]*ll[i]*ll[j] for i in range(3) for j in range(3))/2]
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)];outer_jet=None
            for index,(rr,d) in enumerate(zip(hosts,baselines)):
                l,c=linear(rr,pr)/d,extnorm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
                if index+1==len(outer):outer_jet=jet.copy()
            a=mul(jet,child_coeff(pr))
            b=outer_jet
            return [a[i]/den-parent0*b[i]/outer_den for i in range(3)]
        point_coef=coefficients(project(q,r))[2];bare=ctx.mpf(0)
        baseline_value=(K0/inner_den-parent0)/outer_den
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            shifted_den=[1+dot(full(rr),full(rr)) for rr in hosts]
            outer_product=ctx.fprod(shifted_den[:len(outer)]);inner_product=ctx.fprod(shifted_den[len(outer):])
            ell,geo=geometry(inners,triangles,project(q,r));scalars=[v[0]+sign*v[1]+v[2] for v in geo]
            value=(child([scalars[0],ell,scalars[1],scalars[2]])/inner_product-parent0)/outer_product
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
    return ctx,auxiliary*radmeasure*(ctx.mpf(8)/3*spsi**4)*(3*(1-z)*ctx.sqrt(1-Z)*total_residual+total_mean)
