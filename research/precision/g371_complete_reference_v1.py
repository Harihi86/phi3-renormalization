"""Independent original routes, direct child marginals and formal Taylor jets."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='seven-pair',quadrature=False):
    assert binding['id']=='L4E3G371' and mode in ('seven-pair','seven-scale','seven-direct')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    radial=lambda x:(x/(1-x))**2
    if mode=='seven-scale':
        scale=radial(p[0]);t,v=scale*p[1],scale*(1-p[1]);jac=scale*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    radmeasure=t*t*v*v*jac/4
    angular=8*(1-p[2])*(1-p[3])*ctx.sinpi(p[6])**2
    u=ctx.matrix([ctx.sqrt(t*p[2])*ctx.cospi(2*p[4]),ctx.sqrt(t*p[2])*ctx.sinpi(2*p[4]),ctx.sqrt(t*(1-p[2])),0])
    vvec=ctx.matrix([ctx.sqrt(v*p[3])*ctx.cospi(2*p[5]),ctx.sqrt(v*p[3])*ctx.sinpi(2*p[5]),ctx.sqrt(v*(1-p[3]))*ctx.cospi(p[6]),ctx.sqrt(v*(1-p[3]))*ctx.sinpi(p[6])])
    dot=lambda a,b:(a.T*b)[0]
    external=lambda rr:ctx.matrix([rr[2]-ctx.mpf('.5')*rr[3],ctx.sqrt(3)*rr[3]/2,0,0])
    total=ctx.mpf(0)
    for q,r,excluded in ((u,vvec,2),(u,vvec-u,1),(vvec-u,u,0)):
        tq,tr,qr=dot(q,q),dot(r,r),dot(q,r)
        prop=[1+tq,1+tr,1+dot(q+r,q+r)]
        weight=prop[excluded]**3/sum(x**3 for x in prop)
        route=lambda rr:rr[0]*q+rr[1]*r
        baseline=[1+dot(route(rr),route(rr)) for rr in binding['host_routes']]
        t0=[-H.edge_log(ctx,x)/2 for x in (tq,tr)]
        value=ctx.mpf(0)
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            ratio=ctx.fprod(d/(1+dot(full(rr),full(rr))) for d,rr in zip(baseline,binding['host_routes']))
            factors=[]
            for rows in binding['child_pair_differences']:
                a,b=[dot(full(rr),full(rr)) for rr in rows[:2]]
                if quadrature:factors.append(-ctx.quad(lambda x:(1-x)*M.triangle_log_marginal(ctx,x,a,b,ctx.mpf(1)),[0,ctx.mpf('.5'),1]))
                else:factors.append(-H.triangle_log_complete(ctx,a,b))
            value+=(ratio*ctx.fprod(factors)-ctx.fprod(t0))/2
        if mode!='seven-direct':
            moments=[]
            for x in (tq,tr):
                if quadrature:
                    ell=lambda u:u*(1-u);br=[0,ctx.mpf('.5'),1]
                    u1=ctx.quad(lambda u:-(1-u)*ell(u)/(1+x*ell(u)),br)
                    cp=u1/2;cf=ctx.quad(lambda u:-(1-u)**3/(6*(1+x*ell(u))),br)
                    i2=ctx.quad(lambda u:ell(u)**2/(1+x*ell(u))**2,br)
                else:
                    e1=ctx.diff(lambda z:H.edge_log(ctx,z),x);e2=ctx.diff(lambda z:H.edge_log(ctx,z),x,2)
                    u1=-e1/2;cp=-e1/4;cf=-(1-(x+3)*e1)/12;i2=-e2
                moments.append((u1,cp,cf,i2))
            def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
            def coefficient(projections):
                jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
                def linear(rr):return 2*(rr[0]*(rr[2]*projections[0]+rr[3]*projections[1])+rr[1]*(rr[2]*projections[2]+rr[3]*projections[3]))
                norm=lambda rr:rr[2]**2+rr[3]**2-rr[2]*rr[3]
                for rr,d in zip(binding['host_routes'],baseline):
                    l,c=linear(rr)/d,norm(rr)/d;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
                for i,rows in enumerate(binding['child_pair_differences']):
                    u1,cp,cf,i2=moments[i];a,b=linear(rows[0])/2,linear(rows[1])/2
                    child=[t0[i],u1*(a+b),cp*(norm(rows[0])+norm(rows[1]))+cf*norm(rows[2])+i2*(a*a+a*b+b*b)/3]
                    jet=mul(jet,child)
                return jet[2]
            project=lambda q,r:[q[0],(-q[0]+ctx.sqrt(3)*q[1])/2,r[0],(-r[0]+ctx.sqrt(3)*r[1])/2]
            point2=coefficient(project(q,r));frames=[]
            # Twelve frames reproduce the exact quadratic rotation average.
            for i in range(6):
                for sign in (-1,1):
                    qq=[ctx.mpf(0)]*6;rr=[ctx.mpf(0)]*6
                    qq[i]=ctx.sqrt(tq);rr[i]=qr/ctx.sqrt(tq);rr[(i+1)%6]=sign*ctx.sqrt(tr-qr*qr/tq)
                    frames.append(coefficient(project(qq,rr)))
            value+=sum(frames)/12-point2
        total+=weight*value/ctx.fprod(baseline)
    return ctx,total*radmeasure*angular
