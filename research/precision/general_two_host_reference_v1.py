"""Independent original-route product with explicit G370 loop reflection."""
import mpmath as mp
import triangle_complete_general_v1 as H
import nested_analytic_marginals_v1 as M


def evaluate(point,binding,digits=640,mode='general-pair',quadrature=False):
    assert binding['id'] in ('L4E3G244','L4E3G370') and mode in ('general-pair','general-scale')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    radial=lambda x:(x/(1-x))**2
    if mode=='general-scale':
        scale=radial(p[0]);t,v=scale*p[1],scale*(1-p[1]);jac=scale*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    radmeasure=t*t*v*v*jac/4;angular=8*(1-p[2])*(1-p[3])*ctx.sinpi(p[6])**2
    u=ctx.matrix([ctx.sqrt(t*p[2])*ctx.cospi(2*p[4]),ctx.sqrt(t*p[2])*ctx.sinpi(2*p[4]),ctx.sqrt(t*(1-p[2])),0])
    vv=ctx.matrix([ctx.sqrt(v*p[3])*ctx.cospi(2*p[5]),ctx.sqrt(v*p[3])*ctx.sinpi(2*p[5]),ctx.sqrt(v*(1-p[3]))*ctx.cospi(p[6]),ctx.sqrt(v*(1-p[3]))*ctx.sinpi(p[6])])
    dot=lambda a,b:(a.T*b)[0]
    external=lambda rr:ctx.matrix([rr[2]-ctx.mpf('.5')*rr[3],ctx.sqrt(3)*rr[3]/2,0,0])
    def triangle(args):
        a,b,c=args
        if quadrature:return -ctx.quad(lambda x:(1-x)*M.triangle_log_marginal(ctx,x,a,b,c),[0,ctx.mpf('.5'),1])
        return -H.triangle_log(ctx,a,b,c)
    total=ctx.mpf(0)
    for cq,cr,excluded in ((u,vv,2),(u,vv-u,1),(vv-u,u,0)):
        prop=[1+dot(cq,cq),1+dot(cr,cr),1+dot(cq+cr,cq+cr)]
        weight=prop[excluded]**3/sum(x**3 for x in prop)
        q,r=cq,(-cr if binding['id']=='L4E3G370' else cr)
        route=lambda rr:rr[0]*q+rr[1]*r
        baseline=[1+dot(route(rr),route(rr)) for rr in binding['host_routes']]
        t0=[triangle([dot(route(rr),route(rr)) for rr in rows]) for rows in binding['child_pair_differences']]
        value=ctx.mpf(0)
        for sign in (-1,1):
            full=lambda rr:sign*route(rr)+external(rr)
            ratio=ctx.fprod(d/(1+dot(full(rr),full(rr))) for d,rr in zip(baseline,binding['host_routes']))
            factors=[triangle([dot(full(rr),full(rr)) for rr in rows]) for rows in binding['child_pair_differences']]
            value+=(ratio*ctx.fprod(factors)-ctx.fprod(t0))/2
        total+=weight*value/ctx.fprod(baseline)
    return ctx,total*radmeasure*angular
