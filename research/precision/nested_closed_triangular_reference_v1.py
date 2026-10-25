"""Independent routed parent scale followed by the triangle endpoint transport."""
import mpmath as mp
import nested_closed_marginals_reference_v2 as B
M=B.M


def evaluate(point,binding,digits=640,mode='triangle-control',quadrature=False):
    assert mode in ('triangle-control','causal-control')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    s,z,theta,x0,y,v=p;outer_jac=ctx.mpf(1)
    if mode=='causal-control':
        x=x0*x0/(x0*x0+(1-x0)**2)
        a,b=binding['parent_boundary_shifts'];loops=[a[0],b[0],a[0]-b[0]]
        weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
        ell=sum(weight*l*l for weight,l in zip(weights,loops));h=1/ctx.sqrt(ell)
        denominator=1-s+h*s;outer_jac=h/(denominator*denominator);p[0]=h*s/denominator
    t=(p[0]/(1-p[0]))**2;L=ctx.log1p(t/2);eta=ctx.expm1(v*L)/t;jac=L*ctx.exp(v*L)/t
    total=ctx.mpf(0)
    for u in (eta,1-eta):
        pp=p.copy();pp[5]=u
        total+=B.B.direct(ctx,pp,binding,quadrature)+B.correction(ctx,pp,binding)
    return ctx,outer_jac*jac*total
