"""Conditional radial chart from original routed child variances."""
import mpmath as mp
import nested_closed_marginals_reference_v2 as B
M=B.M


def evaluate(point,binding,digits=640,mode='radial-control',quadrature=False):
    assert mode in ('radial-control','radial-mixture-control')
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    s,z,theta,x0,y,u=p;x=x0*x0/(x0*x0+(1-x0)**2)
    a,b=binding['parent_boundary_shifts'];loops=[a[0],b[0],a[0]-b[0]]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    ell=sum(w*l*l for w,l in zip(weights,loops));tri=u*(1-u);h=(ell*tri)**ctx.mpf('-.25')
    den=1-s+h*s;w=h*s/den;jac=h/(den*den)
    def density(q):return B.B.direct(ctx,q,binding,quadrature)+B.correction(ctx,q,binding)
    pp=p.copy();pp[0]=w;scaled=jac*density(pp)
    if mode=='radial-control':return ctx,scaled
    ratio0=h/(h*(1-s)+s)**2;ratio1=den*den/h
    return ctx,density(p)/(1+ratio0)+scaled*ratio1/(1+ratio1)
