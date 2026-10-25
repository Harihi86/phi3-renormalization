"""Floating G35 diagnostic; no bounds or precision promotion."""
import numpy as np
import two_loop_vertex_rank_one_v1 as F
import triangle_vertex_vector_diagnostic_v3 as COMMON
polynomial=COMMON.polynomial
original=COMMON.original
condition=COMMON.condition


def mapped(points,permutation,transport):
    assert transport in ('uniform','endpoint') and tuple(permutation) in F.PERMUTATIONS
    points=np.asarray(points,dtype=np.longdouble)
    assert points.ndim==2 and points.shape[1]==8 and np.isfinite(points).all()
    assert ((points>0)&(points<1)).all()
    s,v,f0,f1,f2,r0,t,z0=points.T;den=1+s+s*v;weights=[None]*3
    for i,j in enumerate(permutation):weights[j]=(1/den,s/den,s*v/den)[i]
    w0,w1,w2=weights;psi=w0*w1+w0*w2+w1*w2
    fs=(f0,f1,f2);hh=[w*f*(1-f) for w,f in zip(weights,fs)];kappa=w0*w1*w2/psi
    x,y=f0-f1,f0-f2;a=hh[0]+hh[1]+kappa*x*x;cross=hh[0]+kappa*x*y
    determinant=hh[0]*hh[1]+hh[0]*hh[2]+hh[1]*hh[2]+kappa*(hh[0]*(y-x)**2+hh[1]*y*y+hh[2]*x*x)
    if transport=='endpoint':
        rd=r0*r0+(1-r0)**2;r=r0*r0/rd;rc=(1-r0)**2/rd;dr=2*r0*(1-r0)/rd**2
        z=z0*z0;y0=(1-z0)*(1+z0);dy=2*z0
    else:r,rc,dr,z,y0,dy=r0,1-r0,np.ones_like(r0),z0,1-z0,np.ones_like(z0)
    host=(r,rc*t,rc*(1-t));q=host[0]*host[1]+host[0]*host[2]+host[1]*host[2]
    n=a*((cross/a-host[0]-host[2]-(cross/a-host[0])/2)**2+3*(cross/a-host[0])**2/4)
    d=determinant/a;u=a;cc=1+u;sigma=y0/(cc*z)
    d1,d2,d3=u+z,1+u*z,1+z*(u+q*cc)
    x0=q*cc*z/d2;x1=y0*(n*cc*z/d1+d)/d3
    ca=u*(d1*d1+d1*cc*z+(cc*z)**2)/d1**3*q*cc/d2
    cb=cc**3*z*z/d1**3*(n*cc*z/d1+d)/d3
    def l(x):
        out=np.ones_like(x);np.divide(np.log1p(x),x,out=out,where=x!=0);return out
    common=rc*dr*dy;measure=v/(1+v+s*v)**3*common
    positive,negative=measure*ca*l(x0),measure*cb*l(x1)
    alphas=[sigma*weights[i]*fs[i] for i in range(3)]+[sigma*weights[i]*(1-fs[i]) for i in range(3)]+list(host)
    jac=w0*w1*w2*s/den**3*common*sigma**5/(cc*z*z)
    assert np.isfinite(positive).all() and np.isfinite(negative).all()
    return {'value':positive-negative,'positive':positive,'negative':negative,
        'alphas':np.asarray(alphas),'jacobian':jac,'mass_sum':1+sigma}
