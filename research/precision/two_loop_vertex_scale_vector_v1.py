"""Fast seven-dimensional diagnostic only; native Arb remains the evaluator."""
import numpy as np
from scipy.special import spence
import two_loop_vertex_rank_one_v1 as F


def K(value):
    v=np.asarray(value,dtype=np.float64);assert (v>0).all()
    delta=v-1;near=np.abs(delta)<=.5;result=np.empty_like(v)
    if near.any():
        d=delta[near];p=d.copy();total=np.full_like(d,np.pi**2/6-1.25)
        for n in range(1,65):
            total+=2*(-1)**(n+1)*p/(n*n*(n+1)*(n+2));p*=d
        result[near]=total
    if (~near).any():
        z=v[~near]
        result[~near]=np.pi**2/6-spence(z)-z/(2*(1-z))+z*(2-3*z)*np.log(z)/(2*(1-z)**2)
    return result


def parameters(points,perm,transport='endpoint'):
    assert transport in ('uniform','endpoint')
    s,v,f0,f1,f2,r0,t=np.asarray(points,dtype=np.longdouble).T
    den=1+s+s*v;w=[None]*3
    for i,j in enumerate(perm):w[j]=(1/den,s/den,s*v/den)[i]
    psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2];kap=w[0]*w[1]*w[2]/psi
    fs=(f0,f1,f2);h=[a*f*(1-f) for a,f in zip(w,fs)]
    x,y=f0-f1,f0-f2;u=h[0]+h[1]+kap*x*x;c=h[0]+kap*x*y
    det=h[0]*h[1]+h[0]*h[2]+h[1]*h[2]+kap*(h[0]*(y-x)**2+h[1]*y*y+h[2]*x*x)
    if transport=='endpoint':
        rd=r0*r0+(1-r0)**2;r=r0*r0/rd;rc=(1-r0)**2/rd;dr=2*r0*(1-r0)/rd**2
    else:r,rc,dr=r0,1-r0,np.ones_like(r0)
    hs=(r,rc*t,rc*(1-t));q=hs[0]*hs[1]+hs[0]*hs[2]+hs[1]*hs[2]
    xx,yy=c/u-hs[0]-hs[2],c/u-hs[0];n=u*((xx-yy/2)**2+3*yy*yy/4);d=det/u
    return {'Q':q,'u':u,'N':n,'d':d,'jacobian':v/(1+v+s*v)**3*rc*dr}


def evaluate(data):
    q,u,n,d=[np.asarray(data[k],dtype=np.float64) for k in ('Q','u','N','d')]
    a=1+q;x=(1+d)/(a*u);h=n/(a*u);plus=(1+x+h+np.sqrt((x-1)**2+2*h*(1+x)+h*h))/2
    minus=x/plus;la=np.log1p(q)
    parts=np.array([la*la/2,(np.log(u)+1.5)*la,-K(plus),-K(minus),K(1/u),np.full_like(q,np.pi**2/6-1.25)])
    return {'value':np.sum(parts,axis=0)*data['jacobian'],'term_scale':np.sum(np.abs(parts),axis=0)*data['jacobian']}
