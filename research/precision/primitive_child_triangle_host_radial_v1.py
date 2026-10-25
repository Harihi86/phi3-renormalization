"""Floating complete radial antiderivative; diagnostic use only."""
import numpy as np
from scipy.special import spence


def phi(t):
    t=np.asarray(t,dtype=np.longdouble);assert (t>0).all();z=t-1;near=abs(z)<np.longdouble('.1')
    result=np.empty_like(t);pi=np.arccos(np.longdouble(-1));q=t[~near]
    result[~near]=pi*pi/6-np.asarray(spence(np.asarray(q,dtype=float)),dtype=np.longdouble)+q*(2-3*q)*np.log(q)/(2*(1-q)**2)-q/(2*(1-q))
    v=z[near];value=np.full(len(v),pi*pi/6-np.longdouble(5)/4,dtype=np.longdouble)
    for k in range(1,25):value+=(-1)**(k-1)*2*v**k/(k*k*(k+1)*(k+2))
    result[near]=value;return result


def rational_moments(h):
    h=np.asarray(h,dtype=np.longdouble);z=h-1;near=abs(z)<np.longdouble('.1');out=[]
    current=np.empty_like(h);current[~near]=np.log(h[~near])/z[~near]
    for n in range(1,5):
        if n>1:current[~near]=(np.longdouble(1)/(n-1)-current[~near])/z[~near]
        current[near]=sum(((-z[near])**k/(n+k) for k in range(25)),np.zeros(near.sum(),dtype=np.longdouble))
        out.append(current.copy())
    return h*(out[0]+out[1]+out[2]),out[3]


def evaluate(q,ell,De,Ge,meanDe,control=True):
    a=1+q;b=1+ell*a+De;c=ell+Ge
    disc=(1-ell*a+De)**2+4*a*(ell*De-Ge);assert (disc>=0).all()
    root=np.sqrt(disc);plus=(b+root)/(2*a);minus=2*c/(b+root)
    la=np.log1p(q);pi=np.arccos(np.longdouble(-1))
    result=la*np.log(ell)+la*la/2+3*la/2-phi(plus/ell)-phi(minus/ell)+phi(1/ell)+pi*pi/6-np.longdouble(5)/4
    if control:
        A,B=rational_moments(ell);result+=(np.longdouble(1)/4-q)*A-(meanDe-De)*B
    assert np.isfinite(result).all();return result
