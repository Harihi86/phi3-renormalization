"""Stable full-integral diagnostics for the six-state overlapping core.

Every sector and Taylor term is retained. This vector floating implementation
is for diagnostics only; formal acceptance requires interval arithmetic.
"""
import numpy as np


def h(v):
    v=np.asarray(v,dtype=np.longdouble);out=np.empty_like(v);small=v<np.longdouble('.2')
    z=v[small];term=z*z/2;total=term.copy()
    for n in range(2,35):
        term=-term*z*(n-1)/(n+1);total+=term
    out[small]=total;z=v[~small];out[~small]=(1+z)*np.log1p(z)-z
    assert (out>=0).all();return out


def core(t,a,b,x,y,w,kind,log_aux):
    U=a*b+a+b;U0=b*(a+1);U1=a*(b+1);M=1+a+b
    X=x*(1-x);Y=y*(1-y);Q0=b*Y;Q1=a*X;d=a*X+a*b/U*(x-y)**2;Q=Q0+d
    eps=a/U0;den=3*eps+3*eps*eps+eps**3
    if kind=='logarithm':
        z0=t*Q0/M;delta=t*d/M
        increment=delta*np.log1p(z0)+(1+z0)*h(delta/(1+z0))
        bracket=increment-den*h(z0)-(U/U1)**3*h(t*Q1/M)
        return a*b*M/U**3*bracket
    assert kind=='remainder'
    u=w;jac=np.ones_like(w)
    if log_aux:
        scale=t*Q/M;ell=np.log1p(scale);u=np.expm1(w*ell)/scale;jac=ell*np.exp(w*ell)/scale
    c=u*t
    increment=d*(2*M*Q0+M*d+c*Q0*Q)/((M+c*Q0)*(M+c*Q))
    bracket=increment-den*Q0*Q0/(M+c*Q0)-(U/U1)**3*Q1*Q1/(M+c*Q1)
    return t*t*(1-u)*a*b/U**3*bracket*jac


def outer(t,nu,gid):
    A=1+t
    if gid=='L3E2G7':
        root=np.sqrt(t*t+4);delta=2/(t+root);e=A+delta
        return delta*(3*A*A+6*A*delta+(3+t)*delta*delta)/(3*A**3*e**3)
    assert gid=='L3E3G38'
    root=np.sqrt(t*t+4*t*nu*(1-nu)+4);delta=2*(t*nu*(1-nu)+1)/(root+t)
    return -delta*(2*A+delta)/(A*A*(A+delta)**2)


def evaluate(points,gid,kind='logarithm',orbit=False,log_aux=False,squared_a=False):
    pts=np.asarray(points,dtype=np.longdouble)
    assert pts.shape[1]==7 and ((pts>0)&(pts<1)).all()
    s,r,v,x,y,w,nu=pts.T;t=(s/(1-s))**2;A=1+t
    a=r/(1-r);da=1/(1-r)**2
    if squared_a:a=a*a;da=2*r/(1-r)**3
    b=a/v;sector=2*a/v**2*da
    pairs=[(x,y)] if not orbit else [(x,y),(1-x,y),(y,x),(1-y,x)]
    density=sum((core(t,a,b,xx,yy,w,kind,log_aux) for xx,yy in pairs))/len(pairs)
    measure=t*t/(2*A*A)*2*s/(1-s)**3*sector
    result=measure*outer(t,nu,gid)*density
    assert np.isfinite(result).all();return result
