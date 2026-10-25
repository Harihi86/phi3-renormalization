"""Exact log-domain BPH overlap calculation with analytic counterterm moments.

The original two-sector kernel is kept in overlap.py. Here the counterterm
has a separate exact radial primitive, and every common large logarithm is
integrated analytically. Finite corrections retain the same Taylor scheme.
No massless finite constant or numerical fit replaces the massive amplitude.
"""
import math
from functools import lru_cache
from .finite import _gauss_legendre
from .primitive import _smooth_nodes as nodes
from .nested_triangle_extreme import _li2_unit, _log_f, _log_F

def delta_f(logz):
    if logz>=0:
        z=math.exp(-logz);s=math.log1p(z)
        return s+z*(logz+s)
    return _log_f(logz)-logz+1

def delta_F(logz):
    if logz>=0:
        z=math.exp(-logz);s=math.log1p(z)
        return -_li2_unit(z/(1+z))-.5*s*s-delta_f(logz)
    return _log_F(logz)-(.5*logz*logz-logz+1+math.pi**2/6)

@lru_cache(maxsize=64)
def overlap_scaled(log_ratio, order):
    L,n=log_ratio,order
    from .overlap import _h
    nd=nodes(n);phase=[(u*(1-u),v*(1-v),(u-v)**2,wu*wv) for u,wu in nd for v,wv in nd]
    raw=[]
    for x,wx in nd:
        for w,ww in nd:
            r=1-w;c=1+x*r;iv=w/(x+w);pref=2*w/c**3;df=r*(3+3*x*r+x*x*r*r)
            inner=[]
            for a,b,dist,weight in phase:
                T=a+r*dist/c;offset=math.log(b*iv);logz=L+offset
                if logz>=0:
                    low=math.log1p(math.exp(-logz));ell=offset+low;fr=offset-1+delta_f(logz);logistic=1/(1+math.exp(-logz))
                else:
                    ell=math.log1p(math.exp(logz))-L;fr=_log_f(logz)-L;logistic=math.exp(logz)/(1+math.exp(logz))
                alpha=x*T/b*logistic
                inner.append(weight*(T*(ell+alpha*_h(alpha))-df*b*fr))
            raw.append(wx*ww*pref*math.fsum(inner))
    raw= L/36+math.fsum(raw)
    base=-L*L/12+19*L/36+math.pi**2/72-137/108
    main=math.fsum(w*u*(1-u)*delta_F(L+math.log(u*(1-u)/2)) for u,w in nd)
    edge=L+48;rest=[]
    for j in range(math.ceil(edge/4)):
        left=4*j;right=min(edge,left+4)
        for x,w in _gauss_legendre(n):
            t=left+(right-left)*x;v=.5*math.exp(-t);d=1-v;g=(d*d+d-1)/(d*d*d)
            inner=math.fsum(wu*u*(1-u)*delta_f(L+math.log(u*(1-u)/2)-t) for u,wu in nd)
            rest.append(w*(right-left)*v*g*inner)
    return raw+base-main-math.fsum(rest)

