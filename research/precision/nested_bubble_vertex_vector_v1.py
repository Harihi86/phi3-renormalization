"""Full-domain six-coordinate diagnostics for the exact nested reduction.

The primitive and spectral kernels have different auxiliary coordinates.
They agree after that coordinate is integrated, not at each cube point.
"""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble


def evaluate(points,kind='primitive',squared_x=False):
    assert kind in ('primitive','spectral')
    s,z,theta,x0,y,u=np.asarray(points,dtype=np.longdouble).T
    assert ((points>0)&(points<1)).all()
    t=(s/(1-s))**2;A=1+t;X=A+1
    if squared_x:
        den=x0*x0+(1-x0)**2;x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    else:x=x0;jac=np.ones_like(x)
    co,si=np.cos(np.pi*theta),np.sin(np.pi*theta)
    p2=(-co+np.sqrt(np.longdouble(3))*si)/2;q=np.sqrt(t*z)
    l1=-2*q*co;l2=2*q*p2
    V0=x*(1-x)*t;linear=2*x*(1-x)*q*(y*p2-(1-y)*co)
    constant=x*(1-x)+(1-x)**2*y*(1-y)
    if kind=='primitive':
        ell=np.log1p(V0);argument=(1-x+np.expm1(u*ell))/x
        K0=-ell*bubble(argument)/(1+argument)
    else:
        a=(1-u*u)/4;d=1-a;base=x+a*(1-x);w=u**4/12
        K0=w*((1+d)/d**2*np.log1p(V0)-x/d*V0/(1+V0)-np.log1p(a*V0/base)/(a*d**2))
    total=np.zeros_like(t)
    for sign in (-1,1):
        DV=constant+sign*linear;V=V0+DV
        assert (V>=0).all()
        denominator=((X+sign*l1)/A)*((X+sign*l2)/A)
        change=-(2/A+1/(A*A)+sign*(X/A)*(l1+l2)/A+(l1/A)*(l2/A))/denominator
        if kind=='primitive':
            shift=np.log1p(DV/(1+V0));new=(1+V0)*np.exp(u*shift)
            argument=(new-x)/x
            delta=-shift*bubble(argument)/(1+argument)
        else:
            delta=w*((1+d)/d**2*np.log1p(DV/(1+V0))-x/d*DV/((1+V0)*(1+V))-np.log1p(a*DV/(base+a*V0))/(a*d**2))
        total+=(change*K0+delta/denominator)/2
    measure=2/(s*(1-s))*(t/A)**3*(1-z)*(1-x)*jac
    result=measure*total;assert np.isfinite(result).all()
    return result


def orbit(points,kind='primitive',order=1,squared_x=False):
    assert order in (1,4,8)
    result=np.zeros(len(points),dtype=np.longdouble)
    for i in range(order):
        p=np.array(points,dtype=np.longdouble,copy=True);p[:,2]=(p[:,2]+np.longdouble(i)/order)%1
        result+=evaluate(p,kind,squared_x)/order
    return result
