"""Full-domain six-coordinate diagnostics for the exact nested reduction.

The primitive and spectral kernels have different auxiliary coordinates.
They agree after that coordinate is integrated, not at each cube point.
"""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble


def evaluate(points,kind='primitive',squared_x=False,control=False,spectral_log=False):
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
        a=(1-u*u)/4;w=u**4/12
        if spectral_log:
            b=1-x+V0;ell=np.log1p(b/(4*x));a=x/b*np.expm1(u*ell)
            w=(1-4*a)**np.longdouble('1.5')/6*x/b*ell*np.exp(u*ell)
        d=1-a;base=x+a*(1-x)
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
    result=measure*total
    if control:
        q1=q*co;q2=q*p2;alpha=x*(1-x);argument=(1-x+V0)/x
        B=bubble(argument);first=-x*B/(1+V0)**2
        second=2*x*B/(1+V0)**3-bubble_derivative(argument)/(1+V0)**2
        point=(4*(q1*q1+q2*q2-q1*q2)/(A*A)-2/A)*K0+constant*first+linear*linear/2*second-2*(q2-q1)*linear/A*first
        average=(5*t/(3*A*A)-2/A)*K0+(alpha+(1-x)**2/6-t*alpha/A)*first+5*t*alpha*alpha/18*second
        common=2/(s*(1-s))*(t/A)**3*(1-x)*jac
        result+=common*(average/2-(1-z)*point)
    assert np.isfinite(result).all()
    return result


def orbit(points,kind='primitive',order=1,squared_x=False,control=False,spectral_log=False):
    assert order in (1,4,8)
    result=np.zeros(len(points),dtype=np.longdouble)
    for i in range(order):
        p=np.array(points,dtype=np.longdouble,copy=True);p[:,2]=(p[:,2]+np.longdouble(i)/order)%1
        result+=evaluate(p,kind,squared_x,control,spectral_log)/order
    return result


def bubble_derivative(t):
    result=np.empty_like(t);small=t<np.longdouble('.5')
    x=t[small];term=x*x/60;total=2*term/x
    for n in range(2,30):
        term=-term*x*(n-1)/(2*(2*n+3));total+=(n+1)*term/x
    result[small]=total
    x=t[~small];r=np.arcsinh(np.sqrt(x)/2)/np.sqrt(x/(x+4))
    result[~small]=((x-2)*r+2)/(3*x)-np.longdouble(5)/18
    assert (result>=0).all();return result
