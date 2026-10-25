"""Original host products, Feynman simplex variance and hypergeometric bubbles.

No native, vector-diagnostic or production-runtime module is imported here.
The Euler representation of 2F1 binds B to its independently proved Stieltjes
measure: t^2/60 * 2F1(1,1;7/2;-t/4).
"""


def scalar(mp,point,mode,control):
    s,z,theta,x0,y,u=point
    t=(s/(1-s))**2;A=1+t;x=x0*x0/(x0*x0+(1-x0)**2);jac=2*x0*(1-x0)/(x0*x0+(1-x0)**2)**2
    qp1=mp.sqrt(t*z)*mp.cospi(theta);qp2=mp.sqrt(t*z)*(-mp.cospi(theta)+mp.sqrt(3)*mp.sinpi(theta))/2
    weights=[x,(1-x)*y,(1-x)*(1-y)];alpha=x*(1-x);V0=alpha*t
    def bubble(v):return v*v/60*mp.hyp2f1(1,1,mp.mpf(7)/2,-v/4)
    def derivative(v):return v/30*mp.hyp2f1(1,1,mp.mpf(7)/2,-v/4)-v*v/840*mp.hyp2f1(2,2,mp.mpf(9)/2,-v/4)
    if mode=='primitive':
        ell=mp.log1p(V0);arg=mp.exp(u*ell)/x-1;K0=-ell*bubble(arg)/(1+arg)
    else:
        b=1-x+V0;ell=mp.log1p(b/(4*x));a=x/b*mp.expm1(u*ell)
        weight=(1-4*a)**mp.mpf('1.5')/6*x/b*ell*mp.exp(u*ell)
        cm=(a-2)/(1-a)**2;cd=1/(1-a);ce=1/(1-a)**2;mass=1-x+x/a
        K0=weight*(-cm*mp.log1p(V0)-cd*x*V0/(1+V0)-ce/a*mp.log1p(V0/mass))
    full=mp.mpf(0);constant=weights[0]*weights[1]+weights[0]*weights[2]+weights[1]*weights[2]
    linear=2*weights[0]*(weights[1]*qp2-weights[2]*qp1)
    for sign in (-1,1):
        variance=weights[0]*weights[1]*(t+2*sign*qp2+1)+weights[0]*weights[2]*(t-2*sign*qp1+1)+weights[1]*weights[2]
        host=1/(A*(A-2*sign*qp1+1)*(A+2*sign*qp2+1))
        if mode=='primitive':
            shift=mp.log((1+variance)/(1+V0));arg=mp.exp(mp.log1p(V0)+u*shift)/x-1
            difference=-shift*bubble(arg)/(1+arg)
        else:
            difference=weight*(-cm*mp.log((1+variance)/(1+V0))-cd*x*(variance/(1+variance)-V0/(1+V0))-ce/a*mp.log((mass+variance)/(mass+V0)))
        full+=((host-1/A**3)*K0+host*difference)/2
    radial=t*t/2*(2*t/(s*(1-s)));simplex=(1-x)*jac
    value=radial*2*(1-z)*simplex*full
    if control:
        arg=(1-x+V0)/x;B=bubble(arg);first=-x*B/(1+V0)**2;second=2*x*B/(1+V0)**3-derivative(arg)/(1+V0)**2
        point=(4*(qp1*qp1+qp2*qp2-qp1*qp2)/A**2-2/A)*K0+constant*first+linear**2/2*second-2*(qp2-qp1)*linear/A*first
        average=(5*t/(3*A*A)-2/A)*K0+(alpha+(1-x)**2/6-t*alpha/A)*first+5*t*alpha**2/18*second
        value+=radial*simplex/A**3*(average-2*(1-z)*point)
    return value


def evaluate(mp,point,mode,orbit,control):
    p=[mp.mpf(q.numerator)/q.denominator for q in point];total=mp.mpf(0)
    for i in range(orbit):
        q=list(p);q[2]=(q[2]+mp.mpf(i)/orbit)%1
        total+=scalar(mp,q,mode,control)/(2*orbit);q[5]=1-q[5]
        total+=scalar(mp,q,mode,control)/(2*orbit)
    return total
