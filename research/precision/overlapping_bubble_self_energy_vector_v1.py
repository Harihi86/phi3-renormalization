"""Complete five-coordinate weighted-kite diagnostic, six scale sectors."""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble


def density(a,b,x,y,u,control):
    U=a*b+a+b;U0=b*(a+1);U1=a*(b+1);M=1+a+b;d=b*(1-y)
    X=x*(1-x);Y=y*(1-y);Q0=b*Y;Q1=a*X;Q=Q1+Q0+a*b/U*(x-y)**2
    def J(q):
        mass=M+u*q;return d*q*q*(1-u)*bubble((mass-d)/d)/(mass*mass)
    result=J(Q)/U**3-J(Q0)/U0**3-J(Q1)/U1**3
    if control:
        A=a*a*(b+1)/U;B=a-2*a*b*y/U;C=b*Y+a*b*y*y/U
        mean=A*A/5-A*B/2+(B*B-2*A*C)/3+B*C+C*C
        point=Q*Q/U**3-Q0*Q0/U0**3-Q1*Q1/U1**3
        average=mean/U**3-Q0*Q0/U0**3-a*a/(30*U1**3)
        f=d*bubble((M-d)/d)/(M*M)
        result+=f*(average/2-(1-u)*point)
    return a*b*result


def evaluate(points,control=False,antithetic=False):
    p=np.asarray(points,dtype=np.longdouble);assert p.shape[1]==5 and ((p>0)&(p<1)).all()
    r,v,x,y,u=p.T;sectors=[(r*v,r,r),(r,1/v,1/v**2),(1/r,1/(r*v),1/(r**3*v**2))]
    result=np.zeros_like(r)
    for a,b,jac in sectors:
        for aa,bb in ((a,b),(b,a)):
            val=density(aa,bb,x,y,u,control)
            if antithetic:val=(val+density(aa,bb,x,y,1-u,control))/2
            result+=jac*val
    assert np.isfinite(result).all();return result
