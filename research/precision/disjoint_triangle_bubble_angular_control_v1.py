"""Analytic angular control of the full mixed triangle/bubble integrand."""
from fractions import Fraction as Q
import numpy as np
import sympy as S
import disjoint_triangle_bubble_qmc_vector_v1 as BASE


def dot(a,b):return a[0]*b[0]+a[1]*b[1]-(a[0]*b[1]+a[1]*b[0])/2


def verify(binding):
    eps,A,M,C,L,l1,l2=S.symbols('eps A M C L l1 l2',positive=True)
    hp=1/(A*A*(A+eps*eps+eps*l1)*(A+eps*eps+eps*l2))
    hm=hp.subs(eps,-eps)
    paired=(hp*S.log(M+eps*eps*C+eps*L)+hm*S.log(M+eps*eps*C-eps*L))/2-S.log(M)/A**4
    actual=S.diff(paired,eps,2).subs(eps,0)/2
    expected=((l1*l1+l2*l2+l1*l2)/A**2-2/A)*S.log(M)+C/M-L*L/(2*M*M)-(l1+l2)*L/(A*M)
    assert S.simplify(actual-expected/A**4)==0
    t,y,v,q1,q2=S.symbols('t y v q1 q2');xs=[y,(1-y)*v,(1-y)*(1-v)]
    cv=lv=S.Integer(0);sv={};norms={}
    for pair in binding['triangle_pairs']:
        i,j=pair['simplex_pair'];c,a,b=map(S.Rational,pair['momentum'])
        norm=a*a-a*b+b*b;cv+=xs[i]*xs[j]*norm;lv+=xs[i]*xs[j]*2*c*(a*q1+b*q2)
        if i==0:sv[j]=(c*a,c*b);norms[j]=norm
    host=[tuple(map(S.Rational,row)) for row in binding['host_external_shifts']]
    assert all(dot(row,row)==1 for row in host) and set(sv)=={1,2}
    lhost=[2*(a*q1+b*q2) for a,b in host];alpha=y*(1-y)
    cbar=S.integrate(cv,(v,0,1));lbar=S.integrate(lv,(v,0,1));l2bar=S.integrate(S.expand(lv*lv),(v,0,1))
    assert S.simplify(cbar-alpha*(norms[1]+norms[2])/2-(1-y)**2/6)==0
    integrated=((lhost[0]**2+lhost[1]**2+lhost[0]*lhost[1])/A**2-2/A)*S.log(M)+cbar/M-l2bar/(2*M*M)-(lhost[0]+lhost[1])*lbar/(A*M)
    def average(poly):
        pp=S.Poly(S.expand(poly),q1,q2);moments={(0,0):1,(1,0):0,(0,1):0,(2,0):t/6,(1,1):-t/12,(0,2):t/6}
        assert all(powers in moments for powers,_ in pp.terms())
        return sum(co*moments[powers] for powers,co in pp.terms())
    sh=dot(host[0],host[0])+dot(host[1],host[1])+dot(host[0],host[1])
    ss=dot(sv[1],sv[1])+dot(sv[2],sv[2])+dot(sv[1],sv[2])
    cross=dot(tuple(a+b for a,b in zip(*host)),tuple(a+b for a,b in zip(sv[1],sv[2])))
    desired=(2*t*sh/(3*A*A)-2/A)*S.log(M)+cbar/M-t*alpha**2*ss/(9*M*M)-t*alpha*cross/(3*A*M)
    assert S.simplify(average(integrated)-desired)==0
    return {'paired_external_degree_two_exact':True,'simplex_C_L_L2_integrals_exact':True,
        'six_dimensional_tensor_projection_exact':True,'host_norm_sum':str(sh),'child_norm_sum':str(ss),'mixed_dot':str(cross),
        'conditional_control':'full density minus point control plus angular average/2',
        'angular_pdf':'2*(1-z) dz dtheta, theta in [0,1]','full_integrand_retained':True}


def evaluate(binding,points,chart,radial_power=2):
    result=BASE.evaluate(binding,points,chart,radial_power)
    s,z,theta,u=np.asarray(points,dtype=np.longdouble).T;t=(s/(1-s))**radial_power
    if chart==0:
        ell=np.log1p(t/2);y=np.expm1(u*ell)/t;jac=ell*np.exp(u*ell)/t
    else:y=(1+u)/2;jac=np.full_like(y,np.longdouble(1)/2)
    co,si=np.cos(np.pi*theta),np.sin(np.pi*theta);second=(-co+np.sqrt(np.longdouble(3))*si)/2
    vectors={};norms={}
    for row in binding['triangle_pairs']:
        i,j=row['simplex_pair'];c,a,b=map(int,row['momentum'])
        if i==0:vectors[j]=(c*a,c*b);norms[j]=a*a-a*b+b*b
    host=[tuple(map(int,row)) for row in binding['host_external_shifts']]
    h=host[0][0]*co+host[0][1]*second;k=host[1][0]*co+host[1][1]*second
    p1=vectors[1][0]*co+vectors[1][1]*second;p2=vectors[2][0]*co+vectors[2][1]*second
    alpha=y*(1-y);M=1+t*alpha;A=1+t;cbar=alpha*(norms[1]+norms[2])/2+(1-y)**2/6
    point=(4*t*z*(h*h+k*k+h*k)/(A*A)-2/A)*np.log1p(t*alpha)+cbar/M
    point-=2*t*z*alpha*alpha*(p1*p1+p2*p2+p1*p2)/(3*M*M)
    point-=2*t*z*alpha*(h+k)*(p1+p2)/(A*M)
    sh=dot(host[0],host[0])+dot(host[1],host[1])+dot(host[0],host[1])
    ss=dot(vectors[1],vectors[1])+dot(vectors[2],vectors[2])+dot(vectors[1],vectors[2])
    cross=dot(tuple(a+b for a,b in zip(*host)),tuple(a+b for a,b in zip(vectors[1],vectors[2])))
    average=(2*t*sh/(3*A*A)-2/A)*np.log1p(t*alpha)+cbar/M-t*alpha*alpha*ss/(9*M*M)-t*alpha*cross/(3*A*M)
    common=radial_power/(s*(1-s))*(1-y)*jac*(t/A)**3*(BASE.bubble(t)/A)
    return result+common*((1-z)*point-average/2)
