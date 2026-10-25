"""Full mixed integrand minus its exact second-order angular control."""
import numpy as np
from disjoint_triangle_bubble_angular_control_v1 import BASE,dot

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
    return result+common*(1-z)*point


def orbit(binding,points,chart,order,radial_power=2):
    assert order in (1,4,8,16)
    value=np.zeros(len(points),dtype=np.longdouble)
    for i in range(order):
        p=np.array(points,dtype=np.longdouble,copy=True);p[:,2]=(p[:,2]+np.longdouble(i)/order)%1
        value+=evaluate(binding,p,chart,radial_power)/order
    return value
