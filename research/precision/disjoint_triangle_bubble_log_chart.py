"""Resolve the y~1/t simplex layer by an exact conditional coordinate map.

Lower half: y=expm1(u*log1p(t/2))/t, 0<=u<=1.
Upper half: y=(1+u)/2. Both Jacobians are retained without renormalization.
The fixed-order diagnostic still does not certify its finite cubature error.
"""
import math
import numpy as np
import sympy as s
import disjoint_triangle_bubble_pilot as BASE
recipe=BASE.recipe


def verify_map():
    t=s.symbols('t',positive=True);u=s.symbols('u',real=True);L=s.log(1+t/2)
    y=(s.exp(u*L)-1)/t;jac=L*s.exp(u*L)/t
    assert s.simplify(s.diff(y,u)-jac)==0 and s.simplify(y.subs(u,0))==0
    assert s.simplify(y.subs(u,1)-s.Rational(1,2))==0
    assert s.simplify(s.integrate(jac,(u,0,1))-s.Rational(1,2))==0
    assert s.diff((1+u)/2,u)==s.Rational(1,2)
    return {'lower_half_exact':[0,'1/2'],'upper_half_exact':['1/2',1],
            'lower_jacobian':'log1p(t/2)*exp(u*log1p(t/2))/t',
            'upper_jacobian':'1/2','both_positive_at_t_positive':True,
            'original_simplex_factor_1_minus_y_retained':True,
            'weights_not_renormalized':True,'finite_cubature_error_not_certified':True}


def prepare(rec,n):
    assert 2<=n<=24
    node,weight=np.polynomial.legendre.leggauss(n);node=(node+1)/2;weight=weight/2
    return rec,node,weight


def mapped_point(rec,t,z,theta,u,v,chart):
    assert chart in (0,1)
    if chart==0:
        L=math.log1p(t/2);y=math.expm1(u*L)/t;jac=L*math.exp(u*L)/t
    else:y=(1+u)/2;jac=.5
    return BASE.mapped(t,BASE.arrays(rec,z,theta,y,v))*jac


def tile(prepared,left):
    rec,node,weight=prepared;values=[]
    common=weight[:,None,None,None]*weight[None,:,None,None]*weight[None,None,None,:]
    for xx,wx in zip(node,weight):
        t=2.**(left+float(xx));L=math.log1p(t/2)
        ys=np.concatenate((np.expm1(L*node)/t,(1+node)/2))
        wy=np.concatenate((weight*(L/t)*np.exp(L*node),weight/2))
        cached=BASE.arrays(rec,node[:,None,None,None],node[None,:,None,None],ys[None,None,:,None],node[None,None,None,:])
        values.append(float(wx)*float(np.sum(BASE.mapped(t,cached)*common*wy[None,None,:,None],dtype=np.float64)))
    return math.fsum(values)
