"""Exact Schwinger and parameter elimination for the two radial core.

The simple SymPy implementation is an independent low-order oracle for a
future native implementation.  It never integrates divergent radial terms.
"""
from collections import defaultdict
from fractions import Fraction as Q
from math import factorial as fac
from functools import lru_cache
import sympy as s
import theta_two_radial as R

X,Z,Y,V,T,U,W=s.symbols('x z y v t u w')
VARIABLES=(X,Z,Y,V,T,U,W)


@lru_cache(None)
def generators(bi,bj):
    gram=s.Matrix([[1,-s.Rational(1,2)],[-s.Rational(1,2),1]])
    ai,aj,r0,r1=[s.Matrix([v]) for v in R.shifts(bi,bj)]
    dot=lambda a,b:(a*gram*b.T)[0]
    li=X*ai+Y*(1-V)*r0+Y*V*r1
    lj=Z*aj+Y*(1-V)*r0+Y*V*r1
    h=X*dot(ai,ai)+Z*dot(aj,aj)+Y*(1-V)*dot(r0,r0)+Y*V*dot(r1,r1)
    trace=T*dot(li,li)+U*dot(lj,lj)+(W-T-U)*dot(li,lj)
    determinant=(T*U-(W-T-U)**2/4)*(dot(li,li)*dot(lj,lj)-dot(li,lj)**2)
    return tuple(s.Poly(s.expand(v),*VARIABLES,domain=s.QQ) for v in (h,trace,determinant))


def raw_polynomial(n,bi,bj):
    h,trace,determinant=generators(bi,bj)
    value=s.Poly(0,*VARIABLES,domain=s.QQ)
    for j in range(n+1):
        block=s.Poly(0,*VARIABLES,domain=s.QQ)
        for k,c in enumerate(R.invariant_moment(j)):
            block += s.Rational(c.numerator,c.denominator)*trace**(j-2*k)*determinant**k
        value += s.Rational((-1)**(n-j),fac(n-j))*h**(n-j)*block
    return value


def eliminated_terms(n,bi,bj):
    """Return terms (xpower,zpower,ypower,tdegree,udegree,wdegree):QQ.

    Denominator powers are xpower+1,zpower+1,ypower+2; the last +2
    includes the y0,y1 change-of-variables Jacobian. No radial term is
    separately integrated by this representation.
    """
    out=defaultdict(Q)
    for (a,b,c,v,i,j,k),coefficient in raw_polynomial(n,bi,bj).terms():
        weight=Q(int(coefficient.p),int(coefficient.q))
        weight *= Q(fac(a)*fac(b)*fac(c+1),v+1)
        out[a,b,c,i,j,k]+=weight
    return {key:value for key,value in sorted(out.items()) if value}


def common_numerator(n,terms):
    grouped=defaultdict(lambda:defaultdict(Q))
    for (a,b,c,i,j,k),coefficient in terms.items():
        grouped[a,b,c][i,j,k]+=coefficient
    result=s.Poly(0,T,U,W,domain=s.QQ)
    powers={var:[s.Poly((1+var)**k,T,U,W,domain=s.QQ) for k in range(2*n+1)]
            for var in (T,U,W)}
    for (a,b,c),coefficients in grouped.items():
        assert max(a,b,c)<=2*n
        poly=s.Poly.from_dict({key:s.Rational(value.numerator,value.denominator)
                             for key,value in coefficients.items()},(T,U,W),domain=s.QQ)
        result += poly*powers[T][2*n-a]*powers[U][2*n-b]*powers[W][2*n-c]
    return result


def rational_at(terms,t,u,w):
    t,u,w=map(Q,(t,u,w))
    return sum((coefficient*t**i*u**j*w**k/
                ((1+t)**(a+1)*(1+u)**(b+1)*(1+w)**(c+2))
                for (a,b,c,i,j,k),coefficient in terms.items()),Q())
