"""Exact one-radial-branch K2,3 external Taylor coefficients.

The input is topology, not a numerical target.  No original forest column is
used as a second construction.  This module integrates Feynman parameters
and D=6 spherical moments algebraically; the radial integral is still owed.
"""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as Q
from math import factorial as fac, comb


def compositions(total, count):
    if count == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first,count-1):
                yield (first,*rest)


@lru_cache(None)
def symmetric_power_sum(n):
    # y^n+(1-y)^n as a polynomial in a=y(1-y).
    if n==0:return (2,)
    if n==1:return (1,)
    p=list(symmetric_power_sum(n-1));q=symmetric_power_sum(n-2)
    p += [0]*max(0,len(q)+1-len(p))
    for i,c in enumerate(q):p[i+1]-=c
    return tuple(p)


@lru_cache(None)
def uv_polynomial(a,b):
    # Coefficients of u^k v^(a+2b-k) in (u+v)^a(u^2+v^2-uv)^b.
    if b:
        prior=uv_polynomial(a,b-1);out=[0]*(len(prior)+2)
        for i,c in enumerate(prior):
            out[i]+=c;out[i+1]-=c;out[i+2]+=c
        return tuple(out)
    return tuple(comb(a,k) for k in range(a+1))


@lru_cache(None)
def beta_integer(a,b):
    assert a>0 and b>0
    return Q(fac(a-1)*fac(b-1),fac(a+b-1))


@lru_cache(None)
def uv_moment(a,b,l,m):
    if l>m:return uv_moment(a,b,m,l)
    degree=a+2*b
    answer=sum((c*beta_integer(k+l+1,l+1)*beta_integer(degree-k+m+1,m+1)
                for k,c in enumerate(uv_polynomial(a,b))),Q())
    assert answer>0
    return answer


def coefficient(n):
    """R_n(t)=sum c*t^j/(1+t)^(r+1)*J_(m,b)(t).

    J_(m,b)=int_0^1 [y(1-y)]^m/[1+t*y(1-y)]^b dy.
    The amplitude is int_0^infinity t^2 f(t) R_n(t)/2 dt.
    R_0 is deliberately not a renormalized integral.
    """
    assert isinstance(n,int) and 0<=n<=24
    out=defaultdict(Q)
    for a0,a1,a2,d,l,m,v in compositions(n,7):
        j=a0+a1+a2;k=d+l+m+v;r=2*a0+a1+d;s=a1+2*a2+l+m+v
        assert r+s==n+j
        base=Q((-1)**k*fac(r)*fac(s),fac(j+2)*fac(a0)*fac(a1)*fac(a2)*fac(d)*fac(l)*fac(m)*fac(v))
        base*=uv_moment(a1,a2+v,l,m)
        power=1+a1+2*a2+v+min(l,m)
        for h,c in enumerate(symmetric_power_sum(abs(l-m))):
            if c:
                key=(j,r,s+1,power+h)
                assert 1<=key[3]<=key[2]
                out[key]+=base*c
    return {key:value for key,value in sorted(out.items()) if value}


def gaussian_moment_ratios(bubbles):
    m=2*bubbles
    return [Q(comb(m,j)*fac(2),fac(j+2)) for j in range(m+1)]


def vacuum_derivative_bound(bubbles):
    # Integral of the positive Gaussian vacuum weight with one response S.
    # B(t)<=t^2/60 on [0,1]; B(t)<=t*log(2t)/6 on [1,infinity).
    b=bubbles
    low=(Q(b+3,2*b+3)+Q(b,2*b+2))/Q(12*60**b)
    logmoment=lambda n:sum(comb(n,k)*fac(k) for k in range(n+1))
    high=Q((b+1)*logmoment(b)+4*logmoment(b+1),6**b)
    return low+high


def integrated_taylor_remainder(n,bubbles):
    assert bubbles in (1,2) and n>=2*bubbles
    factors=gaussian_moment_ratios(bubbles)
    out=Q(1,(n+1)*4**(n+1))
    for j in range(1,2*bubbles+1):
        out+=factors[j]*fac(j-1)*comb(n,j-1)/4**(n-j+1)
    return out*vacuum_derivative_bound(bubbles)
