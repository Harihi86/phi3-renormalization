"""Exact Taylor recurrence and analytic truncation bound for the squared-bubble primitive.

A(t)=integral_0^t B(s)^2/(1+s)^3 ds. mpmath rounding is not certified here;
the stated tail bounds are exact-real bounds, separate from floating evaluation.
"""
from fractions import Fraction as Q
from functools import lru_cache
import math
import sympy as S


@lru_cache(maxsize=1)
def differential_operator():
    t=S.symbols('t',positive=True)
    p=-(t+2)/(t*(t+4));r=(t-4)/(t*t*(t+4))
    # B''+p B'+r B=0; the symmetric square annihilates B^2.
    square=[4*p*r+2*S.diff(r,t),2*p*p+S.diff(p,t)+4*r,3*p,S.Integer(1)]
    y,v=S.symbols('y v');jet=[y*y]
    for j in range(3):
        expr=jet[-1];jet.append(S.expand(S.diff(expr,t)+S.diff(expr,y)*v-S.diff(expr,v)*(p*v+r*y)))
    assert S.cancel(sum(a*b for a,b in zip(square,jet)))==0
    u=(1+t)**3
    transformed=[S.cancel(sum(square[j]*math.comb(j,k)*S.diff(u,t,j-k) for j in range(k,4))) for k in range(4)]
    denominator=S.lcm([S.denom(v) for v in transformed])
    polys=[S.Poly(S.cancel(denominator*v),t) for v in transformed]
    common=S.gcd_list([p.as_expr() for p in polys]);polys=[S.Poly(S.cancel(p.as_expr()/common),t) for p in polys]
    assert all(c.is_Integer for p in polys for c in p.all_coeffs())
    if polys[-1].LC()<0:polys=[-p for p in polys]
    records=tuple(tuple((int(powers[0]),int(c)) for powers,c in p.terms()) for p in polys)
    assert all(k>=j for j,rows in enumerate(records) for k,c in rows)
    return records


def falling(n,j):
    return math.prod(range(n-j+1,n+1)) if n>=j else 0


@lru_cache(maxsize=16)
def coefficients(count):
    assert type(count) is int and 1<=count<=8192
    operator=differential_operator();f=[Q()]*4+[Q(1,3600)]
    for n in range(5,count+4):
        leading=Q();previous=Q()
        for j,rows in enumerate(operator):
            for k,c in rows:
                index=n+j-k
                if index<j:continue
                factor=c*falling(index,j)
                if index==n:leading+=factor
                else:
                    assert index<n
                    previous+=factor*f[index]
        assert leading
        f.append(-previous/leading)
    out=tuple(f[n+4]/(n+5) for n in range(count))
    assert all(c*((-1)**n)>0 for n,c in enumerate(out))
    return out


def exact_identities():
    # F=2F1(1,1;7/2;-t/4); its coefficients have the defining hypergeometric ratio.
    n=S.symbols('n',integer=True,nonnegative=True)
    assert S.cancel(-(n+1)/(4*n+14)-(-S.Rational(1,4))*(n+1)/(n+S.Rational(7,2)))==0
    b=[Q(1,60)]
    for j in range(23):b.append(-b[-1]*Q(j+1,4*j+14))
    independent=[]
    for n in range(24):
        value=sum((b[i]*b[j]*(-1)**(n-i-j)*math.comb(n-i-j+2,2)
                   for i in range(n+1) for j in range(n-i+1)),Q())/(n+5)
        independent.append(value)
    assert tuple(independent)==coefficients(24)
    t=S.symbols('t',positive=True);N=S.symbols('N',integer=True,nonnegative=True)
    geometric=t**N/(1-t)
    assert S.cancel(t*S.diff(geometric,t)+2*geometric-t**N*((N+2)-(N+1)*t)/(1-t)**2)==0
    # b_n <= (1/60)4^-n, hence the absolute coefficient of t^(n+5)
    # in A is <= binomial(n+2,2)/(2025(n+5)) <= (n+2)/4050.
    assert Q(1,3600)/Q(3,4)**2==Q(1,2025)
    assert S.cancel((N+2)/S.Integer(4050)-(N+2)*(N+1)/(S.Integer(4050)*(N+5))-(4*(N+2)/(S.Integer(4050)*(N+5))))==0
    return {'operator':differential_operator(),'initial_A_coefficient':Q(1,18000),
        'independent_convolution_coefficients':24,'coefficient_recurrence_exact':True,
        'absolute_tail_after_N_terms':'t^(N+5)*((N+2)-(N+1)t)/(4050(1-t)^2)',
        'tail_domain':'0<=t<1','lower_bound':'t^5/(18000(1+t/4)^2(1+t)^3)',
        'lower_bound_reason':'Euler-positive bubble representation gives B(s)>=s^2/(60(1+s/4)); integrate s^4 using endpoint denominators',
        'runtime_domain':'0<=t<=1/16','tail_is_exact_real_only':True,'rounding_certificate_claimed':False}


def evaluate(ctx,t):
    t=ctx.mpf(t);assert 0<=t<=ctx.mpf(1)/16
    if not t:return ctx.mpf(0),{'terms':0,'model_absolute_tail':ctx.mpf(0)}
    target=ctx.power(10,-ctx.dps-20)
    with ctx.extradps(30):
        lower=t**5/(18000*(1+t/4)**2*(1+t)**3)
        power=t;count=1
        while True:
            tail=t**5*power*((count+2)-(count+1)*t)/(4050*(1-t)**2)
            if tail<=target*lower:break
            count+=1;assert count<=8192;power*=t
        value=ctx.mpf(0)
        for c in reversed(coefficients(count)):value=value*t+ctx.mpf(c.numerator)/c.denominator
        value*=t**5;assert value>0
        return value,{'terms':count,'model_absolute_tail':tail,'model_lower_bound':lower}
