"""Independent local Taylor evaluation of A'=B/(1+t)^2, with analytic tails.

The bounds concern exact-real truncation. Floating arithmetic and the closed
anchor are independently checked, not asserted to be interval certified.
"""
from fractions import Fraction as Q
from functools import lru_cache
import math
import mpmath as mp
import nested_analytic_marginals_v1 as M

CENTERS=tuple(Q(2)**j for j in range(-3,5))


@lru_cache(maxsize=16)
def small_coefficients(count):
    assert type(count) is int and 1<=count<=8192
    b=Q(1,60);d1=d2=Q();out=[]
    for n in range(count):
        if n:b=-Q(n,4*n+10)*b
        d=b-2*d1-d2;out.append(d/(n+3));d2,d1=d1,d
    return tuple(out)


def bubble_taylor(center,count,b0,b1,number):
    """ODE recurrence; works for exact rationals or high precision real values."""
    c=number(center);a0=c*c*(c+4);a1=3*c*c+8*c;a2=3*c+4
    p0=-c*(c+2);p1=-2*c-2;r0=c-4
    rows=[b0,b1]
    for n in range(count-2):
        value=(a1*n*(n+1)+p0*(n+1))*rows[n+1]
        value+=(a2*n*(n-1)+p1*n+r0)*rows[n]
        if n:value+=((n-1)*(n-2)-(n-1)+1)*rows[n-1]
        rows.append(-value/(a0*(n+2)*(n+1)))
    return rows[:count]


def derivative_coefficients(center,bubble,number):
    c=number(center);d1=d2=number(0);out=[]
    for b in bubble:
        d=(b-2*(1+c)*d1-d2)/(1+c)**2
        out.append(d);d2,d1=d1,d
    return out


def tail(ctx,center,h,count):
    c=ctx.mpf(center);h=abs(ctx.mpf(h));radius=c+ctx.mpf(1)/2
    assert c>0 and h<=c/3
    rho=h/radius
    maximum=ctx.mpf(8)/105*(2*c+ctx.mpf(1)/2)**2
    return maximum*h*rho**count/((count+1)*(1-rho))


@lru_cache(maxsize=32)
def prepared(center,digits):
    assert center in CENTERS and type(digits) is int and 50<=digits<=4000
    ctx=mp.mp.clone();ctx.dps=digits+60
    number=lambda q:ctx.mpf(q.numerator)/q.denominator if isinstance(q,Q) else ctx.mpf(q)
    c=number(center);count=math.ceil((digits+55)/math.log10(3))
    assert count<=10000
    h=ctx.asinh(ctx.sqrt(c)/2);pref=(c+4)**ctx.mpf('1.5')/(3*ctx.sqrt(c))
    b0=pref*h-4*c/9-ctx.mpf(4)/3
    b1=pref*(c-2)/(c*(c+4))*h+(c+4)/(6*c)-ctx.mpf(4)/9
    bubble=bubble_taylor(center,count,b0,b1,number)
    derivative=derivative_coefficients(center,bubble,number)
    coefficients=tuple(d/(n+1) for n,d in enumerate(derivative))
    anchor=M.bubble_primitive(ctx,c)
    maximum_tail=tail(ctx,c,c/3,count)
    assert maximum_tail<ctx.power(10,-digits-40)
    return ctx,coefficients,anchor


def evaluate(ctx,t):
    t=ctx.mpf(t);assert ctx.isfinite(t) and 0<=t<16 and 50<=ctx.dps<=4000
    if not t:return ctx.mpf(0),{'method':'small','terms':0,'model_absolute_tail':ctx.mpf(0)}
    if t<=ctx.mpf(1)/8:
        with ctx.extradps(60):
            lower=t**3/(180*(1+t/4)*(1+t)**2)
            target=ctx.power(10,-ctx.dps+20)*lower
            count=1
            while t**(count+3)/(45*(1-t))>target:
                count+=1;assert count<=8192
            value=ctx.mpf(0)
            for q in reversed(small_coefficients(count)):value=value*t+ctx.mpf(q.numerator)/q.denominator
            value*=t**3;bound=t**(count+3)/(45*(1-t))
            assert value>0
        return +value,{'method':'small','terms':count,'model_absolute_tail':+bound,
            'tail_excludes_rounding_and_anchor_reference_error':True}
    sign,mantissa,exponent,bits=t._mpf_
    exact_input=Q((-1 if sign else 1)*int(mantissa))*Q(2)**int(exponent)
    chosen=[c for c in CENTERS if abs(exact_input-c)<=c/3]
    # The exact geometric intervals cover [1/12,64/3]. Select using the
    # actual dyadic input; never clip or perturb it at chart boundaries.
    assert chosen
    center=min(chosen,key=lambda c:abs(exact_input-c))
    work,coefficients,anchor=prepared(center,ctx.dps)
    c=work.mpf(center.numerator)/center.denominator;h=work.mpf(t._mpf_)-c
    value=work.mpf(0)
    for coefficient in reversed(coefficients):value=value*h+coefficient
    bound=tail(work,c,h,len(coefficients));assert bound<work.power(10,-ctx.dps-40)
    return ctx.mpf((anchor+h*value)._mpf_),{'method':'local','center':str(center),'terms':len(coefficients),
        'model_absolute_tail':ctx.mpf(bound._mpf_),'tail_excludes_rounding_and_anchor_reference_error':True}
