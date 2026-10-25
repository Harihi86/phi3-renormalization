"""Large-argument Frobenius reference for A'=B/(1+t)^2.

The analytic bound controls series truncation, not floating arithmetic or the
closed-form anchor. The frozen native kernel and earlier references are intact.
"""
from fractions import Fraction as Q
from functools import lru_cache
import math
import mpmath as mp
import nested_analytic_marginals_v1 as M


@lru_cache(maxsize=8)
def exact_coefficients(nmax):
    assert type(nmax) is int and 1<=nmax<=8192
    a=Q(1,6);b=Q(-4,9);alpha=a;beta=b;alpha2=beta2=Q()
    bubble=[(a,b)];divided=[(alpha,beta)];primitive=[]
    for n in range(1,nmax+1):
        an=-Q(2*(2*n-5),n)*a
        bn=Q(2,n)*an-Q(2*(2*n-5),n)*b-Q(10-8*n,n*n)*a
        anew=an-2*alpha-alpha2;bnew=bn-2*beta-beta2
        primitive.append((anew/n,-anew/(n*n)-bnew/n))
        bubble.append((an,bn));divided.append((anew,bnew))
        a,b=an,bn;alpha2,alpha=alpha,anew;beta2,beta=beta,bnew
    return tuple(bubble),tuple(divided),tuple(primitive)


def tail(ctx,z,terms):
    assert 0<z<=ctx.mpf(1)/16
    ratio=4*z
    return (2*abs(ctx.log(z))+10)*ratio**(terms+1)/(1-ratio)


def asymptotic(ctx,z,coefficients):
    logz=ctx.log(z);p=q=ctx.mpf(0)
    for a,b in reversed(coefficients):p=p*z+a;q=q*z+b
    return logz*logz/12+ctx.mpf(4)*logz/9+z*(logz*p+q)


@lru_cache(maxsize=8)
def prepared(digits):
    assert type(digits) is int and 50<=digits<=4000
    ctx=mp.mp.clone();ctx.dps=digits+40
    terms=math.ceil((digits+50)/math.log10(4))
    assert terms<=8192
    exact=exact_coefficients(terms)[2]
    number=lambda q:ctx.mpf(q.numerator)/q.denominator
    coefficients=tuple((number(a),number(b)) for a,b in exact)
    anchor_z=ctx.mpf(1)/16;anchor_tail=tail(ctx,anchor_z,terms)
    assert anchor_tail<ctx.power(10,-digits-35)
    anchor=M.bubble_primitive(ctx,ctx.mpf(16))-asymptotic(ctx,anchor_z,coefficients)
    return ctx,coefficients,anchor,anchor_tail


def evaluate(ctx,t):
    t=ctx.mpf(t);assert ctx.isfinite(t) and t>=16 and 50<=ctx.dps<=4000
    work,coefficients,anchor,anchor_tail=prepared(ctx.dps)
    exact_input=work.mpf(t._mpf_);z=1/exact_input
    truncation=tail(work,z,len(coefficients))+anchor_tail
    assert truncation<work.power(10,-ctx.dps-30)
    value=asymptotic(work,z,coefficients)+anchor
    return ctx.mpf(value._mpf_),{'terms':len(coefficients),'model_absolute_tail':ctx.mpf(truncation._mpf_),
        'anchor_t':16,'anchor_digits':work.dps,'tail_excludes_rounding_and_anchor_reference_error':True}
