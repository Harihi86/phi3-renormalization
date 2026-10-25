"""Unreduced Gaussian reference with independently proved primitive series."""
from fractions import Fraction as Q
import mpmath as mp
import distinct_e3_full_sector_reference_v1 as H
import bubble_primitive_local_series_reference_v1 as LOCAL
import bubble_primitive_large_series_reference_v1 as LARGE

G,A,B=H.G,H.A,H.B
coefficients=H.coefficients


def evaluate_coefficients(pair,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda q:ctx.mpf(q.numerator)/q.denominator
    a,b=pair;values=[]
    assert sum(a.values(),Q())==0
    for t in sorted(set(a)|set(b)):
        assert t>0;x=number(t)
        if t in a:
            primitive=(LARGE if x>=16 else LOCAL).evaluate(ctx,x)[0]
            values.append(number(a[t])*primitive)
        if t in b:values.append(number(b[t])*B.hyper_bubble(ctx,x))
    return ctx,ctx.fsum(values)


def evaluate(point,binding,table,digits=640):
    return evaluate_coefficients(coefficients(point,binding,table),digits)
