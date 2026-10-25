"""Independent 720-order reference with a proved small-argument primitive series.

The exact coefficients still come from unreduced rational Gaussian elimination.
Only the scalar reference evaluator changes. No native kernel is imported.
"""
from fractions import Fraction as Q
import mpmath as mp
import shared_e3_full_sector_reference_v1 as O
import shared_bubble_square_series_reference_v1 as A
S=O.S
coefficients=O.coefficients


def evaluate_coefficients(exact,formulas,digits=640,details=None):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda q:ctx.mpf(q.numerator)/q.denominator
    assert sum(exact.values(),Q())==0 and all(type(t) is Q and t>0 for t in exact)
    values=[];small_constant=Q();small=large=0;maximum_terms=0;tail=[]
    for t,c in sorted(exact.items()):
        value=number(t)
        if t<=Q(1,16):
            primitive,meta=A.evaluate(ctx,value);small+=1;small_constant+=c
            maximum_terms=max(maximum_terms,meta['terms']);tail.append(abs(number(c))*meta['model_absolute_tail'])
        else:
            w=4/(ctx.sqrt(value+4)+ctx.sqrt(value))**2
            primitive=S.evaluate_w(ctx,formulas['0'],w);large+=1
        values.append(number(c)*primitive)
    # A(t)=P(w(t))-P(1). The exact sum of all coefficients is zero.
    if small_constant:values.append(number(small_constant)*S.evaluate_w(ctx,formulas['0'],ctx.mpf(1)))
    value=ctx.fsum(values);bound=ctx.fsum(tail)
    assert bound<ctx.mpf('1e-100')*abs(value)
    if details is not None:details.update({'series_arguments':small,'closed_arguments':large,'maximum_series_terms':maximum_terms,
        'exact_small_argument_coefficient_sum':small_constant,'normalization_evaluated':bool(small_constant),
        'model_relative_tail':ctx.nstr(bound/abs(value),35),'rounding_certificate_claimed':False})
    return ctx,value


def evaluate(point,binding,table,formulas,digits=640):
    return evaluate_coefficients(coefficients(point,binding,table),formulas,digits)
