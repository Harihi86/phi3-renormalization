"""Choose argument-specific truncation; retain the proved full anchor series."""
import bubble_primitive_large_series_reference_v1 as FULL

M=FULL.M


def evaluate(ctx,t):
    t=ctx.mpf(t);assert ctx.isfinite(t) and t>=16 and 50<=ctx.dps<=4000
    work,coefficients,anchor,anchor_tail=FULL.prepared(ctx.dps)
    z=1/work.mpf(t._mpf_);ratio=4*z
    target=work.power(10,-ctx.dps-40)
    factor=(2*abs(work.log(z))+10)/(1-ratio)
    count=max(1,int(work.ceil(work.log(target/factor)/work.log(ratio)))-1)+2
    assert count<=len(coefficients)
    argument_tail=FULL.tail(work,z,count)
    while argument_tail>target:
        count+=1;assert count<=len(coefficients);argument_tail=FULL.tail(work,z,count)
    assert argument_tail<=target and argument_tail+anchor_tail<work.power(10,-ctx.dps-30)
    value=FULL.asymptotic(work,z,coefficients[:count])+anchor
    return ctx.mpf(value._mpf_),{'argument_terms':count,'anchor_terms':len(coefficients),
        'model_absolute_tail':ctx.mpf((argument_tail+anchor_tail)._mpf_),
        'argument_tail':ctx.mpf(argument_tail._mpf_),'anchor_t':16,
        'tail_excludes_rounding_and_anchor_reference_error':True}
