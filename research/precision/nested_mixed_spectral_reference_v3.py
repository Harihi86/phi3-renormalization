"""Independent exact rational endpoint substitution over the full Gaussian reference."""
from fractions import Fraction as Q
import nested_mixed_spectral_reference_v2 as B
M=B.M


def evaluate(point,binding,control=True,digits=300,mode='power2'):
    assert mode in ('power2','power4')
    n=2 if mode=='power2' else 4;p=list(point);p[7]=point[7]**n
    ctx,result=B.evaluate(p,binding,control,digits,mode='balanced')
    jac=Q(n)*point[7]**(n-1)
    return ctx,result*ctx.mpf(jac.numerator)/jac.denominator
