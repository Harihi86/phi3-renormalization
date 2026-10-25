"""Exact outward dyadic versions of the proved Fourier endpoint estimates."""
from fractions import Fraction as Q
from math import factorial
from exact_exp2_enclosures import log2
def up(q,bits=160):
    q=Q(q);assert q>=0;s=1<<bits;return Q((q.numerator*s+q.denominator-1)//q.denominator,s)
def tails():
    epsilon=Q(1,262144);R=Q(12);a=Q(17,18);b=Q(16,27);d=a+18*b*log2()[1]
    low=up(epsilon**2*(d*d+b*d+b*b/2)/48)
    term=total=Q(1)
    for n in range(1,385):term*=Q(48,n);total+=term
    exp_upper=1/total
    high=up(3*exp_upper*sum((Q(factorial(5),factorial(5-j))*R**(5-j)/4**(j+1) for j in range(6)),Q()))
    return {'left':epsilon,'right':R,'lower_endpoint_error':low,'upper_endpoint_error':high,'total':low+high}
