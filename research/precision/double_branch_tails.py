"""Exact absolute endpoint bounds for the signed F1*F2 integral."""
from fractions import Fraction as Q
from math import factorial
from exact_exp2_enclosures import log2
from fourier_endpoint_bounds import up

def tails():
    eps,R=Q(1,262144),Q(12);L=18*log2()[1]
    a=Q(17,18);b=Q(16,27);A=Q(64,243)*Q(54,7)+Q(2,135);B=Q(64,243)*Q(110,21)+Q(1,324);C=Q(64,243)
    coefficients=(a*A,a*B+b*A,a*C+b*B,b*C)
    low=Q()
    for j in range(4):
        low+=sum((coefficients[k]*Q(factorial(k),factorial(k-j))*L**(k-j) for k in range(j,4)),Q())/2**j
    low=up(eps**2*low/48);term=exp48=Q(1)
    for n in range(1,385):term*=Q(48,n);exp48+=term
    high=up(5/exp48*sum((Q(factorial(6),factorial(6-j))*R**(6-j)/4**(j+1) for j in range(7)),Q()))
    return {'left':eps,'right':R,'lower_endpoint_absolute_error':low,'upper_endpoint_absolute_error':high,'total_absolute_error':low+high,'tails_must_expand_both_sides':True}
