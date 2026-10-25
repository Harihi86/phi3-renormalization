"""Absolute rectangular radial tails of the complete truncated kernel."""
from fractions import Fraction as Q
from math import comb, factorial as fac


def ordinary_taylor_majorant(n):
    assert isinstance(n,int) and 1<=n<=24
    return Q(sum(comb(2*j+7,7) for j in range(1,n+1)),4)


def log_tail_polynomial(power,log_bound):
    return sum(Q(fac(power),fac(power-k))*log_bound**(power-k)
               for k in range(power+1))


def bounds(n,log2_lower=-12,log2_upper=80):
    assert -32<=log2_lower<=-1 and 4<=log2_upper<=128
    epsilon=Q(2)**log2_lower;upper=Q(2)**log2_upper;split=upper/4
    log_bound=Q(log2_upper-1)
    coefficient=ordinary_taylor_majorant(n)
    low=Q(101,5000)*coefficient*epsilon**5
    branch=(log_tail_polynomial(2,log_bound)
            +4*log_tail_polynomial(3,log_bound))/(72*split)
    central=16/upper*(Q(1,50400)+log_bound**3/216)
    central+=40/upper**2*(Q(1,50400)+split*log_bound**2/72)
    high=2*coefficient*(branch+central)
    return {'order':n,'log2_radial_domain':(log2_lower,log2_upper),
            'all_four_ordinary_propagator_bound':coefficient,
            'lower_radial_error':low,'upper_radial_error':high,'total':low+high}
