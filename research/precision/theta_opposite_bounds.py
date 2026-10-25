"""Absolute radial tails for a shifted f1 on the first K2,3 path."""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial as fac,comb
from theta_core_radial_series import integrated_taylor_remainder


@lru_cache(None)
def S(n):return sum((Q(i+1,4**(n-i)) for i in range(n+1)),Q()) if n>=0 else Q()


def bounds(terms,n=20,lower_power=-12,upper_power=64):
    eps=Q(2)**lower_power;R=Q(2)**upper_power;low=high=Q()
    for row in terms:
        j,r,k,m=[row[z] for z in ('j','r','b','m')];c=abs(row['coefficient']);beta=Q(fac(m)**2,fac(2*m+1))
        low+=c*beta/Q(7200)*(S(r)*eps**(j+7)/(j+7)+2*S(r-1)*eps**(j+6)/(j+6)+S(r-2)*eps**(j+5)/(j+5))
        if k==m:factor=Q(1);power=k;logs=2
        elif k==m+1:factor=Q(2**(m+2));power=k;logs=3
        else:factor=Q(2**(m+2)*fac(m)*fac(k-m-2),fac(k-1));power=m+1;logs=2
        decay=r+power-j;assert decay>=2
        integral=R**(1-decay)*sum(Q(comb(logs,i)*fac(i)*(upper_power+1)**(logs-i),(decay-1)**(i+1)) for i in range(logs+1))
        high+=c*(2*r*r+2*r+1)*factor*integral/72
    series=integrated_taylor_remainder(n,2);q=low+high+series
    assert max(q.numerator.bit_length(),q.denominator.bit_length())<1000000
    total=Q(-(-q.numerator*2**160//q.denominator),2**160)
    return {'order':n,'lower_endpoint':low,'upper_endpoint':high,'whole_axis_series_remainder':series,'total':total,'log2_domain':(lower_power,upper_power)}
