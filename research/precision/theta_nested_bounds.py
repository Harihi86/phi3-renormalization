"""Signed nested spectral majorant for the one-radial-branch K2,3 core."""
from fractions import Fraction as Q
from math import factorial as fac,comb


def up(q):
    assert q>=0 and max(q.numerator.bit_length(),q.denominator.bit_length())<1000000
    return Q(-(-q.numerator*2**160//q.denominator),2**160)


def bounds(terms,n=20):
    assert n>=2
    J=Q(43)+Q(7,1536)
    series=J*Q(1,4**n)*(Q(1,4*(n+1))+Q(2,3)+Q(n,3))
    low=high=Q();eps=Q(1,1024);R=Q(2**49)
    for row in terms:
        j,r,k,m=[row[z] for z in ('j','r','b','m')];c=abs(row['coefficient'])
        low+=c*Q(fac(m)**2,fac(2*m+1))*Q(5,192)*eps**(5+j)/(5+j)
        if k==m:constant=Q(1);power=k;logs=2
        elif k==m+1:constant=Q(2**(m+2));power=k;logs=3
        else:constant=Q(2**(m+2)*fac(m)*fac(k-m-2),fac(k-1));power=m+1;logs=2
        decay=r+power-j;assert decay>=2
        integral=R**(1-decay)*sum(Q(comb(logs,i)*fac(i)*51**(logs-i),(decay-1)**(i+1)) for i in range(logs+1))
        high+=c*constant*integral/8
    return {'order':n,'lower_endpoint':low,'upper_endpoint':high,'whole_axis_series_remainder':series,'positive_absolute_spectral_Jstar_bound':J,'total':up(low+high+series),'domain':(eps,R)}
