"""Whole radial endpoint bounds, before any numerical integration."""
from fractions import Fraction as Q
from math import comb,factorial as fac
from theta_core_radial_series import beta_integer


def endpoints(rows,bubbles,lower_power=-20,upper_power=128):
    assert bubbles in (1,2) and lower_power<0 and upper_power>0
    eps=Q(1,2**(-lower_power));R=Q(2**upper_power);low=Q();high=Q()
    for row in rows:
        j,r,k,m=(row[key] for key in ('j','r','b','m'));c=abs(row['coefficient'])
        low+=c*beta_integer(m+1,m+1)*eps**(3+2*bubbles+j)/((3+2*bubbles+j)*2*60**bubbles)
        if k==m:
            constant=Q(1);power=k;extra_log=0
        elif k==m+1:
            constant=Q(2**(m+2));power=k;extra_log=1
        else:
            assert k>=m+2
            constant=2**(m+2)*beta_integer(m+1,k-m-1);power=m+1;extra_log=0
        decay=r+power-j;assert decay>=2,(row,decay)
        logs=bubbles+extra_log
        integral=R**(1-decay)*sum(Q(comb(logs,i)*fac(i)*(upper_power+1)**(logs-i),(decay-1)**(i+1)) for i in range(logs+1))
        high+=c*constant*integral/(2*6**bubbles)
    return {'lower':low,'upper':high,'lower_power':lower_power,'upper_power':upper_power}
