"""Independent 1664-bit fixed-point logarithms for native mass validation."""
from fractions import Fraction as Q
from functools import lru_cache
SCALE=2**1664


def ceil_div(a,b):return -(-a//b)


@lru_cache(None)
def unit_log(r):
    assert 1<=r<=2
    z=(r-1)/(r+1);n=z.numerator*SCALE;d=z.denominator;zl=n//d;zh=ceil_div(n,d)
    al=zl*zl//SCALE;ah=ceil_div(zh*zh,SCALE);lo=hi=0;tl,th=zl,zh
    for j in range(640):
        lo+=2*tl//(2*j+1);hi+=ceil_div(2*th,2*j+1)
        tl=tl*al//SCALE;th=ceil_div(th*ah,SCALE)
    hi+=ceil_div(2*th*SCALE,1281*(SCALE-ah))
    return Q(lo,SCALE),Q(hi,SCALE)


@lru_cache(None)
def log_bounds(q):
    q=Q(q);assert q>0;n=q.numerator.bit_length()-q.denominator.bit_length()
    if q<Q(2)**n:n-=1
    r=q/Q(2)**n;assert 1<=r<2
    a,b=unit_log(r);c,d=unit_log(Q(2));shift=(n*c,n*d) if n>=0 else (n*d,n*c)
    return a+shift[0],b+shift[1]


def mass_interval(piece):
    l,h,ml,mr=[piece[k] for k in ('left','right','mid_left','mid_right')];slope=(mr-ml)/(h-l);intercept=ml-slope*l
    rational=slope/2*((h*h-l*l)/2-2*(h-l)+1/(1+h)-1/(1+l))+intercept/2*((h-l)-1/(1+h)+1/(1+l));co=(3*slope-2*intercept)/2
    a,b=log_bounds((1+h)/(1+l));values=(co*a,co*b);return rational+min(values),rational+max(values)
