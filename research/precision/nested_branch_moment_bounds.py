"""Exact outward errors for cumulative Green moments and cut constants."""
from fractions import Fraction as Q
from functools import lru_cache
from exact_exp2_enclosures import log2
from fourier_endpoint_bounds import up

def ceil(q):q=Q(q);return -(-q.numerator//q.denominator)
def negative_coarse(q):q=Q(q);return Q(3,8)**(q.numerator//q.denominator)

@lru_cache(None)
def negative_exp_upper(integer):
    assert integer>=0;term=total=Q(1)
    for j in range(1,769):term*=Q(integer,j);total+=term
    return 1/total

def weights(l,h):
    assert Q(1,262144)<=l<h<=12
    if h<=1:four=[h**3/6,h/2,1/(4*l**3),1/(4*l)]
    elif l>=1:
        e=negative_coarse(4*l);four=[6*h*h,6*h,36*h*e,36*h*h*e]
    else:
        e=Q(3)**ceil(2*h);four=[e/2,e/(2*l),1/(4*l**3),1/(4*l)]
    if h<=Q(1,2):two=[h**9/(10*l**4),h**9/(20*l**4)]
    elif l>=1:
        e=Q(3)**ceil(2*(h-l));two=[4*h**7/l**2*e,12*h**5/l**2*e]
    else:
        e=Q(3)**ceil(2*h);two=[h**7/(3*l**4)*e,h**5/l**4*e]
    return four+two

def errors(piece):
    return [up((piece['right']-piece['left'])*piece['C_error']*w) for w in weights(piece['left'],piece['right'])]

def tails():
    eps,R=Q(1,262144),Q(12);e4=negative_exp_upper(48);e8=negative_exp_upper(96)
    moments=[up(18*e4/R),up(18*e4/R**2),up(54*e8/R**2),up(54*e8/R)]
    common=eps**2*(Q(8,9)+Q(16,27)*(18*log2()[1]+Q(1,2)))
    constants_low=[up(common/20),up(common/40)]
    constants_high=[up(48*e4*(R*R/4+R/8+Q(1,32))),up(36*e4)]
    return {'left':eps,'right':R,'moment_positive_upper_tails':moments,'constant_positive_lower_tails':constants_low,'constant_positive_upper_tails':constants_high}
