"""Outward rational endpoint bounds with the triangle logarithm retained."""
from fractions import Fraction as Q
from exact_exp2_enclosures import log2
from fourier_endpoint_bounds import tails as child_tails,up
def tails():
    prior=child_tails();eps=prior['left'];assert eps==Q(1,262144)
    L=18*log2()[1];d=Q(17,18)+Q(16,27)*L;e=7+4*L;b=Q(16,27)
    lower=up(eps**2*(e*(d*d+b*d+b*b/2)+2*d*d+4*b*d+3*b*b)/48)
    upper=up(prior['upper_endpoint_error']/8)
    return {'left':eps,'right':prior['right'],'lower_endpoint_error':lower,'upper_endpoint_error':upper,'total':lower+upper,'physical_amplitude_sign':-1}
