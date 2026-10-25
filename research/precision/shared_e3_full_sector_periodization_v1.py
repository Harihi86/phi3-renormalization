"""Exact endpoint smoothing of five physical ratios after complete order summation.

This composes a bijective cube transport with an already proved integrand.
It does not assert finite variance or improved Sobol error for that integrand.
"""
from fractions import Fraction as Q
from math import prod


def axis(u,order):
    assert type(u) is Q and 0<u<1 and type(order) is int and order in (1,2)
    if order==1:return u*u*(3-2*u),6*u*(1-u)
    return u**3*(10-15*u+6*u*u),30*u*u*(1-u)**2


def transform(point,order):
    assert len(point)==5
    rows=[axis(u,order) for u in point]
    mapped=[u for u,j in rows];jacobian=prod(j for u,j in rows)
    assert all(0<u<1 for u in mapped) and jacobian>0
    return mapped,jacobian


def weighted_interval(interval,jacobian):
    assert len(interval)==2 and all(type(x) is Q for x in interval)
    assert interval[0]<=interval[1] and type(jacobian) is Q and jacobian>0
    return interval[0]*jacobian,interval[1]*jacobian
