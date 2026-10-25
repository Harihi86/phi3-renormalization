"""Exact measure-preserving averaging of the child fraction reflections."""
from itertools import product
import numpy as np
import primitive_child_triangle_host_vector_v2 as V
R=V.R


def evaluate(points,binding,positive,proposal,control=True,radial_transport=True,simplex_orbit=6,fraction_orbit=8):
    assert fraction_orbit in (1,2,8)
    flips=[(0,0,0)] if fraction_orbit==1 else [(0,0,0),(1,1,1)] if fraction_orbit==2 else list(product((0,1),repeat=3))
    p=np.asarray(points,dtype=np.longdouble);result=np.zeros(len(p),dtype=np.longdouble)
    for pattern in flips:
        x=p.copy()
        for j,flip in enumerate(pattern):
            if flip:x[:,6+j]=1-x[:,6+j]
        result+=V.evaluate(x,binding,positive,proposal,control,radial_transport,simplex_orbit)
    return result/len(flips)
