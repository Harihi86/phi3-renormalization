"""Independent full radial integrals at exact rational interior scale nodes."""
from fractions import Fraction as Q
import g521_mixture_reference_v1 as M


def evaluate(point,binding,control=True,digits=300,shape=10):
    assert shape in (10,11);count=shape-7
    nodes=[Q(2*i+1,2*count) for i in range(count)]
    weights=[Q(3,8),Q(1,4),Q(3,8)] if count==3 else [Q(13,48),Q(11,48),Q(11,48),Q(13,48)]
    def basis(value,i):
        result=Q(1)
        for j,node in enumerate(nodes):
            if j!=i:result*= (value-node)/(nodes[i]-node)
        return result
    ctx,result=M.evaluate(point,binding,control,digits,shape=6)
    for i,node0 in enumerate(nodes):
        for j,node1 in enumerate(nodes):
            pp=list(point);pp[8]=node0;pp[9]=node1
            coefficient=weights[i]*weights[j]-basis(point[8],i)*basis(point[9],j)
            _,value=M.evaluate(pp,binding,control,digits,shape=6)
            result+=ctx.mpf(coefficient.numerator)/coefficient.denominator*value
    return ctx,result
