"""Exact full-domain squared endpoint maps and a reflected base y coordinate."""
from fractions import Fraction as Q


def expand(points):
    mapped=[];weights=[]
    for point in points:
        assert len(point)==5 and all(type(v) in (int,Q) and 0<v<1 for v in point)
        r,v,x,y,u=map(Q,point)
        mapped.extend(((r*r,v*v,x,1-(1-y)**2,u),(r*r,v*v,x,1-y*y,u)))
        weights.extend((4*r*v*(1-y),4*r*v*y))
    assert all(all(0<v<1 for v in p) for p in mapped) and all(w>0 for w in weights)
    return mapped,weights


def collapse(values,weights):
    assert len(values)==len(weights) and len(values)>0 and len(values)%2==0
    assert all(type(w) is Q and w>0 for w in weights)
    assert all(type(lo) is Q and type(hi) is Q and lo<=hi for lo,hi in values)
    result=[]
    for i in range(0,len(values),2):
        result.append(tuple(sum((weights[j]*values[j][k] for j in (i,i+1)),Q()) for k in (0,1)))
    return result
