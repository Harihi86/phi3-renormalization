"""Exact envelopes for a separately proved nonnegative increasing convex function.

This module does not certify numerical node values or establish the shape
assumption. A physical caller must bind both to its original recursive graph.
"""
from bisect import bisect_right
from dataclasses import dataclass
from fractions import Fraction as Q

@dataclass(frozen=True)
class Node:
    x: Q
    lower: Q
    upper: Q
    def __post_init__(self):
        for field in ('x','lower','upper'):object.__setattr__(self,field,Q(getattr(self,field)))
        if self.x<0 or self.lower>self.upper or self.upper<0:raise ValueError('Nonnegative argument and consistent function interval required')

@dataclass(frozen=True)
class Line:
    slope: Q
    intercept: Q
    def __call__(self,x):return self.slope*Q(x)+self.intercept

@dataclass(frozen=True)
class Piece:
    left: Q
    right: Q
    lower: Line
    upper: Line
    def gap_upper(self):return max(self.upper(self.left)-self.lower(self.left),self.upper(self.right)-self.lower(self.right))

class ConvexEnvelope:
    def __init__(self,nodes):
        self.nodes=tuple(nodes)
        if len(self.nodes)<2:raise ValueError('At least two certified nodes required')
        self.coordinates=tuple(n.x for n in self.nodes)
        if any(a>=b for a,b in zip(self.coordinates,self.coordinates[1:])):raise ValueError('Strictly increasing coordinates required')
        maximum_lower=Q(0)
        for node in self.nodes:
            maximum_lower=max(maximum_lower,node.lower)
            if maximum_lower>node.upper:raise ValueError('Node intervals contradict monotonicity')
    def _index(self,x):
        x=Q(x)
        if not self.coordinates[0]<=x<=self.coordinates[-1]:raise ValueError('No extrapolation outside certified nodes')
        return min(bisect_right(self.coordinates,x)-1,len(self.nodes)-2)
    def _lines(self,i):
        a,b=self.nodes[i:i+2];width=b.x-a.x
        slope=(b.upper-a.upper)/width;upper=Line(slope,a.upper-slope*a.x)
        lower=[Line(Q(0),Q(0)),Line(Q(0),a.lower)]
        if i:
            previous=self.nodes[i-1];slope=(a.lower-previous.upper)/(a.x-previous.x)
            lower.append(Line(slope,a.lower-slope*a.x))
        if i+2<len(self.nodes):
            following=self.nodes[i+2];slope=(following.upper-b.lower)/(following.x-b.x)
            lower.append(Line(slope,b.lower-slope*b.x))
        return lower,upper
    def bounds(self,x):
        x=Q(x);lines,upper=self._lines(self._index(x));lo=max(line(x) for line in lines);hi=upper(x)
        if lo>hi:raise ValueError('Node intervals contradict the convex shape assumption')
        return lo,hi
    def pieces(self):
        """Split at all lower-envelope kinks; native jets may not cross them."""
        result=[]
        for i in range(len(self.nodes)-1):
            a,b=self.coordinates[i:i+2];lines,upper=self._lines(i);cuts={a,b}
            for j,first in enumerate(lines):
                for second in lines[j+1:]:
                    if first.slope==second.slope:continue
                    x=(second.intercept-first.intercept)/(first.slope-second.slope)
                    if a<x<b:cuts.add(x)
            ordered=sorted(cuts)
            for left,right in zip(ordered,ordered[1:]):
                lower=max(lines,key=lambda line:line((left+right)/2))
                if lower(left)>upper(left) or lower(right)>upper(right):raise ValueError('Inconsistent convex interpolation intervals')
                assert lower(left)>=0 and lower(right)>=0
                result.append(Piece(left,right,lower,upper))
        return tuple(result)

def positive_weight_error_bound(pieces,weight_mass_upper):
    """Certified mass bounds must refer to the same disjoint pieces."""
    if len(pieces)!=len(weight_mass_upper):raise ValueError('A mass bound is required for every interpolation piece')
    masses=tuple(Q(v) for v in weight_mass_upper)
    if any(v<0 for v in masses):raise ValueError('Positive-kernel mass bounds must be nonnegative')
    return sum((piece.gap_upper()*mass for piece,mass in zip(pieces,masses)),Q())
