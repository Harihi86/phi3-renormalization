"""Exact affine envelopes for the positive decreasing convex spectral C(r)."""
from bisect import bisect_right
from fractions import Fraction as Q
from convex_radial_envelope import Node,Line,Piece
class SpectralEnvelope:
    def __init__(self,nodes):
        self.nodes=tuple(nodes);assert len(nodes)>=2
        self.coordinates=tuple(n.x for n in nodes)
        assert self.coordinates[0]>0 and all(a<b for a,b in zip(self.coordinates,self.coordinates[1:]))
        minimum_upper=nodes[0].upper
        for node in nodes:
            assert node.lower<=minimum_upper;minimum_upper=min(minimum_upper,node.upper)
    def pieces(self):
        out=[]
        for i,(a,b) in enumerate(zip(self.nodes,self.nodes[1:])):
            width=b.x-a.x;slope=(b.upper-a.upper)/width;upper=Line(slope,a.upper-slope*a.x)
            lines=[Line(Q(),Q()),Line(Q(),b.lower)]
            if i:
                prev=self.nodes[i-1];slope=(a.lower-prev.upper)/(a.x-prev.x);lines.append(Line(slope,a.lower-slope*a.x))
            if i+2<len(self.nodes):
                after=self.nodes[i+2];slope=(after.upper-b.lower)/(after.x-b.x);lines.append(Line(slope,b.lower-slope*b.x))
            cuts={a.x,b.x}
            for j,first in enumerate(lines):
                for second in lines[j+1:]:
                    if first.slope!=second.slope:
                        x=(second.intercept-first.intercept)/(first.slope-second.slope)
                        if a.x<x<b.x:cuts.add(x)
            ordered=sorted(cuts)
            for left,right in zip(ordered,ordered[1:]):
                lower=max(lines,key=lambda line:line((left+right)/2))
                assert 0<=lower(left)<=upper(left) and 0<=lower(right)<=upper(right)
                out.append(Piece(left,right,lower,upper))
        return tuple(out)
def inspect(nodes):
    env=SpectralEnvelope(nodes);pieces=[];by_interval={}
    for p in env.pieces():
        gap=p.gap_upper();assert gap>=0;err=gap/2
        # C, Y and Z are nonnegative. No sign assumption on the difference F.
        Cmax=max(p.upper(p.left),p.upper(p.right));pole=max(p.left**-4,p.left**-2)/18
        fmax=Cmax+pole;weight=p.right**5*min(p.right**4/24,p.right**2/3)
        propagated=(p.right-p.left)*weight*(2*fmax*err+err*err)
        i=bisect_right(env.coordinates,(p.left+p.right)/2)-1
        by_interval[i]=by_interval.get(i,Q())+propagated
        pieces.append({'left':p.left,'right':p.right,'lower_slope':p.lower.slope,'lower_intercept':p.lower.intercept,'upper_slope':p.upper.slope,'upper_intercept':p.upper.intercept,
                       'mid_slope':(p.lower.slope+p.upper.slope)/2,'mid_intercept':(p.lower.intercept+p.upper.intercept)/2,
                       'absolute_C_error':err,'density_error_integral_bound':propagated,'node_interval':i})
    assert pieces[0]['left']==nodes[0].x and pieces[-1]['right']==nodes[-1].x
    assert all(a['right']==b['left'] for a,b in zip(pieces,pieces[1:]))
    return {'pieces':pieces,'by_interval':by_interval,'propagated_error':sum(by_interval.values(),Q())}
