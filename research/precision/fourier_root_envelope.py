"""Absolute fourth-derivative interpolation for W; original C enclosure retained."""
from fractions import Fraction as Q
from bisect import bisect_right
from convex_radial_envelope import Node
from cubic_spectral_envelope_v2 import add,scale,mul,bernstein,value,dyadic_upper

M2_BOUND=Q(1,50)

def cubic(nodes,i,left,right):
    assert len(nodes)>=4 and 0<=i<len(nodes)-1 and nodes[i].x<=left<right<=nodes[i+1].x
    start=max(0,min(i-1,len(nodes)-4));indices=list(range(start,start+4));pol=(Q(),);uncertainty=Q();product=(Q(1),)
    for j in indices:
        basis=(Q(1),)
        for k in indices:
            if k!=j:basis=scale(mul(basis,(-nodes[k].x,Q(1))),1/(nodes[j].x-nodes[k].x))
        midpoint=(nodes[j].lower+nodes[j].upper)/2;radius=(nodes[j].upper-nodes[j].lower)/2
        pol=add(pol,scale(basis,midpoint));uncertainty+=radius*max(map(abs,bernstein(basis,left,right)))
        product=mul(product,(-nodes[j].x,Q(1)))
    remainder=M2_BOUND*max(map(abs,bernstein(product,left,right,4)))/24
    return {'midpoint':pol,'indices':indices,'node_uncertainty':uncertainty,'fourth_derivative_remainder':remainder,'W_error':dyadic_upper(uncertainty+remainder)}

def inspect(nodes,cdata):
    assert len(nodes)>=4 and all(a.x<b.x for a,b in zip(nodes,nodes[1:]))
    cp=cdata['pieces'];assert cp and nodes[0].x==cp[0]['left'] and nodes[-1].x==cp[-1]['right']
    assert all(a['right']==b['left'] for a,b in zip(cp,cp[1:]));starts=[p['left'] for p in cp]
    pieces=[];by_interval={}
    for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
        cuts=sorted({a.x,b.x,*[x for x in starts if a.x<x<b.x]});error=Q()
        for left,right in zip(cuts,cuts[1:]):
            j=min(bisect_right(starts,left)-1,len(cp)-1);c=cp[j];assert c['left']<=left<right<=c['right']
            row=cubic(nodes,i,left,right);mc=tuple(c['midpoint']);ec=c['C_error']
            fmax=max(map(abs,bernstein(mc,left,right)))+max(left**-4,left**-2)/18
            if left>=1:
                floor=(2*left).numerator//(2*left).denominator
                fmax=min(fmax,3/left*Q(3,8)**floor+ec)
            propagated=dyadic_upper((right-left)*right**5*fmax*fmax*row['W_error']);error+=propagated
            row.update(left=left,right=right,C_midpoint=mc,C_piece=j,W_node_interval=i,F_midpoint_bound=fmax,W_propagated_error=propagated)
            pieces.append(row)
        by_interval[i]=error
    werror=sum(by_interval.values(),Q());cerror=dyadic_upper(cdata['propagated_error']*M2_BOUND)
    return {'pieces':pieces,'by_interval':by_interval,'W_propagated_error':werror,'C_propagated_error':cerror,'total_propagated_error':werror+cerror}
