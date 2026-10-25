"""Triangle-root W3 interpolation with a logarithmic fourth-derivative bound."""
from fractions import Fraction as Q
from bisect import bisect_right
from convex_radial_envelope import Node
from cubic_spectral_envelope_v2 import bernstein,dyadic_upper
import fourier_root_envelope as E2

def derivative_bound(x):
    assert x>0;k=0
    while x*(1<<k)<1:k+=1
    # 4 log(2)<3: the first four positive terms of exp(3/4) already exceed2.
    return Q(10+3*k),k

def cubic(nodes,i,left,right):
    row=E2.cubic(nodes,i,left,right);bound,k=derivative_bound(min(nodes[j].x for j in row['indices']))
    row['fourth_derivative_remainder']*=50*bound
    row['W_error']=dyadic_upper(row['node_uncertainty']+row['fourth_derivative_remainder'])
    row['fourth_derivative_bound']=bound;row['dyadic_scale_index']=k
    return row

def inspect(nodes,cdata):
    assert len(nodes)>=4 and all(a.x<b.x for a,b in zip(nodes,nodes[1:]))
    cp=cdata['pieces'];xs=[n.x for n in nodes];cs=[p['left'] for p in cp];assert set(cs+[cp[-1]['right']])<=set(xs)
    assert xs[0]==cp[0]['left'] and xs[-1]==cp[-1]['right'];pieces=[];by_interval={};ce=we=Q()
    for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
        left,right=a.x,b.x;j=bisect_right(cs,left)-1;c=cp[j];assert c['left']<=left<right<=c['right']
        row=cubic(nodes,i,left,right);mc=tuple(c['midpoint']);ec=c['C_error'];fm=max(map(abs,bernstein(mc,left,right)))+max(left**-4,left**-2)/18
        if left>=1:
            n=(2*left).numerator//(2*left).denominator;fm=min(fm,3/left*Q(3,8)**n+ec)
        werr=dyadic_upper((right-left)*right**5*fm**2*row['W_error'])
        cerr=dyadic_upper((right-left)*right**5*b.upper*(2*fm*ec+ec**2))
        row.update(left=left,right=right,C_midpoint=mc,C_piece=j,W_node_interval=i,F_midpoint_bound=fm,W_right_upper=b.upper,W_propagated_error=werr,C_propagated_error=cerr)
        pieces.append(row);by_interval[i]=werr;we+=werr;ce+=cerr
    return {'pieces':pieces,'by_interval':by_interval,'W_propagated_error':we,'C_propagated_error':ce,'total_propagated_error':ce+we}
