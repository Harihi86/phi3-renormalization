"""Three CM profiles and propagation through the signed F1*F2 product."""
from fractions import Fraction as Q
import cubic_spectral_envelope_v2 as CM
Node=CM.Node

def profile(nodes):
    return {'nodes':[{'x':n.x,'lower':n.lower,'upper':n.upper} for n in nodes],**CM.inspect(nodes)}

def inspect(first,plus,minus):
    assert [n.x for n in first]==[n.x for n in plus]==[n.x for n in minus]
    profiles=[profile(nodes) for nodes in (first,plus,minus)];pieces=[];total=Q();by_interval={}
    for i,(a,p,m) in enumerate(zip(*(d['pieces'] for d in profiles))):
        left,right=a['left'],a['right'];assert (left,right)==(p['left'],p['right'])==(m['left'],m['right'])
        net=CM.add(p['midpoint'],CM.scale(m['midpoint'],-1));e1=a['C_error'];e2=p['C_error']+m['C_error']
        f1=max(map(abs,CM.bernstein(a['midpoint'],left,right)))+max(left**-4,left**-2)/18
        f2=max(map(abs,CM.bernstein(net,left,right)))+left**-4/180+left**-2/162+left**-1/648
        if left>=1:
            n=(2*left).numerator//(2*left).denominator;exponential=Q(3,8)**n
            f1=min(f1,3*exponential/left+e1);f2=min(f2,5*exponential+e2)
        # F1=M1+d1, F2=M2+d2. The mixed d1*d2 term is required.
        integrand_error=f1*e2+f2*e1+e1*e2
        propagated=CM.dyadic_upper((right-left)*right**5*min(right**4/24,right**2/3)*integrand_error)
        row={'left':left,'right':right,'node_interval':i,'C1_midpoint':a['midpoint'],'C2_net_midpoint':net,'C1_error':e1,'C2_error':e2,'F1_midpoint_bound':f1,'F2_midpoint_bound':f2,'mixed_error_term':e1*e2,'propagated_error':propagated}
        pieces.append(row);by_interval[i]=propagated;total+=propagated
    return {'profiles':profiles,'pieces':pieces,'propagated_error':total,'by_interval':by_interval,'profile_selection_scores_are_not_physical_error_claims':True,'signed_product_not_square':True}
