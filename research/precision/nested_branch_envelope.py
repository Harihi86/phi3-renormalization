"""Three positive profiles with signed recombination and constant uncertainty."""
from fractions import Fraction as Q
from dataclasses import asdict
from convex_radial_envelope import Node
import cubic_spectral_envelope_v2 as CM
from fourier_endpoint_bounds import up

def inspect(first,continuum,mass,constants):
    profiles=[]
    for nodes in (first,continuum,mass):
        value=CM.inspect(nodes);value['nodes']=[asdict(n) for n in nodes];profiles.append(value)
    xs=[n.x for n in first];assert all([n.x for n in nodes]==xs for nodes in (continuum,mass))
    intervals=[constants[key] for key in ('aN','bN')];assert all(l<=h for l,h in intervals)
    midpoint=[(l+h)/2 for l,h in intervals];error=[(h-l)/2 for l,h in intervals]
    assert max(map(abs,intervals[0]))<=Q(17,1620) and max(map(abs,intervals[1]))<=Q(11,1620)
    rows=[];total=Q();by_interval={}
    for i,(p1,pc,pz) in enumerate(zip(*(p['pieces'] for p in profiles))):
        l,h=xs[i:i+2];e1=p1['C_error'];ec=pc['C_error'];ez=pz['C_error'];Y=l**-4;Z=l**-2
        if l>=1:
            power=Q(3,8)**((2*l).numerator//(2*l).denominator);Y=min(Y,12*power/l**2);Z=min(Z,12*power/l)
        eN=ec+(e1+ez)/18+error[0]*Y+error[1]*Z
        m1=max(map(abs,CM.bernstein(p1['midpoint'],l,h)));mc=max(map(abs,CM.bernstein(pc['midpoint'],l,h)));mz=max(map(abs,CM.bernstein(pz['midpoint'],l,h)))
        f1=m1+max(l**-4,l**-2)/18;fN=mc+(m1+mz)/18+abs(midpoint[0])*Y+abs(midpoint[1])*Z
        if l>=1:f1=min(f1,3*power/l+e1);fN=min(fN,power+eN)
        product=up((h-l)*h**5*min(h**4/24,h*h/3)*(f1*eN+fN*e1+e1*eN));total+=product;by_interval[i]=product
        rows.append({'node_interval':i,'left':l,'right':h,'C1_midpoint':p1['midpoint'],'Ccont_midpoint':pc['midpoint'],'Cz_midpoint':pz['midpoint'],'C1_error':e1,'Ccont_error':ec,'Cz_error':ez,'aN_bN_midpoint':midpoint,'aN_bN_error':error,'Y1_bound':Y,'Z1_bound':Z,'FN_error':eN,'F1_midpoint_bound':f1,'FN_midpoint_bound':fN,'mixed_error_term':e1*eN,'propagated_error':product})
    return {'profiles':profiles,'pieces':rows,'aN_bN_intervals':intervals,'aN_bN_midpoint':midpoint,'aN_bN_error':error,'propagated_error':total,'by_interval':by_interval,'profile_selection_scores_are_not_physical_error_claims':True,'signed_product_not_square':True,'shared_C1_errors_not_assumed_independent':True,'constant_intervals_propagated':True}
