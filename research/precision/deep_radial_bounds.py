"""Graph-selected deep nests and exact endpoint allowances for their 2D reduction."""
from fractions import Fraction as Q
from math import factorial
from exact_exp2_enclosures import log2

def parameters(tree):
    assert tree['external_legs'] in (2,3) and tree['loops']==4
    node=tree
    for loops in (4,3,2):
        assert node['loops']==loops and len(node['children'])==1
        node=node['children'][0];assert node['external_legs']==2
    assert node['loops']==1 and not node['children']
    return tree['external_legs']

def tail_polynomial(n,L,rate):
    return sum((Q(factorial(n),factorial(n-r)*rate**(r+1))*L**(n-r) for r in range(n+1)),Q())

def endpoint_bounds(external_legs):
    assert external_legs in (2,3)
    a=Q(2)**-10;to=16 if external_legs==2 else 48;ti=to+1
    To,Ti=Q(2)**to,Q(2)**ti;ln=log2()[1];Lo=(to+2)*ln;Li=(ti+1)*ln
    M1=Q(3 if external_legs==2 else 4);M2=Q(4) if external_legs==2 else 4+2*to*ln
    outer_high=(Q(1,3)+tail_polynomial(3,Lo,3))/To**3 if external_legs==2 else 2*(1+tail_polynomial(3,Lo,1))/To
    return {'log_outer_upper':to,'log_inner_upper':ti,'outer_low':a**5/20,'outer_high':outer_high,
        'inner_low':M1*a**5/120,'inner_high':M2/(12*Ti**3)*(Q(1,3)+tail_polynomial(2,Li,3)),
        'interpolation_domain':(a,Ti),'node_relative_gap_limit':Q(1,2000),'midpoint_relative_error':Q(1,4000)}

def interpolation_propagation(midpoint_integral):
    lo,hi=map(Q,midpoint_integral);assert 0<=lo<=hi
    error=Q(1,4000)
    return ((1-error)*lo,(1+error)*hi)
