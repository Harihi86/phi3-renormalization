"""Complete order sums for two ordered, distinct bubble-weighted edges."""
from fractions import Fraction as Q
from itertools import permutations
from collections import Counter
import shared_e3_full_sector_v1 as M


def forest_signature(binding,p):
    # The remaining spectral representation distinguishes d from e pointwise.
    # Sorting this pair would allow an unproved interchange of the two roles.
    return Counter((state['sign'],tuple(p[e] for e in state['weighted_edges']),
        M.polynomial_signature(state['U'],p),M.polynomial_signature(state['F'],p))
        for state in binding['paired_states'])


def invariant_group(binding,tables):
    assert binding['weighted_kind']=='double'
    identity=tuple(range(6));reference=forest_signature(binding,identity);group=[]
    for p in permutations(range(6)):
        if forest_signature(binding,p)!=reference:continue
        if all(all(table[key][M.mask_permutation(mask,p)]==value for mask,value in table[key].items())
               for table in tables.values() for key in ('gamma','shapes','weights')):group.append(p)
    members=set(group);assert identity in members
    assert all(tuple(p[q[e]] for e in range(6)) in members for p in group for q in group)
    assert all(tuple(p[e] for e in binding['weighted_edges'])==tuple(binding['weighted_edges']) for p in group)
    return group


def polynomial(terms,alpha):
    total=Q()
    for mask,coefficient in terms:
        term=Q(coefficient)
        for e,a in enumerate(alpha):
            if int(mask)&(1<<e):term*=a
        total+=term
    return total


def coefficients(point,binding,table,orbits):
    assert len(point)==6 and all(type(x) is Q and 0<x<1 for x in point)
    assert binding['weighted_kind']=='double';weights=binding['weighted_edges'];assert len(weights)==2 and weights[0]!=weights[1]
    v=point[-1];r=1-v*v;tau=(1-r*r)/4;spectral=v*r**4/6
    a_coeff={};b_coeff={}
    def add(target,t,c):
        assert type(t) is Q and t>0 and type(c) is Q
        target[t]=target.get(t,Q())+c
    for row in M.map_complete(table,point[:5],orbits):
        alpha=row['alpha'];mass=sum(alpha,Q());d,e=[alpha[j] for j in weights]
        t0=(mass-d)/d;s0=(mass+e*(1/tau-1)-d)/d
        c=(tau-2)/(1-tau)**2;h=1/(1-tau);k=1/(tau*(1-tau)**2)
        for state in binding['paired_states']:
            assert state['weighted_edges']==weights
            U=polynomial(state['U'],alpha);F=polynomial(state['F'],alpha);assert U>0 and F>0
            shift=F/(U*d);t1=t0+shift;s1=s0+shift
            scale=Q(state['sign'])*row['weight']*spectral/U**3
            add(a_coeff,t0,scale*c);add(a_coeff,t1,-scale*c)
            add(a_coeff,s0,scale*k);add(a_coeff,s1,-scale*k)
            add(b_coeff,t0,-scale*h*e/(d*(1+t0)**2))
            add(b_coeff,t1,scale*h*e/(d*(1+t1)**2))
    a_coeff={t:c for t,c in a_coeff.items() if c};b_coeff={t:c for t,c in b_coeff.items() if c}
    assert sum(a_coeff.values(),Q())==0
    return a_coeff,b_coeff


def coefficient_rows(pair):
    a,b=pair
    return [{'t':t,'A':a.get(t,Q()),'B':b.get(t,Q())} for t in sorted(set(a)|set(b))]
