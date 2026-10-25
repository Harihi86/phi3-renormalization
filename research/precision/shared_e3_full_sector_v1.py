"""Exact complete Hepp-order sums quotiented only by proved forest-polynomial invariance."""
from fractions import Fraction as Q
from itertools import permutations
from collections import Counter
import math


def mask_permutation(mask,p):return sum(1<<p[e] for e in range(len(p)) if mask&(1<<e))
def polynomial_signature(terms,p):return tuple(sorted((mask_permutation(int(mask),p),Q(coefficient)) for mask,coefficient in terms))
def forest_signature(binding,p):
    return Counter((state['sign'],tuple(sorted(p[e] for e in state['weighted_edges'])),
        polynomial_signature(state['U'],p),polynomial_signature(state['F'],p)) for state in binding['paired_states'])


def invariant_group(binding,tables):
    assert binding['weighted_kind']=='shared';n=6;identity=tuple(range(n))
    reference=forest_signature(binding,identity);group=[]
    for p in permutations(range(n)):
        if forest_signature(binding,p)!=reference:continue
        if all(all(table[key][mask_permutation(mask,p)]==value for mask,value in table[key].items())
               for table in tables.values() for key in ('gamma','shapes','weights')):group.append(p)
    members=set(group);assert identity in members
    assert all(tuple(p[q[e]] for e in range(n)) in members for p in group for q in group)
    assert all(len(p)==len(set(p))==n for p in group)
    return group


def order_orbits(group,n=6):
    group=[tuple(p) for p in group];remaining=set(permutations(range(n)));result=[]
    while remaining:
        representative=min(remaining);orbit={tuple(p[e] for e in representative) for p in group}
        assert len(orbit)==len(group) and orbit<=remaining
        remaining-=orbit;result.append({'order':representative,'multiplicity':len(orbit)})
    assert sum(row['multiplicity'] for row in result)==math.factorial(n)
    return result


def validate_orbits(group,orbits,n=6):
    covered=set()
    for row in orbits:
        orbit={tuple(p[e] for e in row['order']) for p in group}
        assert len(orbit)==row['multiplicity']==len(group) and not covered&orbit
        covered|=orbit
    assert covered==set(permutations(range(n)))


def shape_map(shape,u):
    assert isinstance(u,Q) and 0<u<1
    if shape in (-2,-4):return u**(-shape),(-shape)*u**(-shape-1)
    assert type(shape) is int and shape>0
    denominator=1+(shape-1)*u
    return shape*u/denominator,Q(shape)/denominator**2


def map_complete(table,point,orbits):
    n=table['edges'];assert n==6 and len(point)==n-1 and all(isinstance(u,Q) and 0<u<1 for u in point)
    result=[]
    for row in orbits:
        order=row['order'];remaining=(1<<n)-1;alpha=[Q()]*n;alpha[order[0]]=Q(1)
        remaining^=1<<order[0];scale=Q(1);jacobian=Q(1);shapes=[]
        for j,u in enumerate(point,1):
            shape=table['shapes'][remaining];shapes.append(shape);ratio,derivative=shape_map(shape,u)
            scale*=ratio;alpha[order[j]]=scale;jacobian*=derivative*ratio**(n-j-1);remaining^=1<<order[j]
        assert remaining==0
        result.append({'order':order,'multiplicity':row['multiplicity'],'shapes':shapes,
            'alpha':alpha,'weight':row['multiplicity']*jacobian})
    return result


def hepp_moment(order,exponents):
    n=len(order);return math.prod(Q(1,n-j+sum(exponents[e] for e in order[j:])) for j in range(1,n))


def chart_moment(exponents):
    return sum((math.prod(Q(1,exponents[e]+1) for e in range(len(exponents)) if e!=pivot)
                for pivot in range(len(exponents))),Q())


def symmetrized_moment(exponents,group,orbits):
    total=Q()
    for row in orbits:
        # Invariant test polynomial is the average over the same edge group.
        terms=[]
        for p in group:
            changed=[0]*len(p)
            for e in range(len(p)):changed[p[e]]=exponents[e]
            terms.append(hepp_moment(row['order'],changed))
        total+=row['multiplicity']*sum(terms,Q())/len(group)
    assert total==chart_moment(exponents)
    return total
