"""Exact Hepp weights for a primitive homogeneous rank-three host."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations,permutations
from collections import Counter
import math


def build(binding,enumerate_orders=True):
    U=binding['host_U'];assert len(binding['host_routes'])==9 and all(co==1 for _,co in U)
    masks=tuple(mask for mask,_ in U);full=(1<<9)-1
    nu={s:min((s&m).bit_count() for m in masks) for s in range(full+1)}
    gamma={s:s.bit_count()-3*nu[s] for s in range(1,full)}
    assert all(g>=1 for g in gamma.values()),'Host has another ultraviolet subdivergence'
    assert nu[0]==0 and nu[full]==3 and all(m.bit_count()==3 for m in masks)
    # Linear-matroid basis exchange guarantees a common leading monomial on
    # every nested chain, not merely unrelated minima for individual subsets.
    exchange=0
    for a in masks:
        for b in masks:
            for i in range(9):
                if a&(1<<i) and not b&(1<<i):
                    assert any(((a^(1<<i))|(1<<j)) in masks for j in range(9) if b&(1<<j) and not a&(1<<j))
                    exchange+=1
    H={0:Q(1)}
    for size in range(1,9):
        for s in range(1,full):
            if s.bit_count()==size:H[s]=sum((H[s^(1<<i)] for i in range(9) if s&(1<<i)),Q())/gamma[s]
    total=sum((H[full^(1<<i)] for i in range(9)),Q())
    histogram=Counter()
    if enumerate_orders:
        for order in permutations(range(9)):
            remaining=full;denominator=1
            for e in order[:-1]:remaining^=1<<e;denominator*=gamma[remaining]
            histogram[denominator]+=1
        assert sum(histogram.values())==math.factorial(9)
        assert sum((Q(count,den) for den,count in histogram.items()),Q())==total
    cap=Q(binding['external_energy_trial_bound']);absolute=total*cap*(9*(1+cap)+8)/6
    return {'gamma':gamma,'nu':nu,'H':H,'normalization':total,'basis_exchange_checks':exchange,
            'sector_orders':math.factorial(9),'denominator_histogram':dict(histogram),
            'full_absolute_bound':absolute,'whole_integral_evaluated':False}


def decode(cube,table):
    assert len(cube)==9 and all(type(x) in (int,Q) and 0<x<1 for x in cube)
    remaining=511;u=Q(cube[0]);order=[];exponents=[]
    H=table['H'];gamma=table['gamma']
    while remaining:
        options=[i for i in range(9) if remaining&(1<<i)]
        denom=sum((H[remaining^(1<<i)] for i in options),Q());left=Q()
        for e in options:
            probability=H[remaining^(1<<e)]/denom
            if u<left+probability:
                u=(u-left)/probability;order.append(e);remaining^=1<<e
                if remaining:exponents.append(gamma[remaining])
                break
            left+=probability
        else:raise AssertionError('Hepp interval partition has a gap')
    assert sorted(order)==list(range(9)) and len(exponents)==8 and 0<=u<1
    return {'order':order,'gamma':exponents,'ratios':list(map(Q,cube[1:])),'auxiliary':u}
