"""Exact full-domain Hepp proposal; power estimates select a density, not a proof of convergence."""
from fractions import Fraction as Q
from itertools import permutations
import math
FULL=255

def build(binding,power=2):
    assert power in (2,4)
    gamma={};powers={};weights={0:1}
    for mask in range(1,FULL):
        estimates=[]
        for state in binding['paired_states']:
            nu=min((mask&int(m)).bit_count() for m,c in state['U'])
            nf=min((mask&int(m)).bit_count() for m,c in state['F'])
            assert nf>=nu
            estimates.append(mask.bit_count()+2*nf-5*nu)
        gamma[mask]=min(estimates);powers[mask]=power if gamma[mask]<=1 else 1
    for size in range(1,8):
        for mask in range(1,FULL):
            if mask.bit_count()==size:weights[mask]=powers[mask]*sum(weights[mask^(1<<edge)] for edge in range(8) if mask&(1<<edge))
    total=sum(weights[FULL^(1<<edge)] for edge in range(8))
    assert 0<total<=math.factorial(8)*power**7
    return {'power':power,'gamma':gamma,'powers':powers,'weights':weights,'normalization':total}

def decode(table,selector):
    assert type(selector) in (Q,int) and 0<=selector<1
    remaining=FULL;order=[];exponents=[];probability=Q(1);u=Q(selector)
    while remaining:
        options=[edge for edge in range(8) if remaining&(1<<edge)]
        total=sum(table['weights'][remaining^(1<<edge)] for edge in options);left=Q()
        for edge in options:
            p=Q(table['weights'][remaining^(1<<edge)],total)
            if u<left+p:
                u=(u-left)/p;probability*=p;order.append(edge);remaining^=1<<edge
                if remaining:exponents.append(table['powers'][remaining])
                break
            left+=p
        else:raise AssertionError('exact Hepp conditional partition has a gap')
    assert probability==Q(math.prod(exponents),table['normalization']) and sorted(order)==list(range(8))
    return order,exponents,probability

def transform(table,point):
    assert len(point)==8 and all(type(q) in (int,Q) and 0<q<1 for q in point)
    order,powers,probability=decode(table,point[0]);alpha=[Q()]*8;alpha[order[0]]=Q(1);scale=Q(1);weight=Q(table['normalization'])
    for j,(u,power) in enumerate(zip(point[1:],powers),1):
        r=u**power;scale*=r;alpha[order[j]]=scale;weight*=r**(8-j)/u
    assert all(0<a<=1 for a in alpha) and weight>0
    return alpha,weight,{'order':order,'powers':powers,'probability':probability}

def exact_moment(exponents):
    assert len(exponents)==8 and all(type(x) is int and x>=0 for x in exponents)
    moment={0:Q(1)}
    for size in range(1,8):
        for mask in range(1,FULL):
            if mask.bit_count()==size:
                denominator=size+sum(exponents[j] for j in range(8) if mask&(1<<j))
                moment[mask]=sum((moment[mask^(1<<j)] for j in range(8) if mask&(1<<j)),Q())/denominator
    value=sum((moment[FULL^(1<<j)] for j in range(8)),Q())
    expected=sum((math.prod(Q(1,exponents[e]+1) for e in range(8) if e!=pivot) for pivot in range(8)),Q())
    assert value==expected;return value

def exact_checks(table):
    count=weight=0
    for order in permutations(range(8)):
        remaining=FULL;mass=1
        for edge in order[:-1]:remaining^=1<<edge;mass*=table['powers'][remaining]
        weight+=mass;count+=1
    assert count==math.factorial(8) and weight==table['normalization']
    partitions=0
    for mask in range(1,FULL+1):
        total=sum(table['weights'][mask^(1<<edge)] for edge in range(8) if mask&(1<<edge))
        probabilities=[Q(table['weights'][mask^(1<<edge)],total) for edge in range(8) if mask&(1<<edge)]
        assert min(probabilities)>0 and sum(probabilities,Q())==1;partitions+=1
    checks=0;cumulative=0
    for edge in range(7):
        cumulative+=table['weights'][FULL^(1<<edge)];boundary=Q(cumulative,table['normalization'])
        for delta in (-Q(1,2**50),Q(),Q(1,2**50)):
            decoded=decode(table,boundary+delta);assert decoded[0][0]==(edge if delta<0 else edge+1);checks+=1
    return {'complete_orders':count,'normalization':weight,'conditional_partitions':partitions,'first_pivot_boundary_checks':checks}
