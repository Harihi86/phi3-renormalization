"""Exact complete-sector proposal with rational maps matched to support scores."""
from fractions import Fraction as Q
from itertools import permutations
import math
FULL=255


def build(gamma,critical_power=2):
    assert critical_power in (2,4) and set(gamma)==set(range(1,FULL))
    assert all(type(v) is int for v in gamma.values())
    shapes={s:(gamma[s] if gamma[s]>=1 else -critical_power) for s in gamma}
    beta={s:(Q(v) if v>0 else Q(1,-v)) for s,v in shapes.items()}
    unscaled={0:Q(1)}
    for size in range(1,8):
        for mask in range(1,FULL):
            if mask.bit_count()==size:
                unscaled[mask]=sum((unscaled[mask^(1<<e)] for e in range(8) if mask&(1<<e)),Q())/beta[mask]
    scale=math.lcm(*(v.denominator for v in unscaled.values()))
    weights={s:int(v*scale) for s,v in unscaled.items()}
    normal=sum(weights[FULL^(1<<e)] for e in range(8))
    assert max([normal,scale,*weights.values()])<2**62 and min(weights.values())>0
    return {'critical_power':critical_power,'gamma':gamma,'shapes':shapes,'weights':weights,
            'normalization':normal,'scale':scale,'exact_normalization':Q(normal,scale)}


def decode(table,selector):
    assert type(selector) in (Q,int) and 0<=selector<1
    remaining=FULL;order=[];shapes=[];probability=Q(1);u=Q(selector)
    while remaining:
        choices=[e for e in range(8) if remaining&(1<<e)]
        total=sum(table['weights'][remaining^(1<<e)] for e in choices);left=Q()
        for edge in choices:
            p=Q(table['weights'][remaining^(1<<edge)],total)
            if u<left+p:
                u=(u-left)/p;probability*=p;order.append(edge);remaining^=1<<edge
                if remaining:shapes.append(table['shapes'][remaining])
                break
            left+=p
        else:raise AssertionError('exact conditional partition has a gap')
    inverse_betas=math.prod(Q(1,v) if v>0 else Q(-v) for v in shapes)
    assert probability==inverse_betas*Q(table['scale'],table['normalization'])
    return order,shapes,probability


def ratio(shape,u):
    assert type(shape) is int and (shape>=1 or shape in (-2,-4)) and 0<u<1
    if shape<0:
        power=-shape;return u**power,power*u**(power-1)
    denominator=1+(shape-1)*u
    return shape*u/denominator,Q(shape)/denominator**2


def map_orbit(table,point):
    assert len(point)==8 and all(type(u) in (Q,int) and 0<u<1 for u in point)
    rows=[]
    for orbit in range(8):
        order,shapes,probability=decode(table,(point[0]+Q(orbit,8))%1)
        alpha=[Q()]*8;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,shape) in enumerate(zip(point[1:],shapes),1):
            r,derivative=ratio(shape,u);scale*=r;alpha[order[j]]=scale
            jacobian*=derivative*r**(7-j)
        rows.append({'order':order,'shapes':shapes,'probability':probability,'alpha':alpha,'weight':jacobian/probability})
    return rows


def exact_checks(table):
    normal=Q();count=0
    for order in permutations(range(8)):
        remaining=FULL;mass=Q(1)
        for edge in order[:-1]:
            remaining^=1<<edge;shape=table['shapes'][remaining]
            mass*=Q(1,shape) if shape>0 else Q(-shape)
        normal+=mass;count+=1
    assert count==math.factorial(8) and normal==Q(table['normalization'],table['scale'])
    assert table['weights'][0]==table['scale'] and min(table['weights'].values())>0
    for mask in range(1,FULL+1):
        options=[e for e in range(8) if mask&(1<<e)]
        total=sum(table['weights'][mask^(1<<e)] for e in options)
        if mask!=FULL:
            shape=table['shapes'][mask];beta=Q(shape) if shape>0 else Q(1,-shape)
            assert Q(table['weights'][mask])*beta==total
        else:assert total==table['normalization']
        assert sum((Q(table['weights'][mask^(1<<e)],total) for e in options),Q())==1
    return {'complete_orders':count,'conditional_partitions':255,'normalization':normal}
