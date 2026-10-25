"""An exact positive rectangular proposal: one selector and seven local coordinates."""
from fractions import Fraction as Q
from itertools import combinations
import math
UNITS=2**24
DIM=7

def leaves(model):
    assert type(model) is dict and set(model)=={'tree'}
    result=[]
    def walk(node,low,high,depth):
        assert type(node) is dict and depth<=12
        if set(node)=={'mass'}:
            assert type(node['mass']) is int and node['mass']>0
            result.append({'low':low,'high':high,'mass':node['mass']});return
        assert set(node)=={'axis','left','right'} and type(node['axis']) is int and 0<=node['axis']<DIM
        axis=node['axis'];middle=(low[axis]+high[axis])/2
        lhi=high.copy();lhi[axis]=middle;rlo=low.copy();rlo[axis]=middle
        walk(node['left'],low.copy(),lhi,depth+1);walk(node['right'],rlo,high.copy(),depth+1)
    walk(model['tree'],[Q()]*DIM,[Q(1)]*DIM,0)
    assert 1<=len(result)<=64 and sum(row['mass'] for row in result)==UNITS
    assert sum((math.prod(b-a for a,b in zip(row['low'],row['high'])) for row in result),Q())==1
    return result

def transform(model,point):
    assert len(point)==8 and all(type(x) in (int,Q) and 0<x<1 for x in point)
    cumulative=0
    for row in leaves(model):
        cumulative+=row['mass']
        if point[0]<Q(cumulative,UNITS):
            widths=[b-a for a,b in zip(row['low'],row['high'])]
            mapped=[a+w*x for a,w,x in zip(row['low'],widths,point[1:])]
            weight=math.prod(widths)*Q(UNITS,row['mass'])
            assert len(mapped)==7 and all(0<x<1 for x in mapped) and weight>0
            return mapped,weight
    raise AssertionError('positive complete selector partition required')

def exact_checks(model):
    rows=leaves(model);powers=set()
    for axis in range(DIM):
        for degree in range(5):p=[0]*DIM;p[axis]=degree;powers.add(tuple(p))
    for a,b in combinations(range(DIM),2):p=[0]*DIM;p[a]=p[b]=1;powers.add(tuple(p))
    powers.add((1,)*DIM)
    for power in sorted(powers):
        integral=Q()
        for row in rows:
            volume=math.prod(b-a for a,b in zip(row['low'],row['high']))
            probability=Q(row['mass'],UNITS);weight=volume/probability
            local=math.prod((b**(k+1)-a**(k+1))/(Q(k+1)*(b-a)) for a,b,k in zip(row['low'],row['high'],power))
            integral+=probability*weight*local
        assert integral==math.prod(Q(1,k+1) for k in power)
    cumulative=0;boundaries=0
    for row in rows[:-1]:
        cumulative+=row['mass'];boundary=Q(cumulative,UNITS)
        for delta in (-Q(1,2**40),Q(),Q(1,2**40)):
            point=[boundary+delta]+[Q(2,5)]*7;mapped,weight=transform(model,point);assert all(0<x<1 for x in mapped) and weight>0;boundaries+=1
    return {'positive_leaves':len(rows),'exact_polynomial_moments':len(powers),'selector_boundary_points':boundaries,
            'full_domain_volume':Q(1),'all_positive_probabilities':True,'dimension':8,'physical_coordinates':7}
