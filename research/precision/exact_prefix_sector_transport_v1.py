"""Sum initial edge choices exactly; retain the normalized conditional suffix sampler."""
from fractions import Fraction as Q
from itertools import permutations
import math
import exact_rational_sector_transport_v1 as A


def table_with_edges(table,n=8):
    result=dict(table);result['edges']=result.get('edges',n);return result


def strata(table,depth):
    table=table_with_edges(table);n=table['edges']
    assert type(depth) is int and 1<=depth<n
    full=(1<<n)-1;rows=[]
    for prefix in permutations(range(n),depth):
        remaining=full;left=Q();mass=Q(1)
        for edge in prefix:
            choices=[e for e in range(n) if remaining&(1<<e)]
            total=sum(table['weights'][remaining^(1<<e)] for e in choices)
            before=sum(table['weights'][remaining^(1<<e)] for e in choices if e<edge)
            left+=mass*Q(before,total)
            mass*=Q(table['weights'][remaining^(1<<edge)],total)
            remaining^=1<<edge
        rows.append({'prefix':list(prefix),'left':left,'mass':mass})
    assert len(rows)==math.factorial(n)//math.factorial(n-depth)
    return rows


def validate_strata(table,depth,rows):
    assert rows==strata(table,depth)
    edge=Q()
    for row in rows:
        assert row['left']==edge and row['mass']>0
        edge+=row['mass']
    assert edge==1
    for power in range(5):
        assert sum(((r['left']+r['mass'])**(power+1)-r['left']**(power+1) for r in rows),Q())/(power+1)==Q(1,power+1)
    return {'prefix_count':len(rows),'selector_moments':5,'exact_partition':True}


def mapped_rows(table,point,depth):
    table=table_with_edges(table);n=table['edges']
    assert len(point)==n and all(type(u) in (Q,int) and 0<u<1 for u in point)
    rows=[]
    for stratum in strata(table,depth):
        selector=stratum['left']+stratum['mass']*point[0]
        order,shapes,probability=A.decode(table,selector)
        assert order[:depth]==stratum['prefix']
        alpha=[Q()]*n;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,shape) in enumerate(zip(point[1:],shapes),1):
            ratio,derivative=A.ratio(shape,u);scale*=ratio;alpha[order[j]]=scale
            jacobian*=derivative*ratio**(n-j-1)
        rows.append({**stratum,'alpha':alpha,'order':order,'shapes':shapes,
            'probability':probability,'weight':jacobian*stratum['mass']/probability})
    return rows


def enumerate_order_masses(table):
    """Independent products over every full order, without the conditional decoder."""
    table=table_with_edges(table);n=table['edges'];full=(1<<n)-1
    totals={1:{},2:{}};count=0;normal=Q(table['scale'],table['normalization'])
    for order in permutations(range(n)):
        remaining=full;probability=normal
        for edge in order[:-1]:
            remaining^=1<<edge;shape=table['shapes'][remaining]
            probability*=Q(1,shape) if shape>0 else Q(-shape)
        assert probability>0
        for depth in totals:
            key=order[:depth];totals[depth][key]=totals[depth].get(key,Q())+probability
        count+=1
    assert count==math.factorial(n)
    for depth,masses in totals.items():
        assert sum(masses.values(),Q())==1
        assert masses=={tuple(r['prefix']):r['mass'] for r in strata(table,depth)}
    return {'orders':count,'prefix1':len(totals[1]),'prefix2':len(totals[2])}
