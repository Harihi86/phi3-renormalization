"""Exact pointwise relative envelope gaps, with no sampled extremum."""
from fractions import Fraction as Q
from bisect import bisect_right

def inspect_envelope(envelope,left,right,tolerance):
    left,right,tolerance=map(Q,(left,right,tolerance))
    if not (envelope.coordinates[0]<=left<right<=envelope.coordinates[-1]) or tolerance<=0:raise ValueError('A contained positive-width domain and positive tolerance are required')
    failures=set();maximum=Q();records=[]
    for p in envelope.pieces():
        a=max(left,p.left);b=min(right,p.right)
        if b<=a:continue
        lows=[p.lower(x) for x in (a,b)]
        if min(lows)<=0:raise ValueError('Relative bounds require a strictly positive lower envelope')
        gaps=[(p.upper(x)-p.lower(x))/p.lower(x) for x in (a,b)]
        if min(gaps)<0:raise ValueError('Inconsistent envelope')
        # A ratio of affine functions with positive denominator is monotone or
        # constant. Its extrema over the entire closed piece are at endpoints.
        bound=max(gaps);maximum=max(maximum,bound)
        index=bisect_right(envelope.coordinates,(a+b)/2)-1
        if bound>tolerance:failures.add(index)
        records.append({'left':a,'right':b,'relative_gap_upper':bound,'node_interval':index})
    assert records and records[0]['left']==left and records[-1]['right']==right
    assert all(a['right']==b['left'] for a,b in zip(records,records[1:]))
    return {'maximum_relative_gap':maximum,'failing_intervals':sorted(failures),'pieces':records,'all_passed':not failures}

def dyadic_input(x):
    x=Q(x);n=x.numerator;d=x.denominator
    if n<=0 or d&(d-1):raise ValueError('Positive dyadic required')
    exponent=-(d.bit_length()-1)
    while n%2==0:n//=2;exponent+=1
    if not(1<=n<=2**24 and -30<=exponent<=60):raise ValueError('Outside native input domain')
    assert n*Q(2)**exponent==x
    return n,exponent
