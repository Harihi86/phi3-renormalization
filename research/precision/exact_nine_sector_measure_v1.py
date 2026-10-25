"""Independent integer enumeration of every nine-edge order and conditional partition."""
from fractions import Fraction as Q
from itertools import permutations
import math
import exact_rational_sector_transport_v1 as H


def verify(table):
    assert table['edges']==9
    full=511;shapes=table['shapes'];weights=table['weights']
    assert set(shapes)==set(range(1,full)) and set(weights)==set(range(full))
    assert all(type(s) is int and (s>=1 or s in (-2,-4)) for s in shapes.values())
    assert all(type(w) is int and w>0 for w in weights.values())
    common=math.lcm(*(s for s in shapes.values() if s>0))**8
    numerator=count=0
    for order in permutations(range(9)):
        mask=full;weight=common
        for edge in order[:-1]:
            mask^=1<<edge;shape=shapes[mask]
            if shape>0:
                assert weight%shape==0;weight//=shape
            else:weight*=-shape
        assert weight>0;numerator+=weight;count+=1
    normal=Q(numerator,common)
    assert count==362880 and normal==Q(table['normalization'],table['scale'])
    assert weights[0]==table['scale']
    for mask in range(1,full+1):
        choices=[e for e in range(9) if mask&(1<<e)]
        total=sum(weights[mask^(1<<e)] for e in choices)
        assert 0<total<2**62
        if mask==full:assert total==table['normalization']
        else:
            shape=shapes[mask];beta=Q(shape) if shape>0 else Q(1,-shape)
            assert weights[mask]*beta==total
        assert sum((Q(weights[mask^(1<<e)],total) for e in choices),Q())==1
    return {'all_order_count':count,'normalization':normal,'conditional_partitions':511,
            'integer_enumeration_common_denominator':common,'integer_enumeration_numerator':numerator}
