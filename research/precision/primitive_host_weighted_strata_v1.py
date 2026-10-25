"""Exact conditional Hepp probabilities for each position of the weighted edge."""
from fractions import Fraction as Q
from functools import lru_cache


def build(table,edge=8):
    H={int(k):v for k,v in table['H'].items()};gamma={int(k):v for k,v in table['gamma'].items()}
    bit=1<<edge
    @lru_cache(None)
    def weight(S,k):
        if not S:
            assert k==-1;return Q(1)
        if not S&bit:
            assert k==-1;return H[S]*gamma[S]
        assert 0<=k<S.bit_count()
        options=[e for e in range(9) if S&(1<<e) and ((e==edge)==(k==0))]
        return sum((weight(S^(1<<e),-1 if e==edge else k-1)/gamma.get(S^(1<<e),1) for e in options),Q())
    masses=[weight(511,k) for k in range(9)]
    assert sum(masses,Q())==table['normalization'] and all(m>0 for m in masses)
    # Materialize all valid conditional states for transparent independent replay.
    values={}
    for S in range(512):
        for k in (range(S.bit_count()) if S&bit else (-1,)):
            values[f'{S}:{k}']=weight(S,k)
    return {'edge':edge,'masses':masses,'normalization':table['normalization'],'weights':values}
