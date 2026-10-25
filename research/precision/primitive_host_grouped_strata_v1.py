"""Exact first-group and weighted-group-rank conditional sector masses."""
from fractions import Fraction as Q
from functools import lru_cache


def build(table,kind):
    assert kind in ('first','rank','joint')
    edge=next(i for i,group in enumerate(table['groups']) if 8 in group);bit=1<<edge
    H={int(k):v for k,v in table['H'].items()};gamma={int(k):v for k,v in table['gamma'].items()};Z=table['normalization']
    @lru_cache(None)
    def weight(S,k):
        if not S:assert k==-1;return Q(1)
        if k==-1:return Z if S==63 else H[S]*gamma[S]
        assert S&bit and 0<=k<S.bit_count()
        options=[e for e in range(6) if S&(1<<e) and ((e==edge)==(k==0))]
        return sum((weight(S^(1<<e),-1 if e==edge else k-1)/gamma.get(S^(1<<e),1) for e in options),Q())
    strata=[]
    if kind=='first':
        strata=[{'first':e,'position':-1,'probability':H[63^(1<<e)]/Z} for e in range(6)]
    elif kind=='rank':strata=[{'first':-1,'position':k,'probability':weight(63,k)/Z} for k in range(6)]
    else:
        for e in range(6):
            for k in range(6):
                if (e==edge)!=(k==0):continue
                child=63^(1<<e);q=weight(child,-1 if e==edge else k-1)/gamma[child]/Z
                strata.append({'first':e,'position':k,'probability':q})
    assert sum((s['probability'] for s in strata),Q())==1 and all(s['probability']>0 for s in strata)
    values={}
    for S in range(64):
        for k in [-1]+(list(range(S.bit_count())) if S&bit else []):values[f'{S}:{k}']=weight(S,k)
    return {'kind':kind,'edge':edge,'strata':strata,'weights':values,'normalization':Z}
