"""Graph-derived Hepp-sector importance distribution without factorial sectors.

For primitive E3 graphs the exact gap weights reproduce the proven Hepp
envelope. For subdivergent graphs positive clipped gaps define an importance
distribution only; this is not an analytic remainder certificate.
"""
from functools import lru_cache
from math import log
from .automatic_bounds import _cycle_rank


@lru_cache(maxsize=128)
def sector_distribution(g):
    size=len(g.edges);full=(1<<size)-1
    gaps=[0]*(full+1);values=[0.]*(full+1);values[0]=1.
    primitive=len(g.external)==3
    for mask in range(1,full):
        exact=mask.bit_count()-3*_cycle_rank(g,mask)
        if exact<=0:primitive=False
        gaps[mask]=max(1,exact)
        values[mask]=sum(values[mask^(1<<e)] for e in range(size) if mask&(1<<e))/gaps[mask]
    values[full]=sum(values[full^(1<<e)] for e in range(size))
    import numpy as np
    cumulative=np.zeros((full+1,size))
    for mask in range(1,full+1):
        weights=np.array([values[mask^(1<<e)] if mask&(1<<e) else 0. for e in range(size)])
        cumulative[mask]=np.cumsum(weights/weights.sum())
        last_active=mask.bit_length()-1
        cumulative[mask,last_active:]=1.
    return np.asarray(gaps),cumulative,values[full],primitive


def hepp_samples(g,uniforms):
    import numpy as np
    uniforms=np.asarray(uniforms)
    n=len(g.edges);count=len(uniforms)
    gaps,cdf,bound,primitive=sector_distribution(g)
    mask=np.full(count,(1<<n)-1,dtype=int);order=np.empty((count,n),dtype=int)
    for position in range(n-1,0,-1):
        choice=np.sum(uniforms[:,n-1-position,None]>cdf[mask],axis=1)
        order[:,position]=choice;mask^=1<<choice
    order[:,0]=np.rint(np.log2(mask)).astype(int)
    prefix=np.zeros(count,dtype=int);powers=np.empty((count,n-1))
    for r in range(1,n):
        prefix|=1<<order[:,r-1];powers[:,r-1]=gaps[prefix]
    logt=np.log(uniforms[:,n-1:])/powers
    ordered_log=np.zeros((count,n));ordered_log[:,:-1]=np.cumsum(logt[:,::-1],axis=1)[:,::-1]
    alpha=np.empty((count,n));alpha[np.arange(count)[:,None],order]=np.exp(ordered_log)
    total=alpha.sum(axis=1);alpha/=total[:,None]
    logweight=log(bound)-n*np.log(total)+np.sum((np.arange(1,n)[None,:]-powers)*logt,axis=1)
    return alpha,np.exp(logweight),{'hepp_normalization':bound,'primitive_envelope_proved':primitive,
        'gap_rule':'d(A)=|A|-3h(A); sampling gap=max(1,d(A)); exact recursive prefix dynamic program',
        'sector_count_represented':'all edge orderings, sampled by subset dynamic programming'}
