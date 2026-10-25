"""Diagnostic-only full-domain Hepp transport; no interval or acceptance claims."""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble


def parameters(points,table):
    p=np.asarray(points,dtype=np.longdouble);assert p.ndim==2 and p.shape[1]==9 and ((p>0)&(p<1)).all()
    H={int(k):np.longdouble(str(float(v))) for k,v in table['H'].items()}
    gamma={int(k):v for k,v in table['gamma'].items()}
    n=len(p);remaining=np.full(n,511,dtype=np.int64);aux=p[:,0].copy();order=np.zeros((n,9),dtype=np.int64)
    exponents=np.ones((n,8),dtype=np.int64)
    for depth in range(9):
        for mask in np.unique(remaining):
            where=np.flatnonzero(remaining==mask);edges=[e for e in range(9) if mask&(1<<e)]
            weights=np.array([H[mask^(1<<e)] for e in edges]);prob=weights/weights.sum();cdf=np.cumsum(prob);cdf[-1]=1
            selected=np.searchsorted(cdf,aux[where],side='right');assert (selected<len(edges)).all()
            left=np.concatenate((np.array([0],dtype=np.longdouble),cdf[:-1]))
            aux[where]=(aux[where]-left[selected])/prob[selected]
            chosen=np.array(edges)[selected];order[where,depth]=chosen;remaining[where]=mask^(1<<chosen)
            if depth<8:exponents[where,depth]=[gamma[int(s)] for s in remaining[where]]
    assert ((aux>=0)&(aux<1)).all()
    ratios=p[:,1:]**(1/exponents.astype(np.longdouble));alpha=np.ones_like(p)
    products=np.cumprod(ratios,axis=1)
    alpha[np.arange(n)[:,None],order[:,1:]]=products
    sizes=np.arange(8,0,-1);nu=(sizes[None,:]-exponents)//3
    assert np.array_equal(sizes[None,:]-exponents,3*nu)
    tropical=np.prod(ratios**nu,axis=1)
    return alpha,aux,tropical


def poly(rows,alpha):
    value=np.zeros(len(alpha),dtype=np.longdouble)
    for mask,co in rows:value+=co*np.prod(alpha[:,[e for e in range(9) if mask&(1<<e)]],axis=1)
    return value


def evaluate(points,binding,table,antithetic=True):
    a,u,tropical=parameters(points,table);U=poly(binding['host_U'],a);F=poly(binding['host_F'],a);Q=F/U
    M=a.sum(axis=1);d=a[:,binding['weighted_host_edge']];remaining=M-d
    # Sum unweighted masses directly when d dominates, avoiding M-d cancellation.
    remaining=a[:,:8].sum(axis=1)
    def f(aux):
        mass=M+aux*Q;return -d*Q*bubble((remaining+aux*Q)/d)/(mass*mass)
    value=f(u)
    if antithetic:value=(value+f(1-u))/2
    value*=np.longdouble(str(float(table['normalization'])))*(tropical/U)**3
    assert np.isfinite(value).all() and (value<=0).all()
    return value
