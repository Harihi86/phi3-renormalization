"""Diagnostic six-scale host transport with an exact external-fraction control."""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble


def parameters(points,table):
    p=np.asarray(points,dtype=np.longdouble);n=len(p);assert p.shape[1]==9
    H={int(k):np.longdouble(str(float(v))) for k,v in table['H'].items()}
    gamma={int(k):v for k,v in table['gamma'].items()};nu={int(k):v for k,v in table['nu'].items()}
    remaining=np.full(n,63,dtype=np.int64);aux=p[:,0].copy();order=np.zeros((n,6),dtype=np.int64)
    exponents=np.ones((n,5),dtype=np.int64);leading=np.zeros((n,5),dtype=np.int64)
    for depth in range(6):
        for mask in np.unique(remaining):
            where=np.flatnonzero(remaining==mask);edges=[e for e in range(6) if mask&(1<<e)]
            w=np.array([H[mask^(1<<e)] for e in edges]);prob=w/w.sum();cdf=np.cumsum(prob);cdf[-1]=1
            chosen=np.searchsorted(cdf,aux[where],side='right');assert (chosen<len(edges)).all()
            left=np.concatenate((np.array([0],dtype=np.longdouble),cdf[:-1]));aux[where]=(aux[where]-left[chosen])/prob[chosen]
            edge=np.array(edges)[chosen];order[where,depth]=edge;remaining[where]=mask^(1<<edge)
            if depth<5:
                exponents[where,depth]=[gamma[int(s)] for s in remaining[where]]
                leading[where,depth]=[nu[int(s)] for s in remaining[where]]
    assert ((aux>=0)&(aux<1)).all()
    ratios=p[:,1:6]**(1/exponents.astype(np.longdouble));a=np.ones((n,6),dtype=np.longdouble)
    a[np.arange(n)[:,None],order[:,1:]]=np.cumprod(ratios,axis=1)
    return a,p[:,6:].copy(),aux,np.prod(ratios**leading,axis=1)


def polynomial(rows,alpha):
    result=np.zeros(len(alpha),dtype=np.longdouble)
    for mask,co in rows:result+=co*np.prod(alpha[:,[e for e in range(alpha.shape[1]) if mask&(1<<e)]],axis=1)
    return result


def evaluate(points,binding,table,control=True,transport=True,reflection=False):
    a,fractions,u,tropical=parameters(points,table);n=len(points);U=polynomial(table['U'],a)
    M=a.sum(axis=1);scale=np.longdouble(str(float(table['normalization'])))*(tropical/U)**3
    def evaluate_fractions(base):
        alpha=np.ones((n,9),dtype=np.longdouble);x=base.copy();jac=np.ones(n,dtype=np.longdouble)
        for row in table['fractions']:
            j=row['index']
            if row['weighted'] and transport:jac*=2*x[:,j];x[:,j]=x[:,j]**2
        for g,edges in enumerate(table['groups']):
            if len(edges)==1:alpha[:,edges[0]]=a[:,g]
            else:
                row=next(row for row in table['fractions'] if row['group']==g);z=x[:,row['index']]
                alpha[:,edges[0]]=a[:,g]*z;alpha[:,edges[1]]=a[:,g]*(1-z)
        Q=polynomial(binding['host_F'],alpha)/U;d=alpha[:,binding['weighted_host_edge']]
        other=alpha[:,[e for e in range(9) if e!=binding['weighted_host_edge']]].sum(axis=1)
        def term(aux):
            mass=M+aux*Q;return -d*Q*bubble((other+aux*Q)/d)/(mass*mass)
        value=(term(u)+term(1-u))/2
        if control:
            average=np.zeros(n,dtype=np.longdouble)
            for mask,co in binding['host_F']:
                monomial=np.full(n,co,dtype=np.longdouble)
                for g,edges in enumerate(table['groups']):
                    count=sum(bool(mask&(1<<e)) for e in edges)
                    if len(edges)==2 and binding['weighted_host_edge'] not in edges:
                        if count==1:monomial*=a[:,g]/2
                        elif count==2:monomial*=a[:,g]**2/6
                    else:
                        for e in edges:
                            if mask&(1<<e):monomial*=alpha[:,e]
                average+=monomial
            value+=d*bubble(other/d)/(M*M)*(Q-average/U)
        return jac*value
    result=evaluate_fractions(fractions)
    if reflection:result=(result+evaluate_fractions(1-fractions))/2
    result*=scale;assert np.isfinite(result).all();return result
