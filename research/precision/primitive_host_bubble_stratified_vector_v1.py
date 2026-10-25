"""Diagnostic exact weighted-edge strata with triangular determinant transport."""
import numpy as np
from disjoint_triangle_bubble_qmc_vector_v1 import bubble
from primitive_host_bubble_vector_v1 import poly


def decoded(points,table,strata,position):
    p=np.asarray(points,dtype=np.longdouble);gamma={int(k):v for k,v in table['gamma'].items()}
    W={key:np.longdouble(str(float(v))) for key,v in strata['weights'].items()};edge=strata['edge'];bit=1<<edge
    n=len(p);remaining=np.full(n,511,dtype=np.int64);aux=p[:,0].copy();order=np.zeros((n,9),dtype=np.int64);exponents=np.ones((n,8),dtype=np.int64)
    for depth in range(9):
        for mask in np.unique(remaining):
            where=np.flatnonzero(remaining==mask);k=position-depth if mask&bit else -1
            edges=[e for e in range(9) if mask&(1<<e) and (k==-1 or ((e==edge)==(k==0)))]
            weights=np.array([W[str(mask^(1<<e))+':'+str(-1 if e==edge or k==-1 else k-1)]/gamma.get(mask^(1<<e),1) for e in edges])
            prob=weights/weights.sum();cdf=np.cumsum(prob);cdf[-1]=1
            chosen=np.searchsorted(cdf,aux[where],side='right');assert (chosen<len(edges)).all()
            left=np.concatenate((np.array([0],dtype=np.longdouble),cdf[:-1]));aux[where]=(aux[where]-left[chosen])/prob[chosen]
            selected=np.array(edges)[chosen];order[where,depth]=selected;remaining[where]=mask^(1<<selected)
            if depth<8:exponents[where,depth]=[gamma[int(v)] for v in remaining[where]]
    assert (order[:,position]==edge).all() and ((aux>=0)&(aux<1)).all()
    assert (exponents[:,0]==2).all() and (exponents[:,-1]==1).all()
    return p,order,exponents,aux


def evaluate_decoded(p,order,gamma,u,binding,table):
    n=len(p)
    middle=p[:,2:8]**(1/gamma[:,1:7].astype(np.longdouble))
    ratios=np.concatenate((np.ones((n,1)),middle,np.ones((n,1))),axis=1)
    base=np.ones_like(p);base[np.arange(n)[:,None],order[:,1:]]=np.cumprod(ratios,axis=1)
    first,last=order[:,0],order[:,-1];coeff=[np.zeros(n,dtype=np.longdouble) for _ in range(4)]
    for mask,co in binding['host_U']:
        value=co*np.prod(base[:,[i for i in range(9) if mask&(1<<i)]],axis=1)
        hasfirst=(mask&(1<<first))!=0;haslast=(mask&(1<<last))!=0
        category=np.where(haslast,np.where(hasfirst,2,3),np.where(hasfirst,0,1))
        for j in range(4):coeff[j]+=np.where(category==j,value,0)
    A,B,C,D=coeff;assert (A>0).all() and all((x>=0).all() for x in (B,C,D))
    # U=t2²[A+B*t2+C*t9+D*t2*t9]. Integrate the Beta(2,1)
    # coordinate first; its t9 marginal is 1/[(A+C*t9)(A+B+(C+D)*t9)²].
    a,b,c,d=A,A+B,C,C+D
    one_third=np.longdouble(1)/3;choose_a=u<one_third
    u=np.where(choose_a,3*u,(3*u-1)/2)
    j=np.where(choose_a,a,b);k=np.where(choose_a,c,d);r=p[:,8]
    v=j+k;root=np.sqrt((1-r)*v*v+r*j*j)
    t9=j*r*(2*j+k)/(root*(v+root))
    Ia=(2*a+c)/(2*a*a*(a+c)**2);Ib=(2*b+d)/(2*b*b*(b+d)**2)
    q9=one_third/((a+c*t9)**3*Ia)+(1-one_third)/((b+d*t9)**3*Ib)
    Ay=A+C*t9;By=B+D*t9;s=np.sqrt(p[:,1]);t2=Ay*s/(Ay+By*(1-s))
    alpha=base.copy();alpha[np.arange(n)[:,None],order[:,1:]]*=t2[:,None]
    alpha[np.arange(n),last]*=t9
    U=poly(binding['host_U'],alpha);F=poly(binding['host_F'],alpha);Q=F/U
    check=t2*t2*(A+B*t2+C*t9+D*t2*t9)
    assert (np.abs(U-check)<=np.longdouble('1e-14')*U).all()
    nu=(np.arange(7,1,-1)[None,:]-gamma[:,1:7])//3
    tropical_rest=np.prod(middle**nu,axis=1)
    weight=np.longdouble(str(float(table['normalization'])))*tropical_rest**3/(Ay*(Ay+By)**2*q9)
    M=alpha.sum(axis=1);edge=alpha[:,8];other=alpha[:,:8].sum(axis=1)
    def kernel(aux):
        mass=M+aux*Q;return -edge*Q*bubble((other+aux*Q)/edge)/(mass*mass)
    value=weight*(kernel(u)+kernel(1-u))/2
    assert np.isfinite(value).all() and (value<=0).all();return value


def evaluate(points,binding,table,strata):
    total=np.zeros(len(points),dtype=np.longdouble)
    for position in range(9):
        p,order,gamma,u=decoded(points,table,strata,position)
        total+=np.longdouble(str(float(strata['masses'][position]/strata['normalization'])))*evaluate_decoded(p,order,gamma,u,binding,table)
    return total
