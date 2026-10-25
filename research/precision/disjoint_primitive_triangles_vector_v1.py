"""Long-double diagnostics of the complete factored cube, never acceptance."""
from itertools import combinations, permutations
from math import factorial
import numpy as np
import forest_difference_algebra_v1 as A


def simplex(p):
    root=np.sqrt(1-p[:,0])
    return np.stack((p[:,0]/(1+root),root*p[:,1],root*(1-p[:,1])),axis=1)


def coefficients(rows,alpha,child,r,groups=None,scales=None,mean=False):
    n=len(alpha);k=len(r);count=alpha.shape[1]
    out=[np.zeros(n,dtype=np.longdouble) for _ in range(1<<k)]
    for mask,co in rows:
        term=np.full(n,co,dtype=np.longdouble);bits=0
        if mean:
            if count==3:term*=np.longdouble(2)/factorial((mask&7).bit_count()+2)
            else:
                for g,edges in enumerate(groups):
                    degree=sum(bool(mask&(1<<e)) for e in edges)
                    if degree:term*=scales[:,g]**degree/(2 if degree==1 else 6)
        else:
            for e in range(count):
                if mask&(1<<e):term*=alpha[:,e]
        for i in range(k):
            degree=0
            for j,(a,b) in enumerate(((0,1),(0,2),(1,2))):
                if mask&(1<<(count+3*i+j)):
                    term*=child[i][:,a]*child[i][:,b];degree+=1
            if degree:
                bits|=1<<i
                if degree>1:term*=r[i]**(degree-1)
        out[bits]+=term
    return out


def divided_log(t,positive,negative):
    t=np.asarray(t,dtype=np.longdouble);out=np.ones_like(t)
    small=(t!=0)&(np.abs(t)<np.longdouble('.25'))
    large=(t!=0)&~small
    out[small]=np.log1p(t[small])/t[small]
    out[large]=(np.log(positive[large])-np.log(negative[large]))/t[large]
    return out


def kernel(U,F,mass,r,meanF=None):
    inv=A.inverse(U,r)
    # Rescale P by M before the logarithm. Constants are log1p, never a
    # subtraction of three nearly equal ordinary logarithms.
    P=A.add(U,A.scale(F,1/mass))
    logs=A.subtract(A.logarithm(P,r,np.log,divided_log),A.logarithm(U,r,np.log,divided_log))
    logs[0]=np.log1p(F[0]/(mass*U[0]))
    result=A.multiply(A.power(inv,3,r),logs,r)
    if meanF is not None:
        result=A.add(result,A.scale(A.multiply(A.subtract(meanF,F),A.power(inv,4,r),r),1/mass))
    return -result[-1]


def parameters(points,binding,positive,order,transport=True,reflect=False):
    p=np.asarray(points,dtype=np.longdouble);n=len(p);L=positive['host_loops'];k=4-L
    routes=np.asarray(binding['host_routes'],dtype=np.longdouble)
    pairs=np.asarray(binding['child_pair_differences'],dtype=np.longdouble)
    groups=positive['host_measure']['groups']
    if L==2:
        # Largest a is one; suffix gammas are 1 and 2. Six equal order
        # masses are summed, so no discontinuous random order selection.
        a=np.ones((n,3),dtype=np.longdouble)
        a[:,order[1]]=p[:,0];a[:,order[2]]=p[:,0]*np.sqrt(p[:,1])
        alpha=np.zeros((n,6),dtype=np.longdouble)
        for g,(e,f) in enumerate(groups):
            alpha[:,e]=a[:,g]*p[:,2+g];alpha[:,f]=a[:,g]*(1-p[:,2+g])
        U0=a[:,0]*a[:,1]+a[:,0]*a[:,2]+a[:,1]*a[:,2]
        tropical=a[:,order[0]]*a[:,order[1]]
        base=5;radial_start=9
    else:
        alpha=simplex(p[:,:2]);a=None;U0=np.ones(n,dtype=np.longdouble)
        tropical=U0;base=2;radial_start=8
    M=alpha.sum(axis=1)
    child=[simplex(p[:,base+2*i:base+2*i+2]) for i in range(k)]
    host=np.einsum('ne,ei,ej->nij',alpha,routes[:,:L],routes[:,:L])
    if L==1:inverse=1/host
    else:
        inverse=np.empty_like(host)
        inverse[:,0,0]=host[:,1,1]/U0;inverse[:,1,1]=host[:,0,0]/U0
        inverse[:,0,1]=inverse[:,1,0]=-host[:,0,1]/U0
    r=[];jac=np.ones(n,dtype=np.longdouble)
    for i,x in enumerate(child):
        ell=np.zeros(n,dtype=np.longdouble)
        for j,(aa,b) in enumerate(combinations(range(3),2)):
            v=pairs[i,j,:L]
            ell+=x[:,aa]*x[:,b]*np.einsum('i,nij,j->n',v,inverse,v)
        assert (ell>0).all()
        u=p[:,radial_start+i] if not reflect else 1-p[:,radial_start+i]
        ri=M*u/(1-u);ji=M/(1-u)**2
        if transport:
            h=np.log(M*ell);mask=h!=0;uh=u[mask]*h[mask];vh=(1-u[mask])*h[mask]
            ri[mask]=M[mask]*np.expm1(uh)/(np.exp(uh)*np.expm1(vh))
            # The cube already divided by all r_i: multiply dr_i/du_i.
            ji[mask]=ri[mask]*h[mask]*np.expm1(h[mask])/(np.expm1(uh)*np.expm1(vh))
        r.append(ri);jac*=ji
    volume=np.longdouble(1)/(2**k)
    if L==1:volume/=2
    return alpha,child,r,M+sum(r),jac,volume,tropical,U0,a,groups


def evaluate(points,binding,positive,control=True,transport=True,radial_reflection=1):
    n=len(points);total=np.zeros(n,dtype=np.longdouble);L=positive['host_loops']
    orders=list(permutations(range(3))) if L==2 else [(0,1,2)]
    assert radial_reflection in (1,2)
    for order in orders:
        for reflect in range(radial_reflection):
            alpha,child,r,mass,jac,volume,tropical,U0,a,groups=parameters(points,binding,positive,order,transport,bool(reflect))
            U=coefficients(positive['effective_U'],alpha,child,r)
            F=coefficients(positive['effective_F'],alpha,child,r)
            mean=coefficients(positive['effective_F'],alpha,child,r,groups,a,True) if control else None
            # Normalize both polynomials by U0, preserving F/U; restore
            # U0^-3 with the tropical importance ratio exactly once.
            U=A.scale(U,1/U0);F=A.scale(F,1/U0)
            if mean is not None:mean=A.scale(mean,1/U0)
            z=np.longdouble(3) if L==2 else np.longdouble(1)
            total+=z*(tropical/U0)**3*kernel(U,F,mass,r,mean)*jac*volume/(len(orders)*radial_reflection)
    assert np.isfinite(total).all();return total
