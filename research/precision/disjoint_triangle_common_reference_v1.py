"""Independent full four-loop matrix reference on the common radial chart."""
from itertools import permutations,combinations
import mpmath as mp


def parameters(point,binding,positive,ctx,order,reflect=False,transport=True):
    p=[ctx.mpf(x.numerator)/x.denominator for x in point];L=positive['host_loops'];k=4-L
    groups=positive['host_measure']['groups'];rows=[ctx.matrix([r]) for r in binding['host_routes']]
    def simplex(u,v):
        root=ctx.sqrt(1-u);return [u/(1+root),root*v,root*(1-v)]
    if L==2:
        a=[ctx.mpf(1)]*3;a[order[1]]=p[0];a[order[2]]=p[0]*ctx.sqrt(p[1]);alpha=[None]*6
        for g,(e,f) in enumerate(groups):alpha[e]=a[g]*p[2+g];alpha[f]=a[g]*(1-p[2+g])
        tropical=a[order[0]]*a[order[1]];base=5;shape=[p[9],1-p[9]]
    else:
        a=None;alpha=simplex(p[0],p[1]);tropical=ctx.mpf(1);base=2;shape=simplex(p[8],p[9])
    child=[simplex(p[base+2*i],p[base+2*i+1]) for i in range(k)]
    host=sum((x*r.T*r for x,r in zip(alpha,rows)),ctx.zeros(L+2));inverse=host[:L,:L]**-1
    ell=[]
    for i in range(k):
        value=ctx.mpf(0)
        for (aa,b),row in zip(combinations(range(3),2),binding['child_pair_differences'][i]):
            loop=ctx.matrix([row[:L]]);value+=child[i][aa]*child[i][b]*(loop*inverse*loop.T)[0]
        ell.append(value)
    M=sum(alpha);trace=sum(y*e for y,e in zip(shape,ell))
    if transport==2:trace=1/sum(1/(y*e) for y,e in zip(shape,ell))
    h=ctx.log(M*trace);u=1-p[10] if reflect else p[10]
    def transform(z):
        if not transport or h==0:return M*z/(1-z)
        return M*ctx.expm1(z*h)/(ctx.exp(z*h)*ctx.expm1((1-z)*h))
    radial=transform(u);derivative=ctx.diff(transform,u)
    jac=radial**(k-1)*derivative/(2 if k==3 else 1)
    r=[radial*y for y in shape];volume=ctx.mpf(1)/(2**k*(2 if L==1 else 1))
    return alpha,child,r,M+radial,jac,volume,tropical,a


def evaluate(point,binding,positive,control=True,reflection=2,digits=300,transport=True):
    ctx=mp.mp.clone();ctx.dps=digits;gram=ctx.matrix([[1,ctx.mpf('-.5')],[ctx.mpf('-.5'),1]])
    L=positive['host_loops'];k=4-L;groups=positive['host_measure']['groups']
    orders=list(permutations(range(3))) if L==2 else [(0,1,2)]
    rows=[ctx.matrix([r]) for r in binding['host_routes']];total=ctx.mpf(0)
    for order in orders:
        for reflect in range(reflection):
            alpha,child,r,mass,jac,volume,tropical,a=parameters(point,binding,positive,ctx,order,bool(reflect),transport)
            original_alpha=[None]*12
            for e,val in zip(binding['host_edges'],alpha):original_alpha[e]=val
            for i,edges in enumerate(binding['children']):
                for e,x in zip(edges,child[i]):original_alpha[e]=r[i]*x
            forest=ctx.mpf(0)
            for state in binding['states']:
                route=[ctx.matrix([[ctx.mpf(v) for v in row]]) for row in state['routing']]
                matrix=sum((v*row.T*row for v,row in zip(original_alpha,route)),ctx.zeros(6))
                loop=matrix[:4,:4];U=ctx.det(loop)/ctx.fprod(r)
                schur=matrix[4:,4:]-matrix[4:,:4]*(loop**-1)*matrix[:4,4:]
                Q=sum(gram[i,j]*schur[j,i] for i in range(2) for j in range(2))
                value=ctx.log1p(Q/mass)/U**3
                if control:
                    center=[ctx.mpf(1)/3]*3 if L==1 else [ctx.mpf(0)]*6
                    if L==2:
                        for g,edges in enumerate(groups):
                            for e in edges:center[e]=a[g]/2
                    cm=sum((v*row.T*row for v,row in zip(center,rows)),ctx.zeros(L+2))
                    for i in set(range(k))-set(state['word']):
                        for (aa,b),rr in zip(combinations(range(3),2),binding['child_pair_differences'][i]):
                            row=ctx.matrix([rr]);cm+=r[i]*child[i][aa]*child[i][b]*row.T*row
                    ci=cm[:L,:L]**-1;cs=cm[L:,L:]-cm[L:,:L]*ci*cm[:L,L:]
                    mean=sum(gram[i,j]*cs[j,i] for i in range(2) for j in range(2))
                    if L==1:mean-=ctx.mpf(1)/(12*ctx.det(cm[:L,:L]))
                    else:mean-=sum(a[g]**2*(rows[edges[0]][:,:L]*ci*rows[edges[0]][:,:L].T)[0]/12 for g,edges in enumerate(groups))
                    value+=(mean-Q)/(mass*U**3)
                forest-=state['sign']*value
            z=ctx.mpf(3) if L==2 else ctx.mpf(1)
            total+=z*tropical**3*forest*jac*volume/(ctx.fprod(r)*len(orders)*reflection)
    return ctx,total
