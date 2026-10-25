"""Independent original four-loop forest matrices, then high-precision radial integration."""
from itertools import permutations
import mpmath as mp
import two_pole_radial_master_v1 as F


def evaluate(point,binding,positive,digits=720):
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    gram=ctx.matrix([[1,ctx.mpf('-.5')],[ctx.mpf('-.5'),1]])
    def simplex(u,v):
        t=ctx.sqrt(1-u);return [u/(1+t),t*v,t*(1-v)]
    child=[simplex(p[5],p[6]),simplex(p[7],p[8])];shape=[p[9],1-p[9]]
    child_ranks=[]
    for mask in (1,2,3):
        vectors=[row[:2] for i,rows in enumerate(binding['child_pair_differences']) if mask&(1<<i) for row in rows if any(row[:2])]
        child_ranks.append(2 if any(a[0]*b[1]!=a[1]*b[0] for a in vectors for b in vectors) else 1)
    total=ctx.mpf(0);vandermonde=ctx.matrix([[ctx.mpf(s)**j for j in range(4)] for s in range(1,5)])**-1
    for order in permutations(range(3)):
        a=[ctx.mpf(1)]*3;a[order[1]]=p[0];a[order[2]]=p[0]*ctx.sqrt(p[1]);alpha=[None]*6
        for j,(e,f) in enumerate(positive['host_measure']['groups']):alpha[e]=a[j]*p[2+j];alpha[f]=a[j]*(1-p[2+j])
        M=sum(alpha);tropical=a[order[0]]*a[order[1]];forest=ctx.mpf(0);u0=None
        for active,rank in zip((1,2,3),child_ranks):
            state=next(row for row in binding['states'] if set(row['word'])==set(i for i in range(2) if not active&(1<<i)))
            routes=[ctx.matrix([[ctx.mpf(x) for x in row]]) for row in state['routing']];us=[];fs=[]
            for scale in range(1,5):
                r=[M*scale*y for y in shape];original=[None]*12
                for edge,x in zip(binding['host_edges'],alpha):original[edge]=x
                for i,edges in enumerate(binding['children']):
                    for edge,x in zip(edges,child[i]):original[edge]=r[i]*x
                matrix=sum((x*row.T*row for x,row in zip(original,routes)),ctx.zeros(6));loop=matrix[:4,:4]
                U=ctx.det(loop)/ctx.fprod(r);schur=matrix[4:,4:]-matrix[4:,:4]*(loop**-1)*matrix[:4,4:]
                Q=sum(gram[i,j]*schur[j,i] for i in range(2) for j in range(2));us.append(U);fs.append(U*Q)
            ud=list(vandermonde*ctx.matrix(us));ff=list(vandermonde*ctx.matrix(fs))
            assert all(abs(x)<ctx.power(10,-digits+100)*max(map(abs,ud)) for x in ud[rank+1:])
            assert all(abs(x)<ctx.power(10,-digits+100)*max(map(abs,ff)) for x in ff[rank+2:])
            if u0 is None:u0=ud[0]
            else:assert abs(ud[0]-u0)<ctx.power(10,-digits+100)*abs(u0)
            dd=[x/ud[0] for x in ud[:rank+1]];nn=[x/(M*ud[0]) for x in ff[:rank+2]]
            value=F.primitive(ctx,nn[0],dd,nn)
            forest+=(1 if active==3 else -1)*value
        total+=(tropical/u0)**3*forest/(8*ctx.fprod(shape))
    return ctx,total
