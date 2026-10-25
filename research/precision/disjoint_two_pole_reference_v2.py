"""Independent original four-loop forest matrices, then high-precision radial integration."""
from itertools import permutations
import mpmath as mp
import two_pole_radial_master_v2 as F


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
            routes=[ctx.matrix([[ctx.mpf(x) for x in row]]) for row in state['routing']];us=[];fs=[];effective=[]
            for scale in range(1,5):
                r=[M*scale*y for y in shape];original=[None]*12
                for edge,x in zip(binding['host_edges'],alpha):original[edge]=x
                for i,edges in enumerate(binding['children']):
                    for edge,x in zip(edges,child[i]):original[edge]=r[i]*x
                matrix=sum((x*row.T*row for x,row in zip(original,routes)),ctx.zeros(6));loop=matrix[:4,:4]
                effective.append(matrix[2:,2:]-matrix[2:,:2]*(matrix[:2,:2]**-1)*matrix[:2,2:])
                U=ctx.det(loop)/ctx.fprod(r);schur=matrix[4:,4:]-matrix[4:,:4]*(loop**-1)*matrix[:4,4:]
                Q=sum(gram[i,j]*schur[j,i] for i in range(2) for j in range(2));us.append(U);fs.append(U*Q)
            ud=list(vandermonde*ctx.matrix(us));ff=list(vandermonde*ctx.matrix(fs))
            assert all(abs(x)<ctx.power(10,-digits+100)*max(map(abs,ud)) for x in ud[rank+1:])
            assert all(abs(x)<ctx.power(10,-digits+100)*max(map(abs,ff)) for x in ff[rank+2:])
            if u0 is None:u0=ud[0]
            else:assert abs(ud[0]-u0)<ctx.power(10,-digits+100)*abs(u0)
            dd=[x/ud[0] for x in ud[:rank+1]];nn=[x/(M*ud[0]) for x in ff[:rank+2]]
            # The first two transformed loops are the two independent child loops.
            E0=2*effective[0]-effective[1];E1=effective[1]-effective[0]
            tolerance=ctx.power(10,-digits+150)
            assert ctx.norm(effective[3]-(E0+4*E1))<tolerance*max(1,ctx.norm(effective[3]))
            assert abs(ctx.det(E0[:2,:2])-ud[0])<tolerance*abs(ud[0])
            def pencil(E,mass):
                K=ctx.zeros(5);K[:2,:2]=E[:2,:2];K[2:4,2:4]=ctx.mpf(3)/4*E[:2,:2]
                for j in range(2):
                    K[j,4]=K[4,j]=E[j,2]-E[j,3]/2
                    K[j+2,4]=K[4,j+2]=ctx.mpf(3)/4*E[j,3]
                K[4,4]=mass+E[2,2]+E[3,3]-E[2,3]
                return K
            def eigenvalues(A,B,count):
                lower=ctx.cholesky(A);inverse=lower**-1;whitened=inverse*B*inverse.T
                values=list(ctx.eigsy((whitened+whitened.T)/2,eigvals_only=True))
                bound=tolerance*max(map(abs,values));assert all(x>=-bound for x in values)
                positive=[x for x in values if x>bound];assert len(positive)==count;return positive
            roots=eigenvalues(pencil(E0,M),pencil(E1,M),2*rank+1)
            looproots=eigenvalues(E0[:2,:2],E1[:2,:2],rank)
            for pole in looproots:
                index=min(range(len(roots)),key=lambda j:abs(roots[j]-pole))
                assert abs(roots[index]-pole)<tolerance*abs(pole);roots.pop(index)
            # Independently reconstruct every coefficient of mass*U+F from the pencil roots.
            product=[ctx.mpf(1)]
            for pole in roots:
                out=[ctx.mpf(0)]*(len(product)+1)
                for j,co in enumerate(product):out[j]+=co;out[j+1]+=co*pole
                product=out
            for j,co in enumerate(product):
                wanted=(dd[j] if j<len(dd) else 0)+(dd[j-1] if j else 0)+(nn[j] if j<len(nn) else 0)
                assert abs((1+nn[0])*co-wanted)<tolerance*max(1,abs(wanted))
            value=F.primitive(ctx,nn[0],dd,nn,roots)
            forest+=(1 if active==3 else -1)*value
        total+=(tropical/u0)**3*forest/(8*ctx.fprod(shape))
    return ctx,total
