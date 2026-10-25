"""Independent mpmath common-radial primitive for the three-child L=1 host."""
from itertools import combinations,permutations
import mpmath as mp


def phi(ctx,t):
    z=t-1
    if not z:return ctx.pi**2/6-ctx.mpf(5)/4
    if abs(z)<ctx.mpf('.25'):
        count=int(ctx.ceil((ctx.dps+12)*ctx.log(10)/(-ctx.log(abs(z)))))+2
        result=ctx.pi**2/6-ctx.mpf(5)/4;power=z
        for j in range(1,count+1):
            result+=(-1)**(j-1)*2*power/(j*j*(j+1)*(j+2));power*=z
        tail=2*abs(power)/((count+1)**2*(count+2)*(count+3)*(1-abs(z)))
        assert tail<ctx.power(10,-ctx.dps-10)
        return result
    return ctx.pi**2/6-ctx.polylog(2,1-t)+t*(2-3*t)*ctx.log(t)/(2*(1-t)**2)-t/(2*(1-t))


def closed(ctx,q,h,de,ge,mean_de,control):
    a=1+q;b=1+h*a+de;c=h+ge;disc=(1-h*a+de)**2+4*a*(h*de-ge)
    plus=(b+ctx.sqrt(disc))/(2*a);minus=2*c/(b+ctx.sqrt(disc));la=ctx.log1p(q)
    value=la*ctx.log(h)+la*la/2+3*la/2-phi(ctx,plus/h)-phi(ctx,minus/h)+phi(ctx,1/h)+phi(ctx,ctx.mpf(1))
    if control:
        # Closed antiderivative of the rational moments; no general hypergeometric solver.
        z=1-h;moments=[]
        for n in range(1,5):
            if not z:moment=ctx.mpf(1)/n
            elif abs(z)<ctx.mpf('.5'):
                count=int(ctx.ceil((ctx.dps+12)*ctx.log(10)/(-ctx.log(abs(z)))))+2
                moment=sum(z**j/ctx.mpf(n+j) for j in range(count+1))
                tail=abs(z)**(count+1)/((n+count+1)*(1-abs(z)))
                assert tail<ctx.power(10,-ctx.dps-10)
            else:moment=(-ctx.log(h)-sum(z**j/ctx.mpf(j) for j in range(1,n)))/z**n
            moments.append(moment)
        value+=(ctx.mpf(1)/4-q)*h*sum(moments[:3])-(mean_de-de)*moments[3]
    return value


def data(point,binding,ctx,shape,permutation=(0,1,2),mixture_branch=-1,host_override=None):
    assert binding['host_loops']==1 and len(binding['children'])==3
    p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    def simplex(u,v):
        s=ctx.sqrt(1-u);return [u/(1+s),s*v,s*(1-v)]
    alpha=simplex(p[0],p[1]) if host_override is None else list(host_override);child=[simplex(p[2+2*i],p[3+2*i]) for i in range(3)];y=simplex(p[8],p[9])
    assert shape in (3,5)
    childjac=ctx.mpf(1)
    for i in range(3):
        u,v=p[2+2*i:4+2*i];theta=ctx.pi*v/2;rem=(1-u)*(1+u)
        child[i]=[u*u,rem*ctx.sin(theta)**2,rem*ctx.cos(theta)**2]
        childjac*=2*ctx.pi*u*rem*ctx.sin(theta)*ctx.cos(theta)
    theta=ctx.pi*p[9]/2;rem=(1-p[8])*(1+p[8]);w=[p[8]**2,rem*ctx.sin(theta)**2,rem*ctx.cos(theta)**2]
    density_factor=ctx.pi/(8*ctx.sqrt(ctx.fprod(w)))
    density_factor*=8*childjac
    w=[w[i] for i in permutation]
    ell=[sum(child[i][a]*child[i][b]*row[0]**2 for (a,b),row in zip(combinations(range(3),2),binding['child_pair_differences'][i])) for i in range(3)]
    assert all(v>0 for v in ell)
    if shape==5:y=w
    else:
        scaled=[a/b for a,b in zip(w,ell)];total=sum(scaled);y=[a/total for a in scaled]
    E=sum(a*b for a,b in zip(ell,y));ratio=ctx.sqrt(ctx.fprod(ell)/E**3)
    density_factor*=(ratio if mixture_branch else 1)/(1+ratio)
    rows=[ctx.matrix([row]) for row in binding['host_routes']];host=sum((a*r.T*r for a,r in zip(alpha,rows)),ctx.zeros(3))
    b=[host[0,1],host[0,2]];bc=[sum(row[0]*row[j] for row in binding['host_routes'])/ctx.mpf(3) for j in (1,2)]
    norm=lambda v:v[0]**2+v[1]**2-v[0]*v[1]
    dot=lambda v,w:v[0]*w[0]+v[1]*w[1]-(v[0]*w[1]+v[1]*w[0])/2
    q=host[1,1]+host[2,2]-host[1,2]-norm(b);parts=[]
    for i in range(3):
        matrix=ctx.zeros(3)
        for (a,c),row in zip(combinations(range(3),2),binding['child_pair_differences'][i]):
            v=ctx.matrix([row]);matrix+=y[i]*child[i][a]*child[i][c]*v.T*v
        parts.append(matrix)
    corners=[]
    for mask in range(1,8):
        matrix=sum((part for i,part in enumerate(parts) if mask&(1<<i)),ctx.zeros(3))
        h=matrix[0,0];c=[matrix[0,1],matrix[0,2]];e=matrix[1,1]+matrix[2,2]-matrix[1,2]
        de=e-2*dot(b,c)+h*norm(b);ge=h*e-norm(c);mean=e-2*dot(bc,c)+h*norm(bc)+h/12
        corners.append((mask,q,h,de,ge,mean))
    return density_factor,corners


def mixture(point,binding,ctx,control=True,host_override=None):
    result=ctx.mpf(0)
    for branch,kind in ((0,5),(1,3)):
        factor,corners=data(point,binding,ctx,kind,mixture_branch=branch,host_override=host_override)
        result+=factor*sum((-1)**(3-mask.bit_count())*closed(ctx,q,h,de,ge,mean,control) for mask,q,h,de,ge,mean in corners)
    return result


def evaluate(point,binding,control=True,digits=300,shape=8):
    from math import factorial
    assert shape in (8,9);degree=shape-6
    ctx=mp.mp.clone();ctx.dps=digits
    u,v=[ctx.mpf(q.numerator)/q.denominator for q in point[:2]]
    a=ctx.sqrt(1-u);alpha=[u/(1+a),a*v,a*(1-v)]
    result=mixture(point,binding,ctx,control)
    def basis_term(value,n):return ctx.fprod(degree*value-j for j in range(n))/factorial(n)
    for i in range(degree+1):
        for j in range(degree+1-i):
            k=degree-i-j;node=[ctx.mpf(t)/degree for t in (i,j,k)]
            basis=ctx.fprod(basis_term(value,n) for value,n in zip(alpha,(i,j,k)))
            if degree==2:weight=ctx.mpf(0) if max(i,j,k)==2 else ctx.mpf(1)/3
            else:weight=ctx.mpf(1)/30 if max(i,j,k)==3 else (ctx.mpf(9)/20 if min(i,j,k)==1 else ctx.mpf(3)/40)
            result+=(weight-basis)*mixture(point,binding,ctx,control,node)
    return ctx,result
