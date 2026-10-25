"""Independent mpmath common-radial primitive for the three-child L=1 host."""
from itertools import combinations
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
        # Independent rational-moment integration via hypergeometric values.
        moments=[ctx.hyp2f1(1,n,n+1,1-h)/n for n in range(1,5)]
        value+=(ctx.mpf(1)/4-q)*h*sum(moments[:3])-(mean_de-de)*moments[3]
    return value


def data(point,binding,ctx):
    assert binding['host_loops']==1 and len(binding['children'])==3
    p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    def simplex(u,v):
        s=ctx.sqrt(1-u);return [u/(1+s),s*v,s*(1-v)]
    alpha=simplex(p[0],p[1]);child=[simplex(p[2+2*i],p[3+2*i]) for i in range(3)];y=simplex(p[8],p[9])
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
    return y,corners


def evaluate(point,binding,control=True,digits=300):
    ctx=mp.mp.clone();ctx.dps=digits;y,corners=data(point,binding,ctx)
    value=sum((-1)**(3-mask.bit_count())*closed(ctx,q,h,de,ge,mean,control) for mask,q,h,de,ge,mean in corners)
    # Three child simplex areas, one host area, one radial-shape area.
    return ctx,value/(32*ctx.fprod(y))
