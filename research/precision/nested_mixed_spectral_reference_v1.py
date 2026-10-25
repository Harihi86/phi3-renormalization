"""Independent original-Gaussian mass differentiation and direct radial integration."""
from itertools import combinations
import mpmath as mp
import g521_mixture_reference_v1 as M


def data(point,binding,ctx):
    p=[ctx.mpf(q.numerator)/q.denominator for q in point]
    r=ctx.sqrt(1-p[0]);alpha=[p[0]/(1+r),r*p[1],r*(1-p[1])]
    child=[];jac=ctx.mpf(1)
    for i in range(2):
        u,v=p[2+2*i:4+2*i];s,c=ctx.sinpi(v/2),ctx.cospi(v/2);rem=1-u*u
        child.append([u*u,rem*s*s,rem*c*c]);jac*=2*ctx.pi*u*rem*s*c
    sn,cs=ctx.sinpi(p[6]/2),ctx.cospi(p[6]/2);y=[sn*sn,cs*cs]
    norm=lambda v:v[0]**2+v[1]**2-v[0]*v[1]
    dot=lambda v,w:v[0]*w[0]+v[1]*w[1]-(v[0]*w[1]+v[1]*w[0])/2
    host=sum((a*ctx.matrix(row)*ctx.matrix([row]) for a,row in zip(alpha,binding['host_routes'])),ctx.zeros(3))
    b=[host[0,1],host[0,2]];bc=[sum(row[0]*row[j] for row in binding['host_routes'])/ctx.mpf(3) for j in (1,2)]
    q=host[1,1]+host[2,2]-host[1,2]-norm(b);parts=[]
    for i in range(2):
        matrix=ctx.zeros(3)
        for (a,c),row in zip(combinations(range(3),2),binding['child_pair_differences'][i]):
            matrix+=y[i]*child[i][a]*child[i][c]*ctx.matrix(row)*ctx.matrix([row])
        parts.append(matrix)
    corners=[]
    for mask in range(1,4):
        matrix=sum((m for i,m in enumerate(parts) if mask&(1<<i)),ctx.zeros(3))
        h=matrix[0,0];cross=[matrix[0,1],matrix[0,2]];energy=matrix[1,1]+matrix[2,2]-matrix[1,2]
        de=energy-2*dot(b,cross)+h*norm(b);ge=h*energy-norm(cross)
        mean=energy-2*dot(bc,cross)+h*norm(bc)+h/12
        corners.append((mask,q,h,de,ge,mean))
    s=p[7];a=(1-s*s)/4;x=child[0][0]
    return {'factor':ctx.pi*s**4*jac/(24*sn*cs),'a':a,'x':x,'y':y[0],
            'mass':1+y[0]*x*(1/a-1),'corners':corners,'q':q}


def evaluate(point,binding,control=True,digits=300):
    ctx=mp.mp.clone();ctx.dps=digits;d=data(point,binding,ctx);a,x,y=d['a'],d['x'],d['y']
    def forest(mu):
        return sum((-1)**(2-mask.bit_count())*M.closed(ctx,q,h/mu,de/mu,ge/mu**2,mean/mu,control)
                   for mask,q,h,de,ge,mean in d['corners'])
    unit=forest(ctx.mpf(1));derivative=ctx.diff(forest,ctx.mpf(1))
    answer=(a-2)/(1-a)**2*unit+forest(d['mass'])/(a*(1-a)**2)-x*y/(1-a)*derivative
    return ctx,d['factor']*answer


def direct_radial(point,binding,digits=90):
    """Integrate the four Gaussian forest corners before using any dilogarithm."""
    ctx=mp.mp.clone();ctx.dps=digits;d=data(point,binding,ctx);a,x,y=d['a'],d['x'],d['y']
    def integrand(r):
        def forest(mu,derivative=False):
            mass=1+mu*r;value=ctx.mpf(0)
            for mask,q,h,de,ge,mean in [(0,d['q'],0,0,0,0),*d['corners']]:
                U=1+h*r;F=q+(h*q+de)*r+ge*r*r;Q=F/U
                v=-r*Q/(mass*(mass+Q)*U**3) if derivative else ctx.log1p(Q/mass)/U**3
                value-=(-1)**(2-mask.bit_count())*v
            return value
        return ((a-2)*forest(1)+forest(d['mass'])/a)/(1-a)**2/r-x*y/(1-a)*forest(1,True)/r
    value=ctx.quad(integrand,[0,ctx.mpf('.1'),1,10,100,ctx.inf])
    return ctx,d['factor']*value
