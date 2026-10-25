"""Independent scalar differentiation, formal products and 72-frame fourth moments."""
import mpmath as mp
import nested_two_bubble_reference_v1 as B
M=B.M


def density(ctx,p,binding,quadrature=False):
    result=B.density(ctx,p,binding,True,quadrature)
    s,z,theta,x0,y=p;t=(s/(1-s))**2;A=1+t
    den=x0*x0+(1-x0)**2;x=x0*x0/den;jac=2*x0*(1-x0)/den**2
    parent=binding['parent_boundary_shifts'];parent=[*parent,[a-b for a,b in zip(*parent)]]
    host=binding['host_distinct_routes'];br=binding['disjoint_bubble_boundary']
    norm=lambda r:r[1]**2+r[2]**2-r[1]*r[2]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    V0=sum(w*r[0]**2*t for w,r in zip(weights,parent));Cp=sum(w*norm(r) for w,r in zip(weights,parent))
    Kprime=lambda v:-x*B.bubble(ctx,(1-x+v)/x)/(1+v)**2
    K=[M.bubble_K(ctx,x,V0)]+[ctx.diff(Kprime,V0,n-1)/ctx.factorial(n) for n in range(1,5)]
    jfun=lambda v:B.bubble(ctx,v)/(1+v)
    J=[ctx.diff(jfun,t,n)/ctx.factorial(n) for n in range(5)]
    def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(5)]
    def compose(a,l,c):
        out=[ctx.mpf(0)]*5;power=[ctx.mpf(1)]+[ctx.mpf(0)]*4;delta=[ctx.mpf(0),l,c,ctx.mpf(0),ctx.mpf(0)]
        for coefficient in a:
            out=[v+coefficient*w for v,w in zip(out,power)];power=mul(power,delta)
        return out
    def coefficient(qx,qy):
        q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
        linear=lambda r:2*r[0]*(r[1]*q1+r[2]*q2)
        product=[ctx.mpf(1)]+[ctx.mpf(0)]*4
        for r in host:
            l,c=linear(r)/A,ctx.mpf(norm(r))/A
            product=mul(product,[ctx.mpf(1),-l,l*l-c,-l**3+2*l*c,l**4-3*l*l*c+c*c])
        lp=sum(w*linear(r) for w,r in zip(weights,parent))
        product=mul(product,compose(K,lp,Cp));product=mul(product,compose(J,linear(br),ctx.mpf(norm(br))))
        return product[4]
    point=coefficient(ctx.sqrt(t*z)*ctx.cospi(theta),ctx.sqrt(t*z)*ctx.sinpi(theta))
    # Independent six-dimensional axis/pair cubature, exact for degree<=4.
    # Twelve axes have weight -1/48; sixty pair directions have weight +1/48.
    mean=ctx.mpf(0);radius=ctx.sqrt(t)
    for i in range(6):
        for sign in (-1,1):
            q=[ctx.mpf(0)]*6;q[i]=sign*radius;mean-=coefficient(q[0],q[1])/48
    for i in range(6):
        for j in range(i+1,6):
            for si in (-1,1):
                for sj in (-1,1):
                    q=[ctx.mpf(0)]*6;q[i]=si*radius/ctx.sqrt(2);q[j]=sj*radius/ctx.sqrt(2);mean+=coefficient(q[0],q[1])/48
    common=2/(s*(1-s))*(t/A)**3*(1-x)*jac
    return result+common*(mean/2-(1-z)*point)


def evaluate(point,binding,digits=640,mode='causal',quadrature=False):
    assert mode in ('causal','causal-balanced') and len(point)==5
    ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
    u,v=p[3:];du=u*u+(1-u)**2;xi=u*u/du;jxi=2*u*(1-u)/du**2;w=v;jw=ctx.mpf(1)
    if mode=='causal-balanced':
        dv=v*v+(1-v)**2;w=v*v/dv;jw=2*v*(1-v)/dv**2
    x=(1-xi)*w;y=xi/(1-x);x0=ctx.sqrt(x)/(ctx.sqrt(x)+ctx.sqrt(1-x))
    oldden=x0*x0+(1-x0)**2;oldjac=2*x0*(1-x0)/oldden**2
    transport_jac=jxi*jw*(1-xi)/((1-x)*oldjac);p[3]=x0;p[4]=y
    assert 0<x0<1 and 0<y<1
    s=p[0];x0=p[3];y=p[4];x=x0*x0/(x0*x0+(1-x0)**2)
    a,b=binding['parent_boundary_shifts'];loops=[a[0],b[0],a[0]-b[0]]
    weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
    ell=sum(w*l*l for w,l in zip(weights,loops));h=1/ctx.sqrt(ell);den=1-s+h*s
    pp=p.copy();pp[0]=h*s/den
    return ctx,transport_jac*h/(den*den)*density(ctx,pp,binding,quadrature)
