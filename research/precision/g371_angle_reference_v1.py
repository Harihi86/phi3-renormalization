"""Original routes, Taylor convolution, twelve frames and direct angular oracle."""
import mpmath as mp
import nested_triangle_complete_v1 as H
import nested_analytic_marginals_v1 as M
import g371_angle_marginal_v1 as A


def evaluate(point,binding,digits=640,mode='six-pair',quadrature=False):
    assert binding['id']=='L4E3G371' and mode in ('six-direct','six-pair','six-scale')
    assert len(point)==6 and all(0<x<1 for x in point)
    ctx=mp.mp.clone();ctx.dps=digits
    p=[ctx.mpf(x.numerator)/x.denominator for x in point]
    radial=lambda x:(x/(1-x))**2
    if mode=='six-scale':
        scale=radial(p[0]);t,v=scale*p[1],scale*(1-p[1]);jac=scale*ctx.diff(radial,p[0])
    else:t,v=radial(p[0]),radial(p[1]);jac=ctx.diff(radial,p[0])*ctx.diff(radial,p[1])
    measure=t*t*v*v*jac*(1-p[2])*(1-p[3])
    q=[ctx.sqrt(t*p[2])*ctx.cospi(2*p[4]),ctx.sqrt(t*p[2])*ctx.sinpi(2*p[4])]
    r=[ctx.sqrt(v*p[3])*ctx.cospi(2*p[5]),ctx.sqrt(v*p[3])*ctx.sinpi(2*p[5])]
    project=lambda q,r:[q[0],(-q[0]+ctx.sqrt(3)*q[1])/2,r[0],(-r[0]+ctx.sqrt(3)*r[1])/2]
    proj=project(q,r)
    external_norm=lambda rr:ctx.mpf(rr[2]**2+rr[3]**2-rr[2]*rr[3])
    def linear(rr,pr):return 2*(rr[0]*(rr[2]*pr[0]+rr[3]*pr[1])+rr[1]*(rr[2]*pr[2]+rr[3]*pr[3]))
    # The two cross edges are identified from original loop routes.
    other=[rr for rr in binding['host_routes'] if rr[0]*rr[1]==0]
    cross=[rr for rr in binding['host_routes'] if rr[0]*rr[1]!=0]
    assert len(other)==4 and cross==[[-1,-1,0,0],[1,1,1,1]]
    base=lambda rr:rr[0]**2*t+rr[1]**2*v
    baseline=[1+base(rr) for rr in other];d4=ctx.fprod(baseline)
    t0=[-H.edge_log(ctx,x)/2 for x in (t,v)];p0=ctx.fprod(t0)
    plane_dot=q[0]*r[0]+q[1]*r[1]
    aa=1+t+v+2*plane_dot;bb=2*ctx.sqrt(t*(1-p[2])*v*(1-p[3]))
    lc=linear(cross[1],proj)
    factors=[]
    for sign in (-1,1):
        norm=lambda rr:base(rr)+sign*linear(rr,proj)+external_norm(rr)
        ratio=ctx.fprod(d/(1+norm(rr)) for d,rr in zip(baseline,other))
        children=[]
        for rows in binding['child_pair_differences']:
            av,bv=[norm(rr) for rr in rows[:2]]
            if quadrature:
                children.append(-ctx.quad(lambda x:(1-x)*M.triangle_log_marginal(ctx,x,av,bv,ctx.mpf(1)),[0,ctx.mpf('.5'),1]))
            else:children.append(-H.triangle_log_complete(ctx,av,bv))
        factors.append((sign,ratio*ctx.fprod(children)))
    # This reference obtains M2..M4 by differentiation of M1, unlike native formulas.
    first=lambda a:2/(a+ctx.sqrt(a*a-bb*bb))
    moments=[(-1)**n*ctx.diff(first,aa,n)/ctx.factorial(n) for n in range(4)]
    m1,m2,m3,m4=moments
    value=sum(f*A.pair(ctx,aa,aa+sign*lc+1,bb)/2 for sign,f in factors)-p0*m2
    control=mode!='six-direct'
    if control:
        scalar=[]
        for x in (t,v):
            if quadrature:
                ell=lambda u:u*(1-u);br=[0,ctx.mpf('.5'),1]
                u1=ctx.quad(lambda u:-(1-u)*ell(u)/(1+x*ell(u)),br)
                cp=u1/2;cf=ctx.quad(lambda u:-(1-u)**3/(6*(1+x*ell(u))),br)
                i2=ctx.quad(lambda u:ell(u)**2/(1+x*ell(u))**2,br)
            else:
                e1=ctx.diff(lambda z:H.edge_log(ctx,z),x);e2=ctx.diff(lambda z:H.edge_log(ctx,z),x,2)
                u1=-e1/2;cp=-e1/4;cf=-(1-(x+3)*e1)/12;i2=-e2
            scalar.append((u1,cp,cf,i2))
        def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(3)]
        def coefficients(pr):
            jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
            for rr,d in zip(other,baseline):
                l,c=linear(rr,pr)/d,external_norm(rr)/d
                jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
            for i,rows in enumerate(binding['child_pair_differences']):
                u1,cp,cf,i2=scalar[i];a,b=linear(rows[0],pr)/2,linear(rows[1],pr)/2
                jet=mul(jet,[t0[i],u1*(a+b),cp*(external_norm(rows[0])+external_norm(rows[1]))+cf*external_norm(rows[2])+i2*(a*a+a*b+b*b)/3])
            return jet
        def average(h):
            f2=[];lf1=[]
            for i in range(6):
                for sign in (-1,1):
                    qq=[ctx.mpf(0)]*6;rr=[ctx.mpf(0)]*6
                    qq[i]=ctx.sqrt(t);rr[i]=h/ctx.sqrt(t);rr[(i+1)%6]=sign*ctx.sqrt(v-h*h/t)
                    pr=project(qq,rr);jet=coefficients(pr)
                    f2.append(jet[2]);lf1.append(linear(cross[1],pr)*jet[1])
            return sum(f2)/12,sum(lf1)/12
        h1=ctx.sqrt(t*v)/2
        af,al=average(ctx.mpf(0));f,l=average(h1);bf,bl=(f-af)/h1,(l-al)/h1
        # Check the affine property independently at the opposite Gram value.
        f,l=average(-h1)
        assert abs(f-(af-bf*h1))<=ctx.power(10,-digits+20)*max(1,abs(f))
        assert abs(l-(al-bl*h1))<=ctx.power(10,-digits+20)*max(1,abs(l))
        jpoint=coefficients(proj);f1,f2=jpoint[1:]
        point2=m2*f2-m3*(lc*f1+p0)+m4*lc*lc*p0
        mean2=af*m2+bf*(m1-(1+t+v)*m2)/2-al*m3-bl*(m2-(1+t+v)*m3)/2-p0*m3+ctx.mpf(2)/3*p0*(m3-m4)
        value+=mean2-point2
    if quadrature:
        def angular(xi):
            xx=aa+bb*xi;h=plane_dot+bb*xi/2
            result=sum(f/(xx*(xx+sign*lc+1))/2 for sign,f in factors)-p0/(xx*xx)
            if control:
                result+=(af+bf*h-f2)/(xx*xx)-(al+bl*h-lc*f1)/(xx**3)
                result+=(ctx.mpf(2)/3*(t+v+2*h)-lc*lc)*p0/(xx**4)
            return result
        value=A.angular_integral(ctx,angular,aa,bb,min(aa-lc+1,aa+lc+1))
    return ctx,value*measure/d4
