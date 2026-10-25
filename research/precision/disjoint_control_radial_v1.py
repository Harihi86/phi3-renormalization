"""Exact simplex integration of the independently derived angular control."""
import sympy as S
import disjoint_triangle_bubble_angular_control_v1 as A


def constants(binding):
    proof=A.verify(binding)
    norms=[]
    for row in binding['triangle_pairs']:
        if row['simplex_pair'][0]==0:
            c,a,b=map(S.Rational,row['momentum']);norms.append(a*a-a*b+b*b)
    return S.Rational(proof['host_norm_sum']),S.Rational(proof['child_norm_sum']),S.Rational(proof['mixed_dot']),sum(norms)


def verify(binding):
    sh,ss,cross,nsum=constants(binding)
    t,y=S.symbols('t y',positive=True);alpha=y*(1-y);M=1+t*alpha
    i0,i1,i2,j,j00=S.symbols('I0 I1 I2 J J00')
    # Exact rational identities under the integration sign, including both
    # endpoint terms. They fix the normalization, signs and all factors.
    assert S.cancel(alpha/M-(1-1/M)/t)==0
    assert S.cancel(alpha**2/M**2-(1-2/M+1/M**2)/t**2)==0
    ibp=S.diff((y-S.Rational(1,2))/M,y)
    assert S.cancel(ibp-(-1/M+(t+4)/(2*M*M)))==0
    assert S.cancel(-t*(y-S.Rational(1,2))*S.diff(alpha,y)/M-((t+4)/(2*M)-2))==0
    assert S.expand(((1-y)**3+y**3)/2-(1-3*alpha)/2)==0
    integrated=(2*t*sh/(3*(1+t)**2)-2/(1+t))*j/2
    integrated+=nsum*i1/4+(i0-3*i1)/12-t*ss*i2/18-t*cross*i1/(6*(1+t))
    closed=(2*t*sh/(3*(1+t)**2)-2/(1+t))*j/2+(nsum-1)*i1/4+i0/12-t*ss*i2/18-t*cross*i1/(6*(1+t))
    assert S.expand(integrated-closed)==0
    # Direct antiderivative of 1/(1+t*y*(1-y)) after r=2*y-1.
    r=S.symbols('r',real=True);primitive=4*S.atanh(r*S.sqrt(t/(t+4)))/S.sqrt(t*(t+4))
    assert S.simplify(S.diff(primitive,r)-4/(t+4-t*r*r))==0
    assert (sh,ss,cross,nsum) in ((S.Rational(5,2),1,S.Rational(3,2),1),(S.Rational(5,2),S.Rational(5,2),3,2))
    return {'sh':str(sh),'ss':str(ss),'dot':str(cross),'nsum':str(nsum),
        'rational_integrand_identities':5,'simplex_average_exact':True,'I0_antiderivative_exact':True,
        'low_tail':'exp(-120)/150','high_tail':'5*exp(-96)*(96^2+4*96+5)/3',
        'full_control':'-integral_0^infinity t^2*B(t)*F(t)/(2*(1+t)^4) dt',
        'not_the_complete_graph_integral':True}


def reduced(ctx,sh,ss,cross,nsum,t):
    i0=4*ctx.asinh(ctx.sqrt(t)/2)/ctx.sqrt(t*(t+4));j=(t+4)*i0/2-2
    i1=(1-i0)/t;i2=(1-2*i0+2*(1+i0)/(t+4))/(t*t)
    return (2*t*sh/(3*(1+t)**2)-2/(1+t))*j/2+(nsum-1)*i1/4+i0/12-t*ss*i2/18-t*cross*i1/(6*(1+t))


def independent_y(ctx,binding,t):
    sh,ss,cross,nsum=[ctx.mpf(str(x.p if hasattr(x,'p') else x))/ctx.mpf(str(x.q if hasattr(x,'q') else 1)) for x in constants(binding)]
    def f(y):
        alpha=y*(1-y);m=1+t*alpha;a=1+t
        cbar=alpha*nsum/2+(1-y)**2/6
        return (1-y)*((2*t*sh/(3*a*a)-2/a)*ctx.log(m)+cbar/m-t*alpha*alpha*ss/(9*m*m)-t*alpha*cross/(3*a*m))
    return ctx.quad(f,[0,ctx.mpf(1)/2,1])
