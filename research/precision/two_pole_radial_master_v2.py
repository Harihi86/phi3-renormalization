"""High-precision partial-fraction radial master, for proof/reference only."""
from math import factorial


def phi(ctx,t,n):
    z=t-1;constant=ctx.pi**2/6-sum(ctx.mpf(1)/(j*j) for j in range(1,n))
    if not z:return constant
    if abs(z)<ctx.mpf('.25'):
        count=int(ctx.ceil((ctx.dps+12)*ctx.log(10)/(-ctx.log(abs(z)))))+2;power=z;result=constant
        for j in range(1,count+1):
            coefficient=ctx.mpf(factorial(n-1)*factorial(j-1))/(j*factorial(n+j-1))
            result+=(-1)**(j-1)*coefficient*power;power*=z
        return result
    value=ctx.pi**2/6-ctx.polylog(2,1-t);moment=ctx.log(t)/(t-1)
    for j in range(1,n):
        value-=t*moment/j;moment=(ctx.mpf(1)/j-t*moment)/(1-t)
    return value


def primitive(ctx,q,denominator,numerator,lambdas):
    # D starts at one; numerator is Q*D. Both use ascending powers of s.
    d=list(denominator);f=list(numerator)
    while len(d)>1 and d[-1]==0:d.pop()
    P=[ctx.mpf(0)]*max(len(d)+1,len(f))
    for j,x in enumerate(d):P[j]+=x;P[j+1]+=x
    for j,x in enumerate(f):P[j]+=x
    while P[-1]==0:P.pop()
    assert len(lambdas)==len(P)-1 and all(abs(ctx.im(x))<ctx.mpf('1e-80') and ctx.re(x)>0 for x in lambdas)
    lambdas=list(map(ctx.re,lambdas));a=P[0];assert abs(a-1-q)<ctx.mpf('1e-80')
    if len(d)==2:
        poles=[d[1]];C=ctx.log(poles[0])+ctx.mpf(3)/2;master=lambda c:phi(ctx,c/poles[0],3)
    else:
        assert len(d)==3;disc=ctx.sqrt(d[1]**2-4*d[2]);hi=(d[1]+disc)/2;lo=2*d[2]/(d[1]+disc);poles=[lo,hi]
        if abs(hi-lo)<ctx.power(10,-ctx.dps/2):
            h=(hi+lo)/2;C=ctx.log(h)+sum(ctx.mpf(1)/j for j in range(1,6));master=lambda c:phi(ctx,c/h,6)
        else:
            partial=[]
            for h,other in ((lo,hi),(hi,lo)):
                partial.extend([(h,3,h**3/(h-other)**3),(h,2,-3*h**3*other/(h-other)**4),(h,1,6*h**3*other**2/(h-other)**5)])
            assert abs(sum(co for h,n,co in partial)-1)<ctx.mpf('1e-80')
            C=sum(co*(ctx.log(h)+sum(ctx.mpf(1)/j for j in range(1,n))) for h,n,co in partial)
            master=lambda c:sum(co*phi(ctx,c/h,n) for h,n,co in partial)
    la=ctx.log(a)
    return la*la/2+C*la-sum(master(lam) for lam in lambdas)+master(ctx.mpf(1))+sum(master(h) for h in poles)
