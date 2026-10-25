"""Complete angular reference with a 100-relative-digit explicit Cauchy bound."""
import mpmath as mp
import double_weighted_cauchy_v1 as M
def evaluate(points,binding=None,digits=640,chart='cauchy-cubic',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    assert binding['id']=='L4E2G20' and chart in ('cauchy-cubic','cauchy-diagonal')
    point=list(map(number,points));kind=chart.rsplit('-',1)[1]
    if not quadrature:
        value,bound,count=M.evaluate(ctx,point,kind,100)
        assert bound<ctx.mpf('1e-95')*abs(value)
        return ctx,value
    def full(s,t):
        result=ctx.quad(lambda f:ctx.mpf(8)/3*ctx.sinpi(f)**4*M.bare(ctx,s,t,ctx.cospi(f)),[0,ctx.mpf(1)/4,ctx.mpf(1)/2,ctx.mpf(3)/4,1])
        mass,kin=M.A.mass_kinetic(ctx,s,t)
        return result-mass-kin-M.A.parents(ctx,s,t)
    convert=lambda v:(v/(1-v)**3,(1+2*v)/(1-v)**4)
    a,da=convert(point[0]);b,db=convert(point[1]);bubble=M.A.G.B.bubble
    if kind=='cubic':return ctx,da*db*a*a*b*b*bubble(ctx,a)*bubble(ctx,b)/(4*(1+a)**2*(1+b)**2)*full(a,b)
    s,t=(a+b)**2,a*a;pref=4*a*(a+b)*da*db*s*s*t*t*bubble(ctx,s)*bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)
    return ctx,pref*(full(s,t)+full(t,s))
