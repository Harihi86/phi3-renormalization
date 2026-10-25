"""Closed two-loop triangle-child radial marginal and independent chain-rule jets."""
import triangle_general_ad_reference_v1 as A

class Jet(A.Jet):
    def __init__(self,ctx,value,gradient=None,hessian=None):
        self.ctx=ctx;self.v=ctx.convert(value)
        self.g=[ctx.mpf(0)]*3 if gradient is None else [ctx.convert(x) for x in gradient]
        self.h=[[ctx.mpf(0)]*3 for _ in range(3)] if hessian is None else [[ctx.convert(x) for x in row] for row in hessian]
    def wrap(self,x):return x if isinstance(x,Jet) else Jet(self.ctx,x)
    def __add__(self,x):
        b=self.wrap(x);return Jet(self.ctx,self.v+b.v,[a+c for a,c in zip(self.g,b.g)],[[self.h[i][j]+b.h[i][j] for j in range(3)] for i in range(3)])
    __radd__=__add__
    def __neg__(self):return Jet(self.ctx,-self.v,[-x for x in self.g],[[-x for x in r] for r in self.h])
    def __sub__(self,x):return self+-self.wrap(x)
    def __rsub__(self,x):return self.wrap(x)+-self
    def __mul__(self,x):
        b=self.wrap(x);return Jet(self.ctx,self.v*b.v,[self.g[i]*b.v+self.v*b.g[i] for i in range(3)],
            [[self.h[i][j]*b.v+self.v*b.h[i][j]+self.g[i]*b.g[j]+self.g[j]*b.g[i] for j in range(3)] for i in range(3)])
    __rmul__=__mul__
    def unary(self,value,first,second):
        return Jet(self.ctx,value,[first*x for x in self.g],[[first*self.h[i][j]+second*self.g[i]*self.g[j] for j in range(3)] for i in range(3)])
    def reciprocal(self):return self.unary(1/self.v,-1/self.v**2,2/self.v**3)
    def __truediv__(self,x):return self*self.wrap(x).reciprocal()
    def __rtruediv__(self,x):return self.wrap(x)*self.reciprocal()
    def __pow__(self,power):
        assert isinstance(power,int)
        if power==0:return Jet(self.ctx,1)
        if power==1:return self
        return self.unary(self.v**power,power*self.v**(power-1),power*(power-1)*self.v**(power-2))

class Context(A.Context):
    def __init__(self,ctx):self.base=ctx;self.j=Jet(ctx,ctx.j);self.pi=Jet(ctx,ctx.pi)
    def wrap(self,x):return x if isinstance(x,Jet) else Jet(self.base,x)
    def mpf(self,x):return Jet(self.base,self.base.mpf(x))
    def mpc(self,x):return Jet(self.base,self.base.mpc(x))
    def log1p(self,x):
        x=self.wrap(x);return x.unary(self.base.log1p(x.v),1/(1+x.v),-1/(1+x.v)**2)

def master(ctx,t):
    base=ctx.base if isinstance(t,Jet) else ctx
    value=t.v if isinstance(t,Jet) else t;z=t-1;radius=abs(value-1)
    if radius<base.mpf('0.25'):
        n=3 if radius==0 else max(3,int(base.ceil((base.dps+30)/(-base.log10(radius))))+3)
        result=ctx.pi**2/6-ctx.mpf(5)/4;power=z
        for k in range(1,n+1):
            result+=ctx.mpf(2 if k%2 else -2)*power/(k*k*(k+1)*(k+2));power*=z
        return result
    return ctx.pi**2/6-ctx.polylog(2,1-t)+t*(2-3*t)*ctx.log(t)/(2*(1-t)**2)-t/(2*(1-t))

def closed(ctx,q,ell,D,H):
    # H=Ge/ell is a positive semidefinite residual quadratic form.
    if not isinstance(q,Jet) and q==D==H==0:return ctx.mpf(0)
    a=1+q;b=1+ell*a+D;c=ell*(1+H)
    disc=(1-ell*a+D)**2+4*ell*a*(D-H)
    root=ctx.sqrt(disc);plus=(b+root)/(2*a);minus=2*c/(b+root);la=ctx.log1p(q)
    return la*ctx.log(ell)+la*la/2+3*la/2-master(ctx,plus/ell)-master(ctx,minus/ell)+master(ctx,1/ell)+ctx.pi**2/6-ctx.mpf(5)/4

def invariant_jets(ctx,q,ell,D,H):
    ad=Context(ctx);args=[]
    for i,x in enumerate((q,D,H)):
        g=[ctx.mpf(0)]*3;g[i]=1;args.append(Jet(ctx,x,g))
    result=closed(ad,args[0],ad.wrap(ell),args[1],args[2])
    return result.v,result.g,result.h

def direct(ctx,q,ell,D,H,derivative=None):
    def density(r):
        u=1+r*ell;v=u**-3;A=1+r+q;d=r*(D+r*ell*H)/u;B=A+d
        if derivative is None:return ((1-v)*ctx.log1p(q/(1+r))-v*ctx.log1p(d/A))/r
        coeff=[ctx.mpf(1),r/u,r*r*ell/u]
        if len(derivative)==1:
            i=derivative[0]
            if i==0:return ((1-v)/A+v*d/(A*B))/r
            return -v*coeff[i]/(r*B)
        i,j=derivative
        if i==j==0:return (-(1-v)/A**2-v*d*(A+B)/(A*A*B*B))/r
        return v*coeff[i]*coeff[j]/(r*B*B)
    stops=[ctx.mpf(0),ctx.mpf(1),1/ell,ctx.inf]
    return ctx.quad(density,stops)
