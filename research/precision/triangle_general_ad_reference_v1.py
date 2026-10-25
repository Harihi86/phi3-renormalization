"""Independent second-order automatic differentiation of the dilogarithm scalar.

This differentiates the scalar expression itself. It imports no simplex
boundary moments or native recurrence and performs no finite differences.
"""
import triangle_complete_general_v1 as Scalar


class Jet:
    def __init__(self,ctx,value,gradient=None,hessian=None):
        self.ctx=ctx;self.v=value
        self.g=[ctx.mpf(0)]*3 if gradient is None else list(gradient)
        self.h=[[ctx.mpf(0)]*3 for _ in range(3)] if hessian is None else [list(r) for r in hessian]
    def wrap(self,other):return other if isinstance(other,Jet) else Jet(self.ctx,other)
    def __add__(self,other):
        other=self.wrap(other);return Jet(self.ctx,self.v+other.v,[a+b for a,b in zip(self.g,other.g)],[[self.h[i][j]+other.h[i][j] for j in range(3)] for i in range(3)])
    __radd__=__add__
    def __neg__(self):return Jet(self.ctx,-self.v,[-x for x in self.g],[[-x for x in r] for r in self.h])
    def __sub__(self,other):return self+-self.wrap(other)
    def __rsub__(self,other):return self.wrap(other)+-self
    def __mul__(self,other):
        b=self.wrap(other);return Jet(self.ctx,self.v*b.v,[self.g[i]*b.v+self.v*b.g[i] for i in range(3)],
            [[self.h[i][j]*b.v+self.v*b.h[i][j]+self.g[i]*b.g[j]+self.g[j]*b.g[i] for j in range(3)] for i in range(3)])
    __rmul__=__mul__
    def unary(self,value,first,second):
        return Jet(self.ctx,value,[first*x for x in self.g],[[first*self.h[i][j]+second*self.g[i]*self.g[j] for j in range(3)] for i in range(3)])
    def reciprocal(self):return self.unary(1/self.v,-1/self.v**2,2/self.v**3)
    def __truediv__(self,other):return self*self.wrap(other).reciprocal()
    def __rtruediv__(self,other):return self.wrap(other)*self.reciprocal()
    def __pow__(self,power):
        assert isinstance(power,int)
        if power==0:return Jet(self.ctx,self.ctx.mpf(1))
        if power==1:return self
        return self.unary(self.v**power,power*self.v**(power-1),power*(power-1)*self.v**(power-2))
    def __bool__(self):return bool(self.v)
    def __eq__(self,other):return self.v==self.wrap(other).v
    def __ne__(self,other):return not self==other
    def __lt__(self,other):return self.v<self.wrap(other).v
    def __le__(self,other):return self.v<=self.wrap(other).v
    def __gt__(self,other):return self.v>self.wrap(other).v
    def __ge__(self,other):return self.v>=self.wrap(other).v
    @property
    def real(self):return Jet(self.ctx,self.ctx.re(self.v),[self.ctx.re(x) for x in self.g],[[self.ctx.re(x) for x in r] for r in self.h])
    @property
    def imag(self):return Jet(self.ctx,self.ctx.im(self.v),[self.ctx.im(x) for x in self.g],[[self.ctx.im(x) for x in r] for r in self.h])


class Context:
    def __init__(self,ctx):self.base=ctx;self.j=Jet(ctx,ctx.j)
    def wrap(self,value):return value if isinstance(value,Jet) else Jet(self.base,value)
    def mpf(self,value):return Jet(self.base,self.base.mpf(value))
    def mpc(self,value):return Jet(self.base,self.base.mpc(value))
    def sqrt(self,value):
        x=self.wrap(value);root=self.base.sqrt(x.v)
        return x.unary(root,1/(2*root),-1/(4*x.v*root))
    def log(self,value):
        x=self.wrap(value);return x.unary(self.base.log(x.v),1/x.v,-1/x.v**2)
    def asinh(self,value):
        x=self.wrap(value);root=self.base.sqrt(1+x.v*x.v)
        return x.unary(self.base.asinh(x.v),1/root,-x.v/root**3)
    def polylog(self,degree,value):
        assert degree==2;x=self.wrap(value)
        if x.v==0:return x.unary(self.base.mpf(0),self.base.mpf(1),self.base.mpf('.5'))
        log=self.base.log(1-x.v)
        return x.unary(self.base.polylog(2,x.v),-log/x.v,1/(x.v*(1-x.v))+log/x.v**2)
    def im(self,value):return self.wrap(value).imag
    def re(self,value):return self.wrap(value).real


def invariant_jets(ctx,a,b,c):
    ad=Context(ctx);args=[]
    for i,x in enumerate((a,b,c)):
        g=[ctx.mpf(0)]*3;g[i]=1;args.append(Jet(ctx,x,g))
    result=Scalar.triangle_log(ad,*args)
    return result.v,result.g,result.h
