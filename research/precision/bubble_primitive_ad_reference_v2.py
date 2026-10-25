"""Scalar chain-rule jets; every constant enters the arbitrary-precision context."""
import triangle_general_ad_reference_v1 as A
import nested_analytic_marginals_v1 as M

class Jet:
    def __init__(self,ctx,value,first=0,second=0):
        self.ctx=ctx;self.v=ctx.convert(value);self.d=ctx.convert(first);self.dd=ctx.convert(second)
    def wrap(self,other):return other if isinstance(other,Jet) else Jet(self.ctx,other)
    def __add__(self,other):
        b=self.wrap(other);return Jet(self.ctx,self.v+b.v,self.d+b.d,self.dd+b.dd)
    __radd__=__add__
    def __neg__(self):return Jet(self.ctx,-self.v,-self.d,-self.dd)
    def __sub__(self,other):return self+-self.wrap(other)
    def __rsub__(self,other):return self.wrap(other)+-self
    def __mul__(self,other):
        b=self.wrap(other);return Jet(self.ctx,self.v*b.v,self.d*b.v+self.v*b.d,self.dd*b.v+2*self.d*b.d+self.v*b.dd)
    __rmul__=__mul__
    def unary(self,value,first,second):return Jet(self.ctx,value,first*self.d,first*self.dd+second*self.d*self.d)
    def reciprocal(self):return self.unary(1/self.v,-1/self.v**2,2/self.v**3)
    def __truediv__(self,other):return self*self.wrap(other).reciprocal()
    def __rtruediv__(self,other):return self.wrap(other)*self.reciprocal()
    def __pow__(self,power):
        assert isinstance(power,int)
        if power==0:return Jet(self.ctx,1)
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
    def real(self):return Jet(self.ctx,self.ctx.re(self.v),self.ctx.re(self.d),self.ctx.re(self.dd))
    @property
    def imag(self):return Jet(self.ctx,self.ctx.im(self.v),self.ctx.im(self.d),self.ctx.im(self.dd))

class Context(A.Context):
    def __init__(self,ctx):self.base=ctx;self.j=Jet(ctx,ctx.j)
    def wrap(self,value):return value if isinstance(value,Jet) else Jet(self.base,value)
    def mpf(self,value):return Jet(self.base,self.base.mpf(value))
    def mpc(self,value):return Jet(self.base,self.base.mpc(value))
    def log1p(self,value):
        x=self.wrap(value);return x.unary(self.base.log1p(x.v),1/(1+x.v),-1/(1+x.v)**2)
    def atan(self,value):
        x=self.wrap(value);d=1+x.v*x.v;return x.unary(self.base.atan(x.v),1/d,-2*x.v/d**2)

def variance_jets(ctx,x,variance):
    ad=Context(ctx);result=M.bubble_K(ad,ad.wrap(x),Jet(ctx,variance,1))
    return result.v,result.d,result.dd
