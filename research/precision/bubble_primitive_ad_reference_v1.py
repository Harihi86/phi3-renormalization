"""Second-order chain-rule differentiation of the scalar dilogarithm primitive."""
import triangle_general_ad_reference_v1 as A
import nested_analytic_marginals_v1 as M

class Context(A.Context):
    def log1p(self,value):
        x=self.wrap(value);return x.unary(self.base.log1p(x.v),1/(1+x.v),-1/(1+x.v)**2)
    def atan(self,value):
        x=self.wrap(value);d=1+x.v*x.v;return x.unary(self.base.atan(x.v),1/d,-2*x.v/d**2)

def variance_jets(ctx,x,variance):
    ad=Context(ctx);z=A.Jet(ctx,variance,[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)])
    result=M.bubble_K(ad,ad.wrap(x),z)
    return result.v,result.g[0],result.h[0][0]
