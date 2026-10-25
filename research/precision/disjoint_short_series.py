"""Short complex series with explicit tails for proving existing geometry bounds.

No tail is dropped. If the short enclosure is insufficient, callers retain the
frozen evaluator. 384-bit interval arithmetic and its integer cap are unchanged.
"""
from fractions import Fraction as Q
import disjoint_complex_enclosures as F
I,Z,A=F.I,F.Z,F.A


def choose_order(radius):
    radius=A.bounded(Q(radius));assert 0<=radius<Q(3,4)
    for bound,order in ((Q(1,64),16),(Q(1,16),24),(Q(1,4),48),(Q(1,2),96)):
        if radius<=bound:return order
    return 128


def tail(radius,order,kind):
    r=A.bounded(Q(radius));assert 0<=r<1 and type(order) is int and 16<=order<=128
    if kind=='log1p':return A.bounded(r**(order+1)/((order+1)*(1-r)))
    assert kind=='atanhc'
    return A.bounded(r**order/((2*order+1)*(1-r)))


def log1p(z,trace=None):
    z=F.G.as_z(z);r=z.norm().hi
    if r>Q(1,2):
        if trace is not None:trace['log_fallbacks']=trace.get('log_fallbacks',0)+1
        return F.log_right(1+z)
    order=choose_order(r)
    # Horner evaluation of z - z^2/2 + ... + (-1)^(N+1) z^N/N.
    result=F.G.as_z(Q((-1)**(order+1),order))
    for k in range(order-1,0,-1):result=result*z+Q((-1)**(k+1),k)
    result=result*z;error=I(-tail(r,order,'log1p'),tail(r,order,'log1p'))
    if trace is not None:
        trace['log_series_calls']=trace.get('log_series_calls',0)+1
        trace['log_series_terms']=trace.get('log_series_terms',0)+order
    return Z(result.re+error,result.im+error)


def atanhc(z,trace=None):
    z=F.G.as_z(z);r=z.norm().hi;assert r<Q(3,4),'atanhc analytic disk premise'
    order=choose_order(r);value=F.atanhc(z,order=order)
    if trace is not None:
        trace['atanhc_series_calls']=trace.get('atanhc_series_calls',0)+1
        trace['atanhc_series_terms']=trace.get('atanhc_series_terms',0)+order
    return value


def verify_identities(mutation=None):
    import sympy as s
    z,x=s.symbols('z x')
    for n in (2,3,16):
        logarithm=sum(s.Rational((-1)**(k+1),k)*z**k for k in range(1,n+1))
        if mutation=='log-sign':logarithm=-logarithm
        assert s.factor(s.diff(logarithm,z)-(1-(-z)**n)/(1+z))==0
        truncated=sum((z*x*x)**j for j in range(n))
        remainder=(z*x*x)**n/(1-z*x*x)
        if mutation=='atanhc-remainder-sign':remainder=-remainder
        assert s.factor(1/(1-z*x*x)-truncated-remainder)==0
        assert s.integrate(truncated,(x,0,1))==sum(z**j/s.Integer(2*j+1) for j in range(n))
    # At the geometric majorants' endpoint, the first omitted monomial must
    # already be covered. A missing tail is an invalid enclosure.
    for r in (Q(1,128),Q(1,32),Q(1,8),Q(3,8),Q(2,3)):
        n=choose_order(r)
        for kind,first in (('log1p',r**(n+1)/(n+1)),('atanhc',r**n/(2*n+1))):
            bound=tail(r,n,kind)
            if mutation=='zero-tail':bound=Q()
            assert bound>=first>0
    return {'log_derivative_geometric_identity_exact':True,'atanhc_integrand_geometric_identity_exact':True,
            'straight_line_log_remainder':'r^(N+1)/((N+1)*(1-r))',
            'atanhc_integral_remainder':'r^N/((2N+1)*(1-r))','all_tails_retained':True,
            'no_geometry_allowance_changes':True}
