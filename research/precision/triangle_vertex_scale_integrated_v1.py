"""Analytic relative-scale integral for L3E3G6; seven remaining variables.

The closed form is exact. Floating values returned here are diagnostics;
they are not interval enclosures or whole-graph certificates.
"""
from fractions import Fraction as Q
import sympy as S
import triangle_vertex_rank_one_v1 as F


def exact_identities():
    v, x = S.symbols('v x', positive=True)
    a = 1-v
    decomposition = -x/a+(2*a-1)/a**2+v*v/a**2/(1-a*x)
    assert S.cancel((1-x)**2/(1-a*x)-decomposition) == 0
    primitive = -x*x/(2*a)+(2*a-1)*x/a**2-v*v*S.log(1-a*x)/a**3
    assert S.cancel(S.diff(primitive,x)-(1-x)**2/(1-a*x)) == 0
    derivative = (1-3*v)/(2*(1-v)**2)-v*v*S.log(v)/(1-v)**3
    assert S.cancel(primitive.subs(x,1)-primitive.subs(x,0)-derivative) == 0
    closed = S.pi**2/6-S.polylog(2,1-v)-v/(2*(1-v))+v*(2-3*v)*S.log(v)/(2*(1-v)**2)
    derived = S.diff(closed,v).subs(S.polylog(1,1-v),-S.log(v))
    assert S.cancel(derived-derivative) == 0
    assert S.limit(closed,v,0,dir='+') == 0
    assert S.limit(closed,v,1) == S.pi**2/6-S.Rational(5,4)
    n = S.symbols('n',integer=True,positive=True)
    assert S.cancel(1/n*(1/n-2/(n+1)+1/(n+2))-2/(n*n*(n+1)*(n+2))) == 0
    # The roots of the bare quadratic are real and positive, including
    # coalescence. Their product is used for the smaller root numerically.
    aa, u, nn, d, sigma = S.symbols('A u N d sigma',positive=True)
    bb = 1+d+aa*u+nn
    cc = u*(1+d)
    assert S.expand(bb*bb-4*aa*cc-((1+d-aa*u)**2+2*nn*(1+d+aa*u)+nn*nn)) == 0
    q = aa-1
    relative = 1+sigma*u
    bare_numerator = (1+sigma)*relative+relative*q+sigma*(nn+relative*d)
    assert S.expand(bare_numerator-(aa+bb*sigma+cc*sigma*sigma)) == 0
    return {'K_derivative_proved_by_rational_partial_fractions':True,
        'K_zero_boundary_exact':True,'K_at_one':'pi^2/6-5/4',
        'K_local_series_coefficient':'2*(-1)^(n+1)/(n^2*(n+1)*(n+2))',
        'bare_quadratic_discriminant_nonnegative':True,
        'bare_original_rank_one_numerator_exact':True,
        'whole_graph_certificate':False}


def K(v, mp):
    assert v > 0
    constant = mp.pi**2/6-mp.mpf(5)/4
    if v == 1:
        return constant
    delta = v-1
    if abs(delta) <= mp.mpf('0.5'):
        result = constant
        power = delta
        for n in range(1,10000):
            result += (-1)**(n+1)*2*power/(n*n*(n+1)*(n+2))
            remainder = 2*abs(power*delta)/((n+1)**2*(n+2)*(n+3)*(1-abs(delta)))
            if remainder < mp.eps*max(abs(result),mp.mpf('1e-100')):
                return result
            power *= delta
        raise ArithmeticError('Finite diagnostic series cap exhausted')
    return mp.pi**2/6-mp.polylog(2,1-v)-v/(2*(1-v))+v*(2-3*v)*mp.log(v)/(2*(1-v)**2)


def scale_integral(q, u, n, d, mp):
    assert q >= 0 and u > 0 and n >= 0 and d >= 0
    aa = 1+q
    x, h = (1+d)/(aa*u), n/(aa*u)
    discriminant = (x-1)**2+2*h*(1+x)+h*h
    plus = (1+x+h+mp.sqrt(discriminant))/2
    minus = x/plus
    la = mp.log1p(q)
    return la*la/2+(mp.log(u)+mp.mpf('1.5'))*la-K(plus,mp)-K(minus,mp)+K(1/u,mp)+K(mp.mpf(1),mp)


def parameters(point, permutation):
    assert len(point) == 7
    chart = F.chart(tuple(point)+(Q(1,2),),permutation)
    s,v,f0,f1,f2,r0,t = point
    a0,a1,a2 = chart['host_weights']
    psi = chart['host_U']
    r = chart['triangle_weights'][0]
    a = r*(1-r)
    mu0 = a1*a2*(f0-f2)/psi-f0
    mu1 = a1*a2*(f1-f2)/psi
    u = a*(a1+a2)/psi
    n = a*F.norm(mu0+t,mu1)
    d = (1-r)*t*(1-t)
    rd = r0*r0+(1-r0)**2
    jacobian = v/(1+v+s*v)**3*(1-r)*2*r0*(1-r0)/rd**2
    # These bounds provide an integrable logarithmic majorant for the
    # absolute sigma integral. No cube region or endpoint is excluded.
    assert 0 < chart['host_Q'] <= 1 and 0 <= d <= Q(1,4)
    assert n/u <= 3 and r*(1-r) <= u <= Q(3,4)/(s*v)
    return {'Q':chart['host_Q'],'u':u,'N':n,'d':d,'jacobian':jacobian}


def evaluate(parameters, mp):
    def number(q):
        return mp.mpf(q.numerator)/q.denominator
    q,u,n,d = [number(parameters[k]) for k in ('Q','u','N','d')]
    return number(parameters['jacobian'])*scale_integral(q,u,n,d,mp)


def scale_density_z(z0, q, u, n, d, mp):
    """Independent one-dimensional reference using the proved finite log pair."""
    if z0 == 0:
        return mp.mpf(0)
    z,y,cc = z0*z0,1-z0*z0,1+u
    d1,d2,d3 = u+z,1+u*z,1+z*(u+q*cc)
    x0,x1 = q*cc*z/d2,y*(n*cc*z/d1+d)/d3
    ca = u*(d1*d1+d1*cc*z+(cc*z)**2)/d1**3*q*cc/d2
    cb = cc**3*z*z/d1**3*(n*cc*z/d1+d)/d3
    l = lambda x: mp.log1p(x)/x if x else mp.mpf(1)
    return 2*z0*(ca*l(x0)-cb*l(x1))
