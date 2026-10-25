"""Alternative independent scalar reference on a positive hypergeometric argument.

Pfaff identity: https://dlmf.nist.gov/15.8.E1 (same c, so normalization cancels).
This module is not substituted into any executed or running proof.
"""
import sympy as S


def exact_identity():
    z=S.symbols('z');g0,g1,g2=S.symbols('g0 g1 g2');x=z/(z-1);w=1/(1-z)
    yp=S.diff(w,z)*g0+w*S.diff(x,z)*g1
    ypp=S.diff(w,z,2)*g0+(2*S.diff(w,z)*S.diff(x,z)+w*S.diff(x,z,2))*g1+w*S.diff(x,z)**2*g2
    original=z*(1-z)*ypp+(S.Rational(7,2)-3*z)*yp-w*g0
    transformed=x*(1-x)*g2+(S.Rational(7,2)-S.Rational(9,2)*x)*g1-S.Rational(5,2)*g0
    multiplier=S.cancel(S.expand(original).coeff(g2)/(x*(1-x)))
    assert S.cancel(original-multiplier*transformed)==0
    assert w.subs(z,0)==1 and x.subs(z,0)==0
    assert yp.subs({z:0,g0:1,g1:S.Rational(5,7)})==S.Rational(2,7)
    t=S.symbols('t',positive=True)
    assert S.cancel(x.subs(z,-t/4)-t/(t+4))==0
    assert S.cancel(t*t*w.subs(z,-t/4)/60-t*t/(15*(t+4)))==0
    return {'ODE_multiplier':str(multiplier),'regular_origin_value':'1','regular_origin_derivative':'2/7',
        'domain':'t>=0; principal branch, 1-z>0 and 0<=x<1','ODE_and_initial_coefficients_exact':True}


def hyper_bubble(ctx,t):
    t=ctx.mpf(t);assert t>=0
    if not t:return ctx.mpf(0)
    x=t/(t+4)
    assert 0<x<1,'reference precision must resolve the positive transformed argument'
    return t*t/(15*(t+4))*ctx.hyp2f1(1,ctx.mpf(5)/2,ctx.mpf(7)/2,x)
