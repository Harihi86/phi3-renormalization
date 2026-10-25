"""Symbolic identities behind the independent real antiderivative evaluator."""
import sympy as s
MUTATIONS=('primitive-sign','primitive-factor','curvature-sign','endpoint-factor','bubble-scale')


def verify(mutation=None):
    a,f0,v=s.symbols('a f0 v',positive=True);b=s.symbols('b',real=True)
    f=f0+b*v-a*v*v;D=s.sqrt(b*b+4*a*f0);f1=f0+b-a
    primitive=(v-b/(2*a))*s.log(f)-2*v+D/a*s.atanh((2*a*v-b)/D)
    if mutation=='primitive-sign':primitive-=2*D/a*s.atanh((2*a*v-b)/D)
    if mutation=='primitive-factor':primitive-=D/(2*a)*s.atanh((2*a*v-b)/D)
    assert s.factor(s.diff(primitive,v)-s.log(f))==0,'Wrong logarithmic primitive'
    fc=f.subs(v,s.Rational(1,2))
    inferred=4*fc-2*(f0+f1)
    if mutation=='curvature-sign':inferred=-inferred
    assert s.expand(inferred-a)==0
    assert s.expand(f1-f0+inferred-b)==0
    assert s.expand(D**2-(2*a*v-b)**2-4*a*f)==0
    assert s.expand((2*f0+b)**2-D**2-4*f0*f1)==0
    upper=(2*a-b)/D;lower=-b/D
    endpoint=D/(2*f0+b)
    if mutation=='endpoint-factor':endpoint*=2
    assert s.factor((upper-lower)/(1-upper*lower)-endpoint)==0
    t=s.symbols('t',positive=True);u=s.sqrt(t/(t+4))
    # Both functions vanish at t=0, proving their equality on t>0.
    assert s.simplify(s.diff(s.atanh(u),t)-s.diff(s.asinh(s.sqrt(t)/2),t))==0
    assert s.limit(s.atanh(u),t,0,dir='+')==0
    assert s.limit(s.asinh(s.sqrt(t)/2),t,0,dir='+')==0
    # For B(t)=int_0^1 [(1+x(1-x)t)log(1+x(1-x)t)-x(1-x)t] dx,
    # B''(t)=int [x(1-x)]^2/[1+t*x(1-x)] dx. Polynomial division gives
    # 1/(6t)-1/t^2+I0/t^2, with I0=4atanh(u)/sqrt(t(t+4)).
    B=(t+4)/3*(s.atanh(u)/u-1)-t/9
    if mutation=='bubble-scale':B*=2
    I0=4*s.atanh(u)/s.sqrt(t*(t+4))
    assert s.simplify(s.diff(B,t,2)-(1/(6*t)-1/t**2+I0/t**2))==0
    assert s.limit(B,t,0,dir='+')==0 and s.limit(s.diff(B,t),t,0,dir='+')==0
    x=s.symbols('x',real=True);w=x*(1-x)
    assert s.factor(w*w/(1+t*w)-(w/t-1/t**2+1/(t**2*(1+t*w))))==0
    # Primitive of the rational I0 integrand, with endpoints yielding I0.
    F=2/s.sqrt(t*(t+4))*s.atanh((2*x-1)*u)
    assert s.simplify(s.diff(F,x)-1/(1+t*w))==0
    assert s.simplify(F.subs(x,1)-F.subs(x,0)-I0)==0
    return {'quadratic_primitive_derivative_exact':True,'endpoint_identity_exact':True,
            'coefficients_reconstructed_from_three_original_values':True,
            'positive_branch_premises_checked_at_runtime':True,
            'bubble_integral_verified_by_second_derivative_and_two_initial_values':True,
            'atanh_asinh_identity_verified_by_derivative_and_initial_value':True,
            'normalized_host_or_moment_recurrence_used':False}
