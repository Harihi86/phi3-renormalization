"""Hermite reduction of B(t)^2/(1+t)^3 to rational functions and Li_1, Li_2, Li_3."""
from fractions import Fraction as Q
import sympy as S
from sympy.integrals.rationaltools import ratint_ratpart
w=S.symbols('w',positive=True)
L=S.symbols('L',real=True)


def rational(value):
    value=S.Rational(value);return Q(int(value.p),int(value.q))


def rational_record(value):
    num,den=S.fraction(S.cancel(value))
    return {'numerator':[rational(v) for v in S.Poly(num,w).all_coeffs()],
            'denominator':[rational(v) for v in S.Poly(den,w).all_coeffs()]}


def hermite(value):
    value=S.cancel(value);num,den=S.fraction(value)
    quotient,remainder=S.div(num,den,w)
    polynomial=S.Poly(quotient,w).integrate().as_expr()
    if remainder:
        a,b=ratint_ratpart(S.Poly(remainder,w),S.Poly(den,w),w)
    else:a,b=S.Integer(0),S.Integer(0)
    a=S.cancel(a+polynomial);b=S.cancel(b)
    assert S.cancel(S.diff(a,w)+b-value)==0
    return a,b


def residues(value):
    value=S.cancel(value);num,den=S.fraction(value)
    def residue(pole):
        if S.simplify(den.subs(w,pole))!=0:return S.Integer(0)
        assert S.simplify(S.diff(den,w).subs(w,pole))!=0
        return S.simplify(num.subs(w,pole)/S.diff(den,w).subs(w,pole))
    real={str(p):rational(residue(S.Integer(p))) for p in (0,1,-1)}
    rho=(1+S.I*S.sqrt(3))/2;coefficient=residue(rho)
    a=S.simplify(S.re(coefficient));b=S.simplify(S.im(coefficient)/S.sqrt(3))
    assert a.is_Rational and b.is_Rational
    assert S.simplify(residue(S.conjugate(rho))-S.conjugate(coefficient))==0
    rebuilt=sum(S.Rational(v.numerator,v.denominator)/(w-int(k)) for k,v in real.items())
    rebuilt+=(2*a*w-a-3*b)/(w*w-w+1)
    assert S.cancel(value-rebuilt)==0
    return {'real_poles':real,'rho_plus':[rational(a),rational(b)]}


def derive(moment=0):
    assert moment in (0,1)
    t=(1-w)**2/w;dt=S.diff(t,w)
    b1=-(1+w)**3/(6*w*(1-w));b0=-4*(w*w+w+1)/(9*w)
    density=S.Poly(S.cancel(t**moment*(b1*L+b0)**2/(1+t)**3*dt),L)
    original={k:S.cancel(density.coeff_monomial(L**k)) for k in (0,1,2)}
    rest=dict(original);parts={};simple={}
    for k in (2,1,0):
        parts[k],simple[k]=hermite(rest[k])
        if k:rest[k-1]=S.cancel(rest[k-1]-k*parts[k]/w)
    reconstructed=S.expand(sum((S.diff(parts[k],w)*L**k+(k*parts[k]*L**(k-1)/w if k else 0)+simple[k]*L**k) for k in (0,1,2)))
    for k in (0,1,2):assert S.cancel(S.Poly(reconstructed,L).coeff_monomial(L**k)-original[k])==0
    data={'moment':moment,'rational_parts':{str(k):rational_record(parts[k]) for k in (0,1,2)},
          'simple_poles':{str(k):residues(simple[k]) for k in (0,1,2)}}
    # Only the L^2 simple pole may occur at w=1; its primitive has a finite limit.
    assert data['simple_poles']['0']['real_poles']['1']==data['simple_poles']['1']['real_poles']['1']==0
    assert all(S.limit(parts[k],w,1).is_finite for k in (0,1,2))
    data['derivative_reconstruction_exact']=True
    data['algebra']={'t_of_w':str(t),'dt_dw':str(dt),'bubble_log_coefficient':str(b1),'bubble_constant':str(b0),
                    'rational_parts':{str(k):str(parts[k]) for k in (0,1,2)},
                    'simple_poles':{str(k):str(simple[k]) for k in (0,1,2)}}
    return data


def exact_primitive_rules():
    a=S.symbols('a',nonzero=True);log=S.log(w);ell=S.log(1-w/a)
    phi=[ell,log*ell+S.polylog(2,w/a),log**2*ell+2*log*S.polylog(2,w/a)-2*S.polylog(3,w/a)]
    for k in range(3):
        derivative=S.diff(phi[k],w).subs(S.polylog(1,w/a),-S.log(1-w/a))
        assert S.simplify(derivative-log**k/(w-a))==0
        assert S.diff(log**(k+1)/(k+1),w)==log**k/w
    q=S.symbols('q',positive=True)
    bubble_q=4*(3*(-L/2)-3*q-q**3)/(9*q*(1-q*q))
    bubble_w=-(1+w)**3*L/(6*w*(1-w))-4*(w*w+w+1)/(9*w)
    assert S.cancel(bubble_q.subs(q,(1-w)/(1+w))-bubble_w)==0
    assert S.cancel(4*q*q/(1-q*q)-(1-(1-q)/(1+q))**2/((1-q)/(1+q)))==0
    assert S.limit(bubble_w.subs(L,S.log(w))/(w-1)**4,w,1)==S.Rational(1,60)
    return {'nonzero_pole_rules':3,'zero_pole_rules':3,'rationalizing_substitution_exact':True,'bubble_small_t_coefficient':Q(1,60)}


def number(ctx,q):return ctx.mpf(q.numerator)/q.denominator


def evaluate_w(ctx,data,value):
    value=ctx.mpf(value);assert 0<value<=1;log=ctx.log(value)
    def poly(coefficients):
        result=ctx.mpf(0)
        for coefficient in coefficients:result=result*value+number(ctx,coefficient)
        return result
    total=ctx.mpf(0)
    for k in (0,1,2):
        part=data['rational_parts'][str(k)]
        total+=poly(part['numerator'])/poly(part['denominator'])*log**k
    def phi(pole,k):
        if value==1 and pole==1:
            assert k==2;return -2*ctx.zeta(3)
        z=value/pole;ell=ctx.log(1-z)
        if k==0:return ell
        li2=ctx.polylog(2,z)
        if k==1:return log*ell+li2
        return log*log*ell+2*log*li2-2*ctx.polylog(3,z)
    for k in (0,1,2):
        row=data['simple_poles'][str(k)]
        for pole,coefficient in row['real_poles'].items():
            if not coefficient:continue
            primitive=log**(k+1)/(k+1) if pole=='0' else phi(int(pole),k)
            total+=number(ctx,coefficient)*primitive
        a,b=row['rho_plus'];coefficient=number(ctx,a)+ctx.j*ctx.sqrt(3)*number(ctx,b)
        if coefficient:total+=2*ctx.re(coefficient*phi((1+ctx.j*ctx.sqrt(3))/2,k))
    assert not ctx.im(total)
    return ctx.re(total)


def primitive(ctx,data,t):
    t=ctx.mpf(t);assert t>=0
    if not t:return ctx.mpf(0)
    # Avoid loss in 1-sqrt(t/(t+4)) at large t.
    value=4/(ctx.sqrt(t+4)+ctx.sqrt(t))**2
    return evaluate_w(ctx,data,value)-evaluate_w(ctx,data,1)
