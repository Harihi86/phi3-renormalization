"""Positive convergent-series reference with exact argument halving.

The analytic tail formula controls truncation only. mpmath rounding is not
certified here, so this is an independent reference, not an interval certificate.
No executed program imports this module until separately bound by a later proof.
"""
import sympy as S


def exact_identity():
    n=S.symbols('n',integer=True,nonnegative=True)
    assert S.simplify((S.rf(S.Rational(5,2),n)/S.rf(S.Rational(7,2),n)).rewrite(S.gamma)-5/(2*n+5))==0
    r=S.symbols('r',positive=True);double=2*r/(1+r*r)
    assert S.cancel(S.diff(double,r)/(1-double**2)-2/(1-r*r))==0
    a=S.symbols('a',positive=True);t=a*(a+4);q2=a/(a+4)
    assert S.cancel(4*q2/(1+q2)**2-t/(t+4))==0
    q,h=S.symbols('q h',positive=True);t=4*q*q/(1-q*q)
    first=t*t/(3*(t+4))*((h-q-q**3/3)/q**5)
    second=4*(h-q-q**3/3)/(3*q*(1-q*q))
    assert S.cancel(first-second)==0
    return {'positive_hypergeometric_coefficient':'5/(2n+5)',
        'bubble_series':'t^2/(3(t+4))*sum_{n>=0} x^n/(2n+5),x=t/(t+4)',
        'half_argument_identity':'t=a(a+4); atanh(q(t))=2atanh(q(a))',
        'tail':'x^N/((2N+c)(1-x)), c=5 for bubble series or c=1 for atanh/q',
        'tail_scope':'truncation in exact real arithmetic; not a floating point enclosure'}


def evaluate(ctx,t):
    t=ctx.mpf(t);assert t>=0
    if not t:return ctx.mpf(0),{'terms':0,'halvings':0,'model_absolute_tail':ctx.mpf(0)}
    with ctx.extradps(30):
        tolerance=ctx.power(10,-ctx.dps+15)
        def series(x,c):
            power=ctx.mpf(1);total=ctx.mpf(0)
            for n in range(100000):
                total+=power/(2*n+c);power*=x;N=n+1
                tail=power/((2*N+c)*(1-x))
                if tail<=tolerance*abs(total):return total,tail,N
            raise ArithmeticError('bounded convergent reference series exhausted')
        if t<=1:
            value,tail,terms=series(t/(t+4),5);prefactor=t*t/(3*(t+4))
            result=prefactor*value;bound=prefactor*tail;halvings=0
        else:
            reduced=t;halvings=0
            while reduced>ctx.mpf(4)/3:
                reduced=reduced/(2*(ctx.sqrt(1+reduced/4)+1));halvings+=1
                assert halvings<=64
            q=ctx.sqrt(reduced/(reduced+4));value,tail,terms=series(q*q,1)
            h=ctx.power(2,halvings)*q*value
            prefactor=(t+4)*ctx.sqrt(t+4)/(3*ctx.sqrt(t))
            result=prefactor*h-4*t/9-ctx.mpf(4)/3
            bound=prefactor*ctx.power(2,halvings)*q*tail
        assert result>0
        return result,{'terms':terms,'halvings':halvings,'model_absolute_tail':bound}


def hyper_bubble(ctx,t):return evaluate(ctx,t)[0]
