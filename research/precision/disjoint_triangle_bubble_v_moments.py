"""Exact integration of the last simplex coordinate in a quadratic logarithm.

Let f(v)=fc*(1+beta*s-gamma*s^2), s=2*v-1, gamma>=0 and f(0),f(1)>0.
Factor f/fc=(1-lambda_1*s)(1-lambda_2*s). Both real roots lambda_i lie
strictly between -1 and 1. With r_n=lambda_1^(2n)+lambda_2^(2n),

 integral_0^1 log(f(v)) dv = log(fc)-sum_{n>=1} r_n/(2*n*(2*n+1)),
 r_0=2, r_1=beta^2+2*gamma,
 r_n=(beta^2+2*gamma)*r_(n-1)-gamma^2*r_(n-2).

For rho=max|lambda_i|<1, the omitted positive sum after n=N is at most
  2*rho^(2N+2)/((2N+2)*(2N+3)*(1-rho^2)).
No sample-based quadrature error estimate is used for this coordinate.
"""
from fractions import Fraction as Q
import sympy as s


def verify(binding,mutation=None):
    t,y,v,q1,q2=s.symbols('t y v q1 q2')
    xs=[y,(1-y)*v,(1-y)*(1-v)];poly=1;cross={};inner=None
    for row in binding['triangle_pairs']:
        i,j=row['simplex_pair'];c,a,b=map(s.Rational,row['momentum'])
        norm=a*a-a*b+b*b
        poly+=xs[i]*xs[j]*(c*c*t+norm+2*c*(a*q1+b*q2))
        if i==0:
            assert c*c==1;cross[j]=(a*c,b*c,norm)
        else:
            assert (i,j)==(1,2) and c==0 and norm==1;inner=norm
    assert set(cross)=={1,2} and inner==1
    da,db=cross[1][0]-cross[2][0],cross[1][1]-cross[2][1]
    assert da*da-da*db+db*db==1,'The two outer shifts differ by a unit external leg'
    fc=s.expand(poly.subs(v,s.Rational(1,2)));slope=s.expand(s.diff(poly,v).subs(v,s.Rational(1,2)))
    curvature=(1-y)**2
    if mutation=='curvature-sign':curvature=-curvature
    beta=slope/(2*fc);gamma=curvature/(4*fc);z=s.symbols('z')
    assert s.cancel(poly.subs(v,(1+z)/2)-fc*(1+beta*z-gamma*z*z))==0
    p,r=s.symbols('lambda1 lambda2');B=-(p+r);G=-p*r
    assert s.expand((1-p*z)*(1-r*z)-(1+B*z-G*z*z))==0
    coefficient=B*B+2*G
    if mutation=='recurrence-sign':coefficient=B*B-2*G
    assert s.expand((p*p+r*r)-coefficient)==0
    assert s.expand(p*p*r*r-G*G)==0
    # Characteristic identity proves the recurrence for every n, rather than
    # fitting it to a few numerical moments.
    assert s.expand(p**4-coefficient*p*p+G*G)==0
    assert s.expand(r**4-coefficient*r*r+G*G)==0
    B,G=s.symbols('B G');moments=[s.Integer(2),B*B+2*G]
    for n in range(2,9):moments.append(s.expand((B*B+2*G)*moments[-1]-G*G*moments[-2]))
    # An independent formal log expansion checks the exact integration factor.
    eps=s.symbols('eps');truncated=sum((-1)**(j+1)*(B*eps*z-G*eps**2*z*z)**j/s.Integer(j) for j in range(1,17))
    for n in range(1,9):
        coef=s.expand(truncated).coeff(eps,2*n)
        actual=s.integrate(coef,(z,-1,1))/2
        denominator=2*n*(2*n+1)
        if mutation=='simplex-measure':denominator*=2
        assert s.expand(actual+moments[n]/denominator)==0
    # Unit external separation gives |beta|<=1/4, gamma<=1/4. Write
    # fc=D+alpha*w^2, D=1+(1-y)/4; Cauchy bounds the slope squared.
    a,w,D=s.symbols('alpha w2 D',nonnegative=True)
    assert s.expand((D+a*w)**2-4*D*a*w-(D-a*w)**2)==0
    assert s.expand(s.Rational(1,4)-y*(1-y)-(y-s.Rational(1,2))**2)==0
    assert s.Rational(2,3)**2-s.Rational(1,4)*s.Rational(2,3)-s.Rational(1,4)==s.Rational(1,36)
    return {'quadratic_from_graph_exact':True,'unit_external_separation_exact':True,
            'power_sum_recurrence_all_orders':True,'formal_moments_checked':8,
            'simplex_integral_measure':'dv on [0,1]','uniform_rho_upper':'2/3',
            'series_tail_sign':'omitted positive sum; subtract from log(fc)',
            'outer_integral_not_certified':True}


def uniform_tail(N):
    assert type(N) is int and 1<=N<=128
    return 2*Q(2,3)**(2*N+2)/((2*N+2)*(2*N+3)*(1-Q(2,3)**2))
