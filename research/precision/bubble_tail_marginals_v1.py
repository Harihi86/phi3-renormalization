"""Closed bubble tails and an exact parent-parameter mean; no inner sampling."""
from math import comb
import nested_analytic_marginals_v1 as M

# Each coefficient is real + i*sqrt(3)*imag, with exact rational strings.
COEFFICIENTS={
 3:(('0','-2/3'),('-2/3','0'),('0','4/27')),
 4:(('0','-14/27'),('-14/27','0'),('0','4/27'),('8/81','0')),
 5:(('0','-4/9'),('-4/9','0'),('0','32/243'),('8/81','0'),('0','-16/729')),
}


def bubble(ctx,t):
 if t==0:return ctx.mpf(0)
 return (t+4)**ctx.mpf('1.5')/(3*ctx.sqrt(t))*ctx.asinh(ctx.sqrt(t)/2)-4*t/9-ctx.mpf(4)/3


def pole_integrals(ctx,r,z):
 """I_j=int_r^1 atanh(v)/(v-z)^j dv; log cut is never crossed."""
 def primitive(sign,v):
  return ctx.log(1+sign*z)*ctx.log(v-z)-ctx.polylog(2,-sign*(v-z)/(1+sign*z))
 result=[ctx.mpf(0),sum(sign*(primitive(sign,1)-primitive(sign,r)) for sign in (-1,1))/2]
 ell={1:ctx.log(1-z)-ctx.log(r-z)}
 J={0:ctx.log(2/(1+r))}
 for k in range(1,5):
  if k>1:ell[k]=((1-z)**(1-k)-(r-z)**(1-k))/(1-k)
  J[k]=(ell[k]-J[k-1])/(1+z)
  v=((r-z)**(-k)-(1-z)**(-k))/(-k)
  result.append(-v*ctx.atanh(r)+sum((1-z)**(m-k-1)*J[m] for m in range(1,k+1))/k)
 return result


def rational_power_integral(ctx,r,z,j,power):
 result=ctx.mpc(0)
 for l in range(power+1):
  exponent=l-j+1
  primitive=(ctx.log(1-z)-ctx.log(r-z)) if exponent==0 else ((1-z)**exponent-(r-z)**exponent)/exponent
  result+=comb(power,l)*z**(power-l)*primitive
 return result


def tails(ctx,a):
 """T_n(a)=int_a^infinity B(u)/(1+u)^n du for n=3,4,5."""
 if not a>=0:raise ValueError('nonnegative bubble-tail lower limit required')
 r=ctx.sqrt(a/(a+4));z=ctx.j/ctx.sqrt(3)
 integrals=pole_integrals(ctx,r,z);pieces={j:integrals[j]-rational_power_integral(ctx,r,z,j,1)-rational_power_integral(ctx,r,z,j,3)/3 for j in range(1,6)}
 return tuple(2*ctx.re(sum((ctx.mpf(re)+ctx.j*ctx.sqrt(3)*ctx.mpf(im))*pieces[j] for j,(re,im) in enumerate(COEFFICIENTS[n],1))) for n in (3,4,5))


def parent_moments(ctx,xi,t):
 """Integrals of K0,x*K',x*(1-x)*K',x^2*K'' over 0<x<1-xi."""
 if not 0<xi<1 or not t>=0:raise ValueError('strict parent variance coordinate and nonnegative radius required')
 X=1-xi;d=1+xi*X*t;a=d/X-1;b=1/X-1
 T3,T4,T5=tails(ctx,a);T30=tails(ctx,b)[0]
 I0=-X*(M.bubble_primitive(ctx,a)-M.bubble_primitive(ctx,b))-d*T3+T30
 I1=-d*T4;Ic=I1+d*d*T5;I2=d*bubble(ctx,a)/(1+a)**4-2*d*T5
 return I0,I1,Ic,I2
