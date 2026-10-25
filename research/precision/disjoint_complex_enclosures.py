"""Independent directed rational complex functions for mixed-graph replay.

No native library or floating-point result participates in an enclosure.
Arithmetic reuses the frozen 384-bit outward-rounded rectangle type and the
1,000,000-bit integer guard. Principal logarithms require the entire right
half-plane premise. Exponential phase winding is retained by periodic sine
and cosine ranges, with exact rational critical-point locations in units of pi.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
import theta_two_radial_geometry_audit_v4 as G
I,Z,A=G.I,G.Z,G.A


def square(z):
    z=G.as_z(z)
    return Z(z.re.square()-z.im.square(),2*z.re*z.im)


def power(z,n):
    assert type(n) is int and 0<=n<=256
    result=G.as_z(1);z=G.as_z(z)
    while n:
        if n&1:result=result*z
        n//=2
        if n:z=square(z)
    return result


@lru_cache(2048)
def atan_series(x):
    x=A.bounded(Q(x));assert abs(x)<=Q(1,2)
    if not x:return I.point(0)
    t=I.point(x);term=t;result=term;xx=t.square()
    for j in range(1,192):
        term=-(term*xx);result=result+term/(2*j+1)
    # Alternating remainder for a real point, bounded uniformly on |x|<=1/2.
    tail=Q(1,(1<<385)*385)
    return I(*G.rounded(result.lo-tail,result.hi+tail))


@lru_cache(1)
def pi():
    return 16*atan_series(Q(1,5))-4*atan_series(Q(1,239))


@lru_cache(2048)
def atan_point(x):
    x=A.bounded(Q(x))
    if x<0:return -atan_point(-x)
    if x<=Q(1,2):return atan_series(x)
    if x<=1:return atan_series(Q(1,2))+atan_series((x-Q(1,2))/(1+x/2))
    return pi()/2-atan_point(1/x)


def atan_real(x):
    x=G.as_i(x)
    return I(atan_point(x.lo).lo,atan_point(x.hi).hi)


@lru_cache(2048)
def trig_pi_point(x,cosine):
    x=A.bounded(Q(x));x-=2*((x+1)//2)
    assert -1<=x<1
    if not x:return I.point(int(cosine))
    angle=I.point(x)*pi();xx=angle.square()
    term=I.point(1) if cosine else angle;result=term
    for j in range(1,80):
        den=(2*j-1)*(2*j) if cosine else (2*j)*(2*j+1)
        term=-(term*xx/den);result=result+term
    degree=160 if cosine else 161
    tail=Q(22,7)**degree/factorial(degree)
    return I(max(Q(-1),result.lo-tail),min(Q(1),result.hi+tail))


def trig_pi_real(x,cosine):
    x=G.as_i(x)
    if x.hi-x.lo>=2:return I(Q(-1),Q(1))
    ends=[trig_pi_point(x.lo,cosine),trig_pi_point(x.hi,cosine)]
    lo=min(v.lo for v in ends);hi=max(v.hi for v in ends)
    shift=Q() if cosine else Q(1,2)
    first=-(-(x.lo-shift).numerator//(x.lo-shift).denominator)
    last=(x.hi-shift).numerator//(x.hi-shift).denominator
    for integer in range(first,last+1):
        value=Q(1 if integer%2==0 else -1)
        lo=min(lo,value);hi=max(hi,value)
    return I(lo,hi)


def exp_real(x):
    x=G.as_i(x);assert -256<=x.lo<=x.hi<=256
    y=x/I(*G.EXP.log2())
    return I(G.EXP.exp2(y.lo)[0],G.EXP.exp2(y.hi)[1])


def exponential(z):
    z=G.as_z(z);amplitude=exp_real(z.re)
    phase=z.im/pi()
    return Z(amplitude*trig_pi_real(phase,True),amplitude*trig_pi_real(phase,False))


def sin_cos_pi(z):
    z=G.as_z(z);s=trig_pi_real(z.re,False);c=trig_pi_real(z.re,True)
    phase=z.im*pi();positive=exp_real(phase);negative=exp_real(-phase)
    ch=(positive+negative)/2;sh=(positive-negative)/2
    return Z(s*ch,c*sh),Z(c*ch,-s*sh)


def log_right(z):
    z=G.as_z(z);assert z.re.lo>0,'Entire logarithm rectangle must lie in the right half-plane'
    norm2=z.re.square()+z.im.square()
    real=I(G.BUB.log_bounds(norm2.lo)[0],G.BUB.log_bounds(norm2.hi)[1])/2
    return Z(real,atan_real(z.im/z.re))


def atanhc(z,order=128):
    z=G.as_z(z);assert type(order) is int and 16<=order<=128
    norm=z.norm().hi;assert norm<Q(3,4),'atanhc analytic disk premise'
    result=G.as_z(Q(1,2*order-1))
    for j in range(order-2,-1,-1):result=result*z+Q(1,2*j+1)
    tail=A.bounded(norm**order/((2*order+1)*(1-norm)))
    error=I(-tail,tail)
    return Z(result.re+error,result.im+error)
