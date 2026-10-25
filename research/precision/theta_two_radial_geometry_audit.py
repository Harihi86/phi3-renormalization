"""Independent rational complex rectangles for saved Cauchy majorants.

No native evaluator is called. Exponentials, logarithms, square roots and
trigonometric functions use directed integer series / integer square roots.
The scalar majorant algebra is checked separately by majorant_audit_v4.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
import exact_exp2_enclosures as EXP
import exact_radial_node_bounds_v2 as BUB
from exact_range_enclosures_v2 import sqrt_bounds
import theta_two_radial_majorant_audit_v4 as A

SCALE=1<<384


def rounded(lo,hi):
    lo,hi=Q(lo),Q(hi);assert lo<=hi
    return Q(lo.numerator*SCALE//lo.denominator,SCALE),Q(-(-hi.numerator*SCALE//hi.denominator),SCALE)


@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q
    def __post_init__(self):
        object.__setattr__(self,'lo',A.bounded(Q(self.lo)))
        object.__setattr__(self,'hi',A.bounded(Q(self.hi)))
        assert self.lo<=self.hi
    @classmethod
    def point(cls,x):return cls(Q(x),Q(x))
    def __add__(self,b):
        b=as_i(b);return I(*rounded(self.lo+b.lo,self.hi+b.hi))
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-as_i(b)
    def __rsub__(self,b):return as_i(b)+-self
    def __mul__(self,b):
        b=as_i(b);p=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return I(*rounded(min(p),max(p)))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=as_i(b);assert b.lo>0 or b.hi<0
        return self*I(*rounded(1/b.hi,1/b.lo))
    def square(self):
        low=0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2)
        return I(*rounded(low,max(self.lo**2,self.hi**2)))
    def sqrt(self):
        assert self.lo>=0
        return I(sqrt_bounds(self.lo)[0],sqrt_bounds(self.hi)[1])


def as_i(x):return x if isinstance(x,I) else I.point(x)


@dataclass(frozen=True)
class Z:
    re: I
    im: I=I.point(0)
    def __add__(self,b):
        b=as_z(b);return Z(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return Z(-self.re,-self.im)
    def __sub__(self,b):return self+-as_z(b)
    def __rsub__(self,b):return as_z(b)+-self
    def __mul__(self,b):
        b=as_z(b);return Z(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=as_z(b);den=b.re.square()+b.im.square();assert den.lo>0
        numerator=self*Z(b.re,-b.im)
        return Z(numerator.re/den,numerator.im/den)
    def norm(self):return (self.re.square()+self.im.square()).sqrt()
    def root(self):
        assert self.re.lo>0
        real=I(sqrt_bounds(self.re.lo)[0],sqrt_bounds((self.norm().hi+self.re.hi)/2)[1])
        return Z(real,self.im/(2*real))


def as_z(x):return x if isinstance(x,Z) else Z(as_i(x))


@lru_cache(4096)
def trig_point(x,cosine):
    x=Q(x);assert abs(x)<=1
    if not x:return I.point(int(cosine))
    term=I.point(1 if cosine else x);total=term;square=I.point(x).square()
    for j in range(1,101):
        den=(2*j-1)*(2*j) if cosine else (2*j)*(2*j+1)
        term=-(term*square/den);total=total+term
    # Taylor remainder with real derivatives bounded by one on the real axis.
    error=Q(1,factorial(201 if cosine else 202))
    return I(*rounded(total.lo-error,total.hi+error))


def complex_exp2(x):
    amplitude=I(EXP.exp2(x.re.lo)[0],EXP.exp2(x.re.hi)[1])
    phase=x.im*I(*EXP.log2());assert -1<=phase.lo<=phase.hi<=1
    sine=I(trig_point(phase.lo,False).lo,trig_point(phase.hi,False).hi)
    low_abs=0 if phase.lo<=0<=phase.hi else min(abs(phase.lo),abs(phase.hi))
    high_abs=max(abs(phase.lo),abs(phase.hi))
    cosine=I(trig_point(high_abs,True).lo,trig_point(low_abs,True).hi)
    return Z(amplitude*cosine,amplitude*sine)


@lru_cache(4096)
def bubble_upper(a):
    a=Q(a);assert a>0
    if a<=1:
        term=I.point(a).square()/60;total=term
        for n in range(2,129):
            term=-(term*a*(n-1)/(2*(2*n+3)));total=total+term
        # Successive absolute ratios are <=a/4<=1/4.
        next_term=term*a*128/(2*261)
        return total.hi+max(abs(next_term.lo),abs(next_term.hi))*Q(4,3)
    return BUB.bubble_bounds(a)[1]


def bubble_majorant(z):
    a=z.re.lo;m=z.norm().hi;assert a>-4
    flat=m*m/(60*(1+min(a,Q())/4))
    return min(flat,m*m*bubble_upper(a)/(a*a)) if a>0 else flat


def geometry(x,y,chart):
    r=complex_exp2(x);w=complex_exp2(y) if chart else y
    jac=r*Z(I(*EXP.log2()))
    if chart:jac=jac*w*Z(I(*EXP.log2()))
    u=r*(r+2*w+2)/4;t=u+w;d=1+r+w;ad=(r+2)*(r+2*w+2)/2
    return {'r':r,'w':w,'t':t,'u':u,'d':d,'ad':ad,'jac':jac}


def verify_box(proof,method,x,y,chart):
    p={k:A.interval(v) for k,v in proof.items()};g=geometry(x,y,chart)
    def upper(name,value):assert p[name][1]>=value,(name,'upper geometry premise')
    def lower(name,value):assert p[name][0]<=value,(name,'lower geometry premise')
    t,u,d,ad=g['t'],g['u'],g['d'],g['ad']
    assert t.re.lo>-1 and u.re.lo>-1 and d.re.lo>0
    lower('lt',(1+t).norm().lo);lower('lu',(1+u).norm().lo)
    upper('mt',t.norm().hi);upper('mu',u.norm().hi)
    upper('bubble_t',bubble_majorant(t));upper('bubble_u',bubble_majorant(u))
    upper('mapped_jacobian',g['jac'].norm().hi)
    lower('real_u_lower',u.re.lo)
    if 'real_t_lower' in p:lower('real_t_lower',t.re.lo)
    if method=='absolute_coefficients':
        lower('ld',d.norm().lo);lower('lad',ad.norm().lo)
        upper('r_mod',g['r'].norm().hi);upper('w_mod',g['w'].norm().hi)
    else:
        assert method=='original_propagators' and t.re.lo>0 and u.re.lo>0
        rt,ru=t.root(),u.root();delta=g['w']/(rt+ru);cross=rt*ru
        lower('A0',1+(delta.re.square()-delta.im.square()).lo)
        lower('A1',2*cross.re.lo);upper('H0',delta.norm().hi**2);upper('H1',2*cross.norm().hi)
        upper('md',d.norm().hi);lower('real_d_lower',d.re.lo)
    return {'method':method,'premises_checked':len(p),'independent_integer_transcendentals':True}


def verify(data,request):
    words=request.split();mode=words.pop(0);bits,chart=int(words.pop(0)),int(words.pop(0))
    assert mode in ('range','cell') and bits<=512
    if mode=='cell':n=int(words.pop(0));assert n<=64
    vals=[Q(int(words[j]),int(words[j+1])) for j in range(0,len(words),2)]
    xl,xr,yl,yr=vals[:4]
    assert chart in (0,1) and -128<=xl<xr<=128 and 0<=yl<yr<=(128 if chart else 1)
    x=Z(I(xl,xr));y=Z(I(yl,yr))
    if mode=='range':boxes=[(x,y)]
    else:
        rx,ry=vals[4:];assert (xr-xl)/2<rx<=1 and (yr-yl)/2<ry<=1
        cx=(xl+xr)/2;cy=(yl+yr)/2
        boxes=[(Z(I(cx-rx,cx+rx),I(-rx,rx)),y),(x,Z(I(cy-ry,cy+ry),I(-ry,ry)))]
    assert len(boxes)==len(data['majorant_methods']);out=[]
    for i,((x,y),method) in enumerate(zip(boxes,data['majorant_methods'])):
        key='propagator_majorants' if method=='original_propagators' else 'coefficient_majorants'
        out.append(verify_box(data[key][i],method,x,y,chart))
    return out
