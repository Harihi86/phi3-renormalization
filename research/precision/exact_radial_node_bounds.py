"""Independent integer-series and rational bounds for radial node witnesses."""
from fractions import Fraction as Q
from functools import lru_cache
from exact_exp2_enclosures import exp2,log2,ceil_div,SCALE
from exact_range_enclosures_v2 import sqrt_bounds,interval

@lru_cache(None)
def log_bounds(x):
    x=Q(x);assert x>0
    n=x.numerator.bit_length()-x.denominator.bit_length()
    if x<Q(2)**n:n-=1
    r=x/Q(2)**n;assert 1<=r<2
    z=(r-1)/(r+1);scaled=z*SCALE;zl=scaled.numerator//scaled.denominator;zh=ceil_div(scaled.numerator,scaled.denominator)
    z2l=zl*zl//SCALE;z2h=ceil_div(zh*zh,SCALE);tl,th=zl,zh;sl=sh=0
    for j in range(160):
        sl+=2*tl//(2*j+1);sh+=ceil_div(2*th,2*j+1)
        tl=tl*z2l//SCALE;th=ceil_div(th*z2h,SCALE)
    sh+=ceil_div(2*th*SCALE,321*(SCALE-z2h))
    ll,lh=log2();nl,nh=(n*ll,n*lh) if n>=0 else (n*lh,n*ll)
    return Q(sl,SCALE)+nl,Q(sh,SCALE)+nh

def bubble_bounds(u):
    u=Q(u);assert u>0
    rl,rh=sqrt_bounds(u);sl,sh=sqrt_bounds(u+4)
    ll,_=log_bounds((rl+sl)/2);_,lh=log_bounds((rh+sh)/2)
    assert rl>0 and ll>=0
    return max(Q(),(u+4)/3*(ll*sl/rh-1)-u/9),(u+4)/3*(lh*sh/rl-1)-u/9

def remainder_bounds(u,P):
    u,P=Q(u),Q(P);assert u>=0 and P>0
    D=u+1;difference=D-P;sl,sh=sqrt_bounds(difference*difference+4*P)
    if difference>0:dl,dh=2*P/(difference+sh),2*P/(difference+sl)
    else:dl,dh=(P-D+sl)/2,(P-D+sh)/2
    El,Eh=(D+P+sl)/2,(D+P+sh)/2
    def numerator(d):return P*d*(3*D*D+6*D*d+(u+3)*d*d)
    assert 0<=dl<=dh and El>0
    return numerator(dl)/(3*D**3*Eh**3),numerator(dh)/(3*D**3*El**3)

def verify_range(cell,P,top):
    assert cell['valid'] and cell['method']=='positive_range'
    i,l=cell['geometry'];assert 0<=i<2**l and 0<=l<=30
    width=top+10;xa=Q(width*i,2**l)-10;xb=Q(width*(i+1),2**l)-10;h=xb-xa
    vals=[]
    for v in cell['witnesses']:
        a,b=interval(v);assert a==b and a>0;vals.append(a)
    ul,uh,ln,B,R=vals
    assert ul<=exp2(xa)[0] and uh>=exp2(xb)[1] and ln>=log2()[1]
    assert B>=bubble_bounds(uh)[1] and R>=remainder_bounds(ul,P)[1]
    required=h*ln*uh**3*B*R/(2*(ul+1)**2)
    lo,hi=interval(cell['central_contribution']);el,eh=interval(cell['error_bound'])
    assert lo<=required/2<=hi and eh>=required/2 and el>=0
    return required
