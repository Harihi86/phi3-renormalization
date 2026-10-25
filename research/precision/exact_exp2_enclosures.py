"""Independent rational enclosures using directed fixed-point integer series."""
from fractions import Fraction as Q
from functools import lru_cache
BITS=384
SCALE=1<<BITS
def ceil_div(a,b):return -(-a//b)
@lru_cache(None)
def log2():
    n=160
    lo=sum(2*SCALE//((2*r+1)*3**(2*r+1)) for r in range(n))
    tail=Q(2,(2*n+1)*3**(2*n+1))/Q(8,9)
    hi=lo+n+ceil_div(tail.numerator*SCALE,tail.denominator)
    return Q(lo,SCALE),Q(hi,SCALE)
@lru_cache(None)
def exp2(x):
    x=Q(x);whole=x.numerator//x.denominator;fraction=x-whole
    if not fraction:return (Q(2)**whole,)*2
    ll,lh=log2();xl=(ll*fraction*SCALE).numerator//(ll*fraction*SCALE).denominator
    temp=lh*fraction*SCALE;xh=ceil_div(temp.numerator,temp.denominator)
    assert 0<=xl<=xh<SCALE
    tl=th=sl=sh=SCALE
    for n in range(1,129):
        tl=tl*xl//(SCALE*n);th=ceil_div(th*xh,SCALE*n);sl+=tl;sh+=th
    sh+=ceil_div(th*xh,SCALE*129-xh)
    factor=Q(2)**whole
    return Q(sl,SCALE)*factor,Q(sh,SCALE)*factor
def verify_inputs(g,values):
    def interval(value):
        a,b,e=map(int,value);return Q(a)*Q(2)**e,Q(b)*Q(2)**e
    ix,lx,iy,ly=g;xa=Q(58*ix,2**lx)-10;xb=Q(58*(ix+1),2**lx)-10;ya=Q(59*iy,2**ly)-10;yb=Q(59*(iy+1),2**ly)-10
    parsed=[]
    for value in values:
        lo,hi=interval(value);assert lo==hi and lo>0;parsed.append(lo)
    tl,th,ul,uh,ln=parsed
    assert tl<=exp2(xa)[0] and th>=exp2(xb)[1] and ul<=exp2(ya)[0] and uh>=exp2(yb)[1] and ln>=log2()[1]
    return tl,th,ul,uh,ln,xb,yb
def density_bound(g,values,j,k):
    tl,th,ul,uh,ln,xb,yb=verify_inputs(g,values)
    Bt=th**2/60;Bu=uh**2/60;K=Q(1);R=th/(ul+1)**2
    if tl>=2:Bt=min(Bt,th*xb*ln/6);K=min(K,4/tl**3)
    if ul>=2:Bu=min(Bu,uh*yb*ln/6)
    if ul>=2*max(1,th):R=min(R,4*th**2/ul**5)
    return ln**2*th*uh*th**2*Bt**j*K/(2*(tl+1)**(j+2))*uh**2*Bu**k*R/(2*(ul+1)**(k+1))

def verify_cell(c,j,k):
    def interval(value):
        a,b,e=map(int,value);return Q(a)*Q(2)**e,Q(b)*Q(2)**e
    assert c['valid'] and c['method']=='positive_range'
    ix,lx,iy,ly=c['geometry'];assert 0<=ix<2**lx and 0<=iy<2**ly
    area=Q(58,2**lx)*Q(59,2**ly);assert interval(c['area'])==(area,area)
    required=density_bound(c['geometry'],c['range_inputs'],j,k)
    lo,hi=interval(c['density_upper']);assert lo==hi and hi>=required>0
    cl,ch=interval(c['central_contribution']);el,eh=interval(c['error_bound'])
    assert cl<=hi*area/2<=ch and 0<=el<=eh and eh>=hi*area/2
    return {'exact_required_density_bound':required,'stored_density_bound':hi,'area':area,'upper_integral':hi*area}
