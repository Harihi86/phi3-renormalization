"""Sharper positive range witnesses; frozen exp2 reference remains unchanged."""
from fractions import Fraction as Q
from math import isqrt
import exact_exp2_enclosures as BASE
exp2,log2,verify_inputs=BASE.exp2,BASE.log2,BASE.verify_inputs
SCALE=1<<384
def interval(value):
    a,b,e=map(int,value);return Q(a)*Q(2)**e,Q(b)*Q(2)**e
def sqrt_bounds(x):
    x=Q(x);assert x>=0;scaled=x.numerator*SCALE*SCALE//x.denominator;n=isqrt(scaled)
    lo=Q(n,SCALE);hi=lo if n*n*x.denominator==x.numerator*SCALE*SCALE else Q(n+1,SCALE)
    return lo,hi
def density_bound(g,values,aux,j,k):
    tl,th,ul,uh,ln,xb,yb=verify_inputs(g,values)
    sl,slhi=interval(aux[0]);du,duhi=interval(aux[1]);assert sl==slhi and du==duhi and sl>0 and du>0
    assert sl<=sqrt_bounds(tl*tl+4)[0]
    difference=ul+1-th;rl,rh=sqrt_bounds(difference*difference+4*th)
    delta_upper=2*th/(difference+rl) if difference>0 else (th-ul-1+rh)/2
    assert du>=delta_upper
    Bt=th**2/60;Bu=uh**2/60;K=Q(1);R=th/(ul+1)**2
    if tl>=2:Bt=min(Bt,th*xb*ln/6);K=min(K,4/tl**3)
    if ul>=2:Bu=min(Bu,uh*yb*ln/6)
    if ul>=2*max(1,th):R=min(R,4*th**2/ul**5)
    Dl,Dh=ul+1,uh+1
    R=min(R,th*du*(3*Dh**2+6*Dh*du+(uh+3)*du**2)/(3*Dl**3*max(Dl,tl)**3))
    qh=4/(sl+tl);Dtlo,Dthi=tl+1,th+1;Hlo=tl+2+sl;Jlo=3*sl+tl+2
    Nh=8*Dthi**2+4*Dthi*qh*(3*Dthi-2)+2*qh**2*(8*Dthi-1)+3*qh**3
    K=min(K,Nh/(Dtlo**2*Hlo**2*Jlo))
    return ln**2*th*uh*th**2*Bt**j*K/(2*(tl+1)**(j+2))*uh**2*Bu**k*R/(2*(ul+1)**(k+1))
def verify_cell(c,j,k):
    assert c['valid'] and c['method']=='positive_range'
    ix,lx,iy,ly=c['geometry'];assert 0<=ix<2**lx and 0<=iy<2**ly
    area=Q(58,2**lx)*Q(59,2**ly);assert interval(c['area'])==(area,area)
    required=density_bound(c['geometry'],c['range_inputs'],c['range_aux'],j,k)
    lo,hi=interval(c['density_upper']);assert lo==hi and hi>=required>0
    cl,ch=interval(c['central_contribution']);el,eh=interval(c['error_bound'])
    assert cl<=hi*area/2<=ch and 0<=el<=eh and eh>=hi*area/2
    return {'exact_required_density_bound':required,'stored_density_bound':hi,'area':area,'upper_integral':hi*area}
