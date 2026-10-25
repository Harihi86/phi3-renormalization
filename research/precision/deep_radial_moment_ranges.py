"""Rational witnesses keeping the root decay together with powers of t."""
from fractions import Fraction as Q
from deep_radial_ranges import geometry
from exact_radial_node_bounds_v2 import remainder_bounds,exp2,log2,sqrt_bounds,interval

def moment_upper(tl,th,legs,n):
    tl,th=Q(tl),Q(th);assert 0<tl<=th and legs in (2,3) and n in (0,1,2)
    bounds=[]
    if tl<=2:
        b=min(th,Q(2))
        bounds.append(b**(3+n)/(2*(b+1)**(4 if legs==2 else 2)))
    if th>=2:
        a=max(tl,Q(2))
        if legs==2:bounds.append(2*a**(n-2)/(a+1)**2)
        else:
            b=th if n==2 else a
            bounds.append(2*b**n/(b+1)**2)
    return max(bounds)

def verify(cell,piece,legs):
    assert cell['valid'] and cell['method']=='positive_range';xa,xb,ya,yb=geometry(cell['geometry'],legs)
    area=(xb-xa)*(yb-ya);assert interval(cell['area'])==(area,area)
    assert len(cell['witnesses'])==13;values=[]
    for i,v in enumerate(cell['witnesses']):
        a,b=interval(v);assert a==b and (a>=0 if i==9 else a>0);values.append(a)
    tl,th,ul,uh,ln,width,M,K,R,sl,m0,m1,m2=values
    left,right,ml,mr=[Q(piece[k]) for k in ('left','right','mid_left','mid_right')]
    assert tl<=exp2(xa)[0] and th>=exp2(xb)[1] and ln>=log2()[1]
    assert ul<=left+(right-left)*ya and uh>=left+(right-left)*yb and width>=right-left>0
    assert M>=max(ml+(mr-ml)*ya,ml+(mr-ml)*yb)
    assert R>=remainder_bounds(ul,th)[1]
    if legs==2:assert sl==0 and K>=remainder_bounds(tl,1)[1]
    else:
        assert sl>0 and sl<=sqrt_bounds(tl*tl+4)[0]
        q=4/(sl+tl);Dl,Dh=tl+1,th+1;H=tl+2+sl;J=3*sl+tl+2
        N=8*Dh**2+4*Dh*q*(3*Dh-2)+2*q*q*(8*Dh-1)+3*q**3
        required=min(Q(1),N/(Dl**2*H**2*J))
        if tl>=2:required=min(required,4/tl**3)
        assert K>=required
    for n,v in enumerate((m0,m1,m2)):assert v>=moment_upper(tl,th,legs,n)
    candidates=[th**3*K*R/(2*(th+1)**2),m0*R,m1/(ul+1)**2]
    if ul>=2*max(Q(1),th):candidates.append(4*m2/ul**5)
    required=area*ln*width*uh**2*M/(2*(uh+1)**2)*min(candidates)
    lo,hi=interval(cell['center']);el,eh=interval(cell['error_bound'])
    assert 0<=lo<=hi and 0<=el<=eh and lo-eh<=0 and hi+eh>=required and eh>=hi
    return {'area':area,'exact_required_upper_integral':required}
