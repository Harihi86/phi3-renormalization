"""Independent rational audit of positive boxes for affine inner envelopes."""
from fractions import Fraction as Q
from exact_radial_node_bounds_v2 import remainder_bounds,exp2,log2,sqrt_bounds,interval

def geometry(g,legs):
    pi,ix,lx,iy,ly=g
    assert all(isinstance(v,int) for v in g) and pi>=0 and 0<=lx<=30 and 0<=ly<=30 and 0<=ix<2**lx and 0<=iy<2**ly
    width=26 if legs==2 else 58;assert legs in (2,3)
    return Q(width*ix,2**lx)-10,Q(width*(ix+1),2**lx)-10,Q(iy,2**ly),Q(iy+1,2**ly)

def verify(cell,piece,legs):
    assert cell['valid'] and cell['method']=='positive_range';xa,xb,ya,yb=geometry(cell['geometry'],legs)
    area=(xb-xa)*(yb-ya);assert interval(cell['area'])==(area,area)
    values=[]
    for i,v in enumerate(cell['witnesses']):
        a,b=interval(v);assert a==b and (a>0 if i<9 else a>=0);values.append(a)
    tl,th,ul,uh,ln,width,M,K,R,sl=values
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
    bound=area*ln*width*th**3*K/(2*(th+1)**2)*uh**2*M*R/(2*(uh+1)**2)
    lo,hi=interval(cell['center']);el,eh=interval(cell['error_bound'])
    assert 0<=lo<=hi and 0<=el<=eh and lo-eh<=0 and hi+eh>=bound and eh>=hi
    return {'area':area,'exact_required_upper_integral':bound}
