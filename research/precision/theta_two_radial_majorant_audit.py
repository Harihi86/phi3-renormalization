"""Exact rational replay of saved majorant algebra and Cauchy cell allowances.

Geometry and bubble enclosures still require independent analytic-function
replay; this module does not claim to replace that part of the final audit.
"""
from fractions import Fraction as Q
from math import isqrt

CAP=1000000
RADII=(Q(1,3),Q(2,5),Q(1,2),Q(3,5),Q(5,8),Q(2,3))


def bounded(q):
    assert q.numerator.bit_length()<=CAP and q.denominator.bit_length()<=CAP
    return q


def interval(data):
    a,b,e=map(int,data);assert a<=b and abs(e)<=CAP
    assert max(a.bit_length(),b.bit_length())+max(e,0)<=CAP
    return bounded(Q(a)*Q(2)**e),bounded(Q(b)*Q(2)**e)


def sqrt_lower(q,bits=1024):
    q=Q(q);assert q>=0
    n=isqrt((q.numerator << (2*bits))//q.denominator)
    return Q(n,1 << bits)


def propagators(data,N):
    p={k:interval(v) for k,v in data.items()}
    lo=lambda k:p[k][0]
    hi=lambda k:p[k][1]
    for key in ('A0','A1','H1','lt','lu','angular','external_series','radial_scale','bound'):
        assert lo(key)>0,(key,'positive premise')
    assert lo('H0')>=0
    radius=[r for r in RADII if lo('external_radius')<=r<=hi('external_radius')]
    assert len(radius)==1;R=radius[0]
    for name,norm,den in (('eta_t','mt','lt'),('eta_u','mu','lu')):
        e=hi(name);assert 0<=lo(name)<=e<1
        difference=e*lo(den)-R*R
        assert difference>=0 and difference*difference>=4*R*R*hi(norm),(name,'shift bound')
    e=hi('eta_central');assert 0<=lo('eta_central')<=e<1
    A0,A1,H0,H1=lo('A0'),lo('A1'),hi('H0'),hi('H1')
    c=A1/H1;b=A0-c*H0
    assert e>0
    vertex=R/(e*c)
    if vertex*vertex>=H0:
        assert e*b-R*R-R*R/(e*c)>=0,'central quadratic minimum'
    else:
        difference=e*A0-R*R
        assert difference>=0 and difference*difference>=4*R*R*H0,'central endpoint'
    d=sqrt_lower(A0*A0+2*A0*A1)
    angular=4/(A0+A1+d)**2
    assert hi('angular')>=angular,'angle bound'
    series=sum(R**(-2*n) for n in range(1,N+1))
    assert hi('external_series')>=series,'external series'
    scale=hi('mt')**2*hi('mu')**2*hi('bubble_t')*hi('bubble_u')*hi('md')*hi('mapped_jacobian')/(4*lo('lt')**2*lo('lu')**2)
    assert hi('radial_scale')>=scale,'radial Jacobian/scale'
    need=hi('radial_scale')*hi('angular')*hi('external_series')/(lo('lt')*lo('lu')*(1-hi('eta_t'))*(1-hi('eta_u'))*(1-hi('eta_central'))**2)
    assert hi('bound')>=need,'propagator product'
    return hi('bound')


def coefficients(data,N,polynomial=None):
    p={k:interval(v) for k,v in data.items()};lo=lambda k:p[k][0];hi=lambda k:p[k][1]
    for key in ('lt','lu','ld','lad','real_u_lower'):assert lo(key)>0
    for key in ('r_mod','w_mod','bubble_t','bubble_u','mapped_jacobian','numerator_absolute'):assert lo(key)>=0
    if polynomial is not None:
        r,w=hi('r_mod'),hi('w_mod')
        rp=[Q(1)];wp=[Q(1)]
        for j in range(max(row['powers'][0] for row in polynomial)):rp.append(bounded(rp[-1]*r))
        for j in range(max(row['powers'][1] for row in polynomial)):wp.append(bounded(wp[-1]*w))
        value=sum(abs(Q(row['coefficient']))*rp[row['powers'][0]]*wp[row['powers'][1]] for row in polynomial)
        assert hi('numerator_absolute')>=value,'absolute polynomial'
    need=hi('mt')**2*hi('mu')**2*hi('bubble_t')*hi('bubble_u')*hi('mapped_jacobian')*hi('numerator_absolute')
    need/=4*lo('lt')**(2*N+3)*lo('lu')**(2*N+3)*lo('ld')**(4*N-2)*lo('lad')**2
    assert hi('bound')>=need,'coefficient product'
    return hi('bound')


def cell(data,N,n,area,qx,qy,polynomial=None):
    assert data['finite'] and data['maximum_order']==N and data['gauss_order']==n
    assert data['gauss_moments_checked']==2*n
    assert interval(data['area'])[0]<=area<=interval(data['area'])[1]
    assert 0<qx<1 and 0<qy<1
    for key,q in (('ratio_x',qx),('ratio_y',qy)):
        a,b=interval(data[key]);assert a<=q<=b
    assert len(data['majorant_methods'])==len(data['propagator_majorants'])==len(data['coefficient_majorants'])==2
    upper=[]
    for j,method in enumerate(data['majorant_methods']):
        if method=='original_propagators':bound=propagators(data['propagator_majorants'][j],N)
        else:
            assert method=='absolute_coefficients'
            bound=coefficients(data['coefficient_majorants'][j],N,polynomial)
        key='majorant_x' if j==0 else 'majorant_y';native=interval(data[key])[1]
        assert native>=bound;upper.append(native)
    needed=2*area*(upper[0]*qx**(2*n)/(1-qx)+upper[1]*qy**(2*n)/(1-qy))
    eh=interval(data['error_bound'])[1];assert eh>=needed
    a,b=interval(data['central']['real']);lo,hi=interval(data['real'])
    assert lo<=a-eh and hi>=b+eh,'final enclosure'
    il,ih=interval(data['imag']);assert il<=0<=ih
    return {'interval':(lo,hi),'error_upper':eh,'methods':data['majorant_methods']}
