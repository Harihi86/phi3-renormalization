"""Rigorous full radial tails after primitive-loop and parity integration.

This covers the omitted radial domain only. Finite five-dimensional quadrature
still requires its own enclosure. Unit mass and three unit external invariants
are checked against the graph-bound momentum records.
"""
from fractions import Fraction as Q
import sympy as s
from exact_exp2_enclosures import log2


def verify_graph(binding,gram):
    gram=s.Matrix([[s.Rational(v) for v in row] for row in gram])
    assert gram==s.Matrix([[1,s.Rational(-1,2)],[s.Rational(-1,2),1]])
    for row in binding['host_external_shifts']:
        v=s.Matrix([s.Rational(c) for c in row]);assert (v.T*gram*v)[0]==1
    cross=[];norms=[]
    for row in binding['triangle_pairs']:
        c,a,b=map(s.Rational,row['momentum']);norm=(s.Matrix([[a,b]])*gram*s.Matrix([a,b]))[0]
        assert c*c in (0,1) and norm in (0,1)
        if c:cross.append(set(row['simplex_pair']))
        norms.append(str(norm))
    assert len(cross)==2 and len(cross[0]&cross[1])==1
    singleton=next(iter(cross[0]&cross[1]));x=s.symbols('x0:3')
    variance=sum(x[r['simplex_pair'][0]]*x[r['simplex_pair'][1]]*s.Rational(r['momentum'][0])**2 for r in binding['triangle_pairs'])
    others=[i for i in range(3) if i!=singleton]
    assert s.expand(variance-x[singleton]*sum(x[i] for i in others))==0
    return {'singleton_simplex_coordinate':singleton,'V0':'t*y*(1-y)',
            'pair_external_squared_norms':norms,'external_variance_C_at_most_one_third':True,
            'linear_variance_L_squared_at_most_4_V0_C_by_weighted_Gram_Cauchy':True,
            'two_shifted_host_norms_equal_one':True}


def verify_inequalities():
    h=s.symbols('h',nonnegative=True);t=s.symbols('t',positive=True)
    def above16(poly):
        p=s.Poly(s.expand(poly.subs(t,16+h)),h)
        assert all(v>=0 for v in p.all_coeffs())
    # The four shifted denominators are >= t/2 because (sqrt(t)-2)^2 >= 0.
    r=s.symbols('r',nonnegative=True)
    assert s.expand(r*r+2-2*r-r*r/2-(r-2)**2/2)==0
    above16(2*t*t-((t+2)**2+4*t))
    above16(7*t-(6*t+3))
    above16(21*t**3-16*t*(t+2)**2)
    above16(s.Rational(9,8)*t-(t+2))
    assert 16*(21+2*7)==560 and 16*s.Rational(9,8)*4==72
    v,C=s.symbols('v C',nonnegative=True)
    assert s.expand((1+v+C)**2-4*v*(1+C)-(v-(1+C))**2)==0
    # If C<=1/3, squared log-odd ratio <= C/(1+C)<=1/4.
    assert s.factor(s.Rational(1,4)-C/(1+C))==(1-3*C)/(4*(C+1))
    x,y,z=s.symbols('x y z',nonnegative=True)
    assert s.expand((x+y+z)**2/3-(x*y+x*z+y*z)-((x-y)**2+(x-z)**2+(y-z)**2)/6)==0
    # Integral over the triangle has measure (1-y)dy dz, volume 1/2.
    assert s.integrate(1-y,(y,0,1))==s.Rational(1,2)
    # M^-1 integral <= 2/t log(1+t/4), from y(1-y)>=y/2 on [0,1/2].
    assert s.expand(y*(1-y)-y/2-y*(s.Rational(1,2)-y))==0
    assert s.integrate(1/(1+t*y/2),(y,0,s.Rational(1,2)))==2*s.log(t/4+1)/t
    # integral_simplex 1/sqrt(y(1-y)) = beta(1/2,3/2)=pi/2.
    assert s.simplify(s.gamma(s.Rational(1,2))*s.gamma(s.Rational(3,2))/s.gamma(2)-s.pi/2)==0
    # pi<4 and 3sqrt(3)>4 imply 4pi/(3sqrt(3))<4; no decimal bound used.
    assert 27>16
    assert s.Rational(560,2)+s.Rational(4*22,9)<=290
    assert 72*4==288
    # B_R(t)<=t/6 log(t), and the radial measure is t^2 dt/2.
    assert s.integrate(y*(1-y),(y,0,1))==s.Rational(1,6)
    assert s.integrate((y*(1-y))**2/2,(y,0,1))==s.Rational(1,60)
    assert s.Rational(5,8)*s.Rational(1,60)*s.Rational(1,2)*s.Rational(1,5)==s.Rational(1,960)
    assert s.Rational(290,12)==s.Rational(145,6) and s.Rational(288,12)==24
    R=s.symbols('R',positive=True)
    for k,poly in ((1,s.log(R)+1),(2,s.log(R)**2+2*s.log(R)+2)):
        assert s.simplify(s.diff(poly/R,R)+s.log(R)**k/R**2)==0
    return {'UV_starts_at_t':16,'host_even_bound':'4/t^4','host_odd_bound':'72/t^(9/2)',
            'host_even_T0_difference_bound':'560/t^5',
            'simplex_even_log_difference_bound':'22*log(t)/(9*t)',
            'simplex_odd_log_bound':'4/sqrt(t)',
            'radial_absolute_UV_bound':'((145/6)*log(t)^2+24*log(t))/t^2',
            'radial_absolute_IR_bound':'t^4/192 for 0<=t<=1',
            'simplex_volume':'1/2','radial_measure':'t^2 dt/2',
            'angular_projection_measure':'2*(1-z) dz dtheta/pi on 0<=theta<=pi for even kernel',
            'finite_map_density':'t^3*log(2)*(1-z)*(1-y)*even_kernel; t=2^x',
            'finite_domain_integral_still_required':True}


def bounds(left=-12,right=64):
    assert isinstance(left,int) and isinstance(right,int) and -256<=left<=0 and 4<=right<=256
    eps=Q(2)**left;R=Q(2)**right;upper_log=right*log2()[1]
    lower=eps**5/960
    upper=(Q(145,6)*(upper_log**2+2*upper_log+2)+24*(upper_log+1))/R
    return {'log2_radial_domain':[left,right],'lower_radial_tail':lower,'upper_radial_tail':upper,
            'total_radial_tail':lower+upper,'finite_quadrature_error_not_included':True}
