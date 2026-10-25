"""Square-root-free five-dimensional cubature diagnostic, never a certificate.

The two child loops are already analytically integrated. q inversion removes
all odd angular terms. The remaining measure is derived from the normalized
S5 two-coordinate disk projection and the full triangle simplex.
"""
from fractions import Fraction as Q
import math
import numpy as np
import sympy as s


def recipe(binding):
    def qs(v):return Q(str(v))
    shifts=[tuple(qs(v) for v in row) for row in binding['host_external_shifts']]
    assert len(shifts)==2
    pairs=[]
    for row in binding['triangle_pairs']:
        c,a,b=map(qs,row['momentum']);norm=a*a-a*b+b*b
        assert c*c in (0,1) and norm in (0,1)
        pairs.append((*row['simplex_pair'],float(c),float(a),float(b),float(norm)))
    assert all(a*a-a*b+b*b==1 for a,b in shifts)
    return {'shifts':[(float(a),float(b)) for a,b in shifts],'pairs':pairs}


def exact_map_identities():
    z,y,v=s.symbols('z y v',real=True);phi=s.symbols('phi',real=True)
    assert s.integrate(2*(1-z),(z,0,1))==1
    assert s.integrate(z*2*(1-z),(z,0,1))*s.integrate(s.cos(phi)**2,(phi,0,s.pi))/s.pi==s.Rational(1,6)
    assert s.integrate(z*z*2*(1-z),(z,0,1))*s.integrate(s.cos(phi)**4,(phi,0,s.pi))/s.pi==s.Rational(1,16)
    jac=s.Matrix([y,(1-y)*v]).jacobian([y,v]).det()
    assert s.expand(jac-(1-y))==0 and s.integrate(jac,(y,0,1),(v,0,1))==s.Rational(1,2)
    r,h,k,E=s.symbols('r h k E')
    assert s.expand((2*r*h)*(2*r*k)-4*r*r*h*k)==0
    # atanh(z)/z is analytic in z^2, including z=0. Series tail is geometric.
    n=12
    assert s.expand(s.series(s.atanh(z),z,0,2*n+1).removeO()-sum(z**(2*j+1)/s.Integer(2*j+1) for j in range(n)))==0
    return {'normalized_disk_measure':True,'sphere_second_and_fourth_moments_exact':True,
            'simplex_jacobian_exact':True,'simplex_volume':'1/2','no_sqrt_of_disk_radius_in_final_kernel':True,
            'atanhc_series':'sum_{j>=0} w^j/(2j+1)',
            'atanhc_real_tail_after_N_terms':'w^N/((2*N+1)*(1-w)) for 0<=w<=1/4',
            'radial_log2_map_density':'t^3 log(2) (1-z) (1-y) G_even',
            'five_dimension_cubature_error_not_certified':True}


def bubble(t):
    assert t>0
    z=t/(t+4)
    if t<=1:
        power=z*z;total=0.
        for j in range(2,65):total+=power/(2*j+1);power*=z
        return (t+4)*total/3
    return (t+4)/3*(math.asinh(math.sqrt(t)/2)/math.sqrt(z)-1)-t/9


def arrays(rec,z,theta,y,v):
    """Precompute angular and simplex factors with NumPy broadcasting."""
    co=np.cos(np.pi*theta);si=np.sin(np.pi*theta);second=(-co+math.sqrt(3)*si)/2
    x=(y,(1-y)*v,(1-y)*(1-v));alpha=0.;C=0.;R=0.
    for i,j,c,a,b,norm in rec['pairs']:
        w=x[i]*x[j];alpha=alpha+w*c*c;C=C+w*norm;R=R+w*c*(a*co+b*second)
    projections=[a*co+b*second for a,b in rec['shifts']]
    return z,y,alpha,C,R,*projections


def mapped(t,cached):
    z,y,alpha,C,R,h,k=cached
    A=1+t;X=t+2;tz=t*z;M=1+t*alpha;E=M+C
    determinant=(X*X-4*tz*h*h)*(X*X-4*tz*k*k)
    U=X*X+4*tz*h*k;V2=4*tz*X*X*(h+k)**2
    He=U/(A*A*determinant)
    Hd=(V2-U*(2*A+1+4*tz*h*k))/(A**4*determinant)
    L2=4*tz*R*R;VL=4*tz*X*(h+k)*R;w=L2/(E*E)
    if np.any(determinant<=0) or np.any(M<=0) or np.any(w<0) or np.any(w>0.250000000001):
        raise ArithmeticError('Pilot left graph-certified real domain')
    delta=(2*M*C+C*C-L2)/(M*M)
    even_log=np.log1p(delta)/2
    # Fixed finite series; its explicit tail is known, but this floating pilot
    # does not turn that tail into an enclosure of the complete cubature.
    atanhc=np.zeros_like(w)+1/65
    for j in range(31,-1,-1):atanhc=atanhc*w+1/(2*j+1)
    cross=-VL/(A*A*determinant*E)*atanhc
    even=-bubble(t)*(Hd*np.log1p(t*alpha)+He*even_log+cross)
    result=t**3*math.log(2)*(1-z)*(1-y)*even
    if np.any(~np.isfinite(result)):raise ArithmeticError('Nonfinite diagnostic kernel')
    return result


def prepare(rec,n):
    assert 2<=n<=24
    node,weight=np.polynomial.legendre.leggauss(n);node=(node+1)/2;weight=weight/2
    cached=arrays(rec,node[:,None,None,None],node[None,:,None,None],node[None,None,:,None],node[None,None,None,:])
    weights=weight[:,None,None,None]*weight[None,:,None,None]*weight[None,None,:,None]*weight[None,None,None,:]
    assert abs(float(weight.sum())-1)<1e-13
    return node,weight,cached,weights


def tile(prepared,left):
    node,weight,cached,weights=prepared
    values=[]
    for x,w in zip(node,weight):values.append(float(w)*float(np.sum(mapped(2.**(left+float(x)),cached)*weights,dtype=np.float64)))
    return math.fsum(values)
