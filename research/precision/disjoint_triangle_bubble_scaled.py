"""Exact identities for host normalization and logarithmic-chart scaling.

The logarithm is continued from real physical coordinates. Re(Nm)>0 and
Re(1+Cs/Nm)>0 fix the residual logarithms; u*L supplies the exponential phase.
No principal-log product rule is assumed across a branch cut.
"""
import sympy as s


def verify_identities(mutation=None):
    t,z,h,k,M,C,R=s.symbols('t z h k M C R')
    A=1+t;X=2+t
    det=(X**2-4*t*z*h**2)*(X**2-4*t*z*k**2)
    U=X**2+4*t*z*h*k;V2=4*t*z*X**2*(h+k)**2
    he=U/(A**2*det)
    hd=(V2-U*(2*A+1+4*t*z*h*k))/(A**4*det)
    ho=-4*t*z*X*(h+k)*R/(A**2*det*(M+C))
    eta=4*t*z/X**2;beta=(2*A+1)/X**2
    da=1-eta*h**2;db=1-eta*k**2;u0=1+eta*h*k
    base=t**3/(A**4*da*db)
    if mutation=='host-power':base*=A
    n=eta*(h+k)**2-u0*(beta+eta*h*k)
    if mutation=='overall-subtraction':n=eta*(h+k)**2+u0*(beta+eta*h*k)
    odd=-(A/X)**2*(4*t*z/X)*(h+k)*R/(M+C)
    if mutation=='odd-sign':odd=-odd
    assert s.cancel(base*n-t**3*hd)==0
    assert s.cancel(base*(A/X)**2*u0-t**3*he)==0
    assert s.cancel(base*odd-t**3*ho)==0
    b,v,n1,n2,n3,r1,r2=s.symbols('b v n1 n2 n3 r1 r2')
    y=(b-1)/t;alpha=y*(1-y)
    norm=n2+v*(n1-n2);projection=r2+v*(r1-r2)
    fullC=alpha*(v*n1+(1-v)*n2)+(1-y)**2*v*(1-v)*n3
    fullR=alpha*(v*r1+(1-v)*r2)
    yb=(1-1/b)/t
    if mutation=='scaled-y':yb=(1+1/b)/t
    nm=1-y+yb
    cs=(1-y)*(yb*norm+(1-y)*v*(1-v)*n3/b)
    rs=(1-y)*yb*projection
    if mutation=='pair-interpolation':rs=(1-y)*yb*(r2+v*(r1+r2))
    ne=nm+cs
    if mutation=='scaled-mass':nm=1+y+yb
    assert s.cancel(b*nm-(1+t*alpha))==0
    assert s.cancel(b*cs-fullC)==0
    assert s.cancel(b*rs-fullR)==0
    assert s.cancel(b*ne-(1+t*alpha+fullC))==0
    assert s.cancel(cs/nm-fullC/(1+t*alpha))==0
    assert s.cancel(rs/ne-fullR/(1+t*alpha+fullC))==0
    return {'host_even_exact':True,'host_overall_difference_exact':True,
            'host_odd_exact':True,'scaled_mass_exact':True,
            'scaled_external_variance_exact':True,'scaled_linear_projection_exact':True,
            'branch_anchor':'real physical u; log(M)=u*L+log(Nm), Re(Nm)>0',
            'no_principal_log_product_identity_assumed':True}


MUTATIONS=('host-power','overall-subtraction','odd-sign','scaled-y','pair-interpolation','scaled-mass')
