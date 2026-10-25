"""Log-domain special functions and exact-A finite BPH triangle quadrature.

The scalar special functions are shared internal numerical utilities. They
use only the standard library and never form exp(logz) at a large logz.
"""
from functools import lru_cache
from math import exp, expm1, frexp, fsum, hypot, inf, isfinite, log, log1p, pi, sqrt

from .finite import _gauss_legendre
from .kinematics import normalized_gram_vectors


def _scale_log_ratio(scale,mass):
    """log(scale/mass) without overflow, invariant under power-of-two scaling."""
    if not scale:return -inf
    a,ea=frexp(scale);b,eb=frexp(mass)
    return log(a/b)+(ea-eb)*log(2.)


def _li2_unit(value):
    """Real Li_2(x) on[0,1], using reflection and a positive <=1/2 series."""
    if not 0<=value<=1:raise ValueError('The unit-interval dilogarithm requires0<=x<=1')
    if value==0:return 0.
    if value==1:return pi*pi/6
    reflected=value>.5
    x=1-value if reflected else value
    terms=[];power=x
    for k in range(1,65):
        term=power/(k*k);terms.append(term)
        if term<2e-18:break
        power*=x
    series=fsum(terms)
    return pi*pi/6-log(value)*log1p(-value)-series if reflected else series


def _softplus(value):
    return value+log1p(exp(-value)) if value>0 else log1p(exp(value))


def _log_f(logz):
    """f(z)=(1+1/z)log(1+z)-1, supplied log(z), including z=0."""
    if logz==-inf:return 0.
    if logz==inf:return inf
    if not isfinite(logz):raise ValueError('logz must be real')
    if logz < -2:
        z=exp(logz);value=0.
        for k in range(32,0,-1):value=z*((-1)**(k+1)/(k*(k+1))+value)
        return value
    logarithm=_softplus(logz)
    return (1+exp(-logz))*logarithm-1


def _log_F(logz):
    """F(z)=-Li_2(-z)-f(z), supplied log(z), stable through logz1400."""
    if logz==-inf:return 0.
    if logz==inf:return inf
    if not isfinite(logz):raise ValueError('logz must be real')
    if logz < -2:
        z=exp(logz);value=0.
        for k in range(32,0,-1):value=z*((-1)**(k+1)/(k*k*(k+1))+value)
        return value
    logarithm=_softplus(logz)
    q=-expm1(-logarithm)
    return _li2_unit(q)+.5*logarithm*logarithm-(1+exp(-logz))*logarithm+1


@lru_cache(maxsize=24)
def _phase_nodes(order,power):
    """Beta-CDF endpoint maps, retaining both complements before rounding."""
    result=[]
    for x,weight in _gauss_legendre(order):
        if power==4:
            cdf=lambda t:t**4*(35+t*(-84+t*(70-20*t)))
            jacobian=140*x**3*(1-x)**3
        elif power==3:
            cdf=lambda t:t**3*(10+t*(-15+6*t))
            jacobian=30*x*x*(1-x)**2
        else:raise ValueError('Endpoint power must be3or4')
        complement=cdf(1-x) if x>=.5 else 1-cdf(x)
        value=cdf(x) if x<=.5 else 1-cdf(1-x)
        result.append((value,complement,weight*jacobian))
    return tuple(result)


@lru_cache(maxsize=32)
def _outer_data(shape,vectors,order,power):
    s1,s2,s3=shape;p1x,p1y,p2x,p2y=vectors
    nodes=_phase_nodes(order,power);result=[]
    for a,outerd,wa in nodes:
        # For s3=0, p2=-p1 and the second outer fraction disappears exactly.
        phases=((.5,.5,1.),) if s3==0 else nodes
        for y,yc,wy in phases:
            b=outerd*y;c=outerd*yc
            W=a*b*s2+a*c*s1+b*c*s3
            result.append((W,c*p1x-b*p2x,c*p1y-b*p2y,wa*wy*outerd))
    return tuple(result)


@lru_cache(maxsize=32)
def _shape_coefficient(shape,power):
    """b=.5*Integral_Delta logW-1/8, independent of the momentum scale."""
    if min(shape)==0:
        # Every legal zero-leg shape, normalized by max(si), is(0,1,1)
        # up to permutation. Its massless triangle logarithmic moment is-1.
        return -.625,-1.,0.
    values=[]
    for order in (128,256):
        nodes=_phase_nodes(order,power);terms=[]
        for a,outerd,wa in nodes:
            for y,yc,wy in nodes:
                b=outerd*y;c=outerd*yc
                W=a*b*shape[1]+a*c*shape[0]+b*c*shape[2]
                terms.append(wa*wy*outerd*log(W))
        values.append(fsum(terms))
    return .5*values[-1]-.125,values[-1],abs(values[-1]-values[-2])


def _j_series(value):
    answer=0.
    for k in range(16,0,-1):
        coefficient=(-1)**(k+1)*(1/(k+1)+1/(k+2))/k
        answer=value*(coefficient+answer)
    return answer


def _j_positive_correction(logR):
    """j(R-1)-1.5logR, without forming a large R."""
    if logR<.05:return _j_series(expm1(logR))-1.5*logR
    eta=exp(-logR);inverse=eta/(-expm1(-logR))
    return (inverse-.5*inverse*inverse)*logR-1.25+.5*inverse


def _j_negative(logQ):
    """j(Q-1), supplied logQ<=0; j(z)=Integral(1+x)log(1+zx)dx."""
    if logQ>-.05:return _j_series(expm1(logQ))
    Q=exp(logQ);complement=-expm1(logQ)
    return Q*(1.5*Q-2)/(complement*complement)*logQ-1.25-.5/complement


def _logaddexp(a,b):
    maximum=max(a,b)
    return maximum+log1p(exp(-abs(a-b)))


def _numpy_functions(np):
    def li2(value):
        value=np.clip(value,0.,1.)
        reflection=value>.5;x=np.where(reflection,1-value,value)
        series=np.zeros_like(x)
        for k in range(48,0,-1):series=x*(1/(k*k)+series)
        with np.errstate(divide='ignore',invalid='ignore'):
            reflected=pi*pi/6-np.log(value)*np.log1p(-value)-series
        reflected=np.where(value==1,pi*pi/6,reflected)
        return np.where(reflection,reflected,series)
    def jseries(value):
        answer=np.zeros_like(value)
        for k in range(16,0,-1):
            coefficient=(-1)**(k+1)*(1/(k+1)+1/(k+2))/k
            answer=value*(coefficient+answer)
        return answer
    def jpos(logR):
        eta=np.exp(-logR);complement=-np.expm1(-logR)
        inverse=np.divide(eta,complement,out=np.zeros_like(eta),where=complement!=0)
        large=(inverse-.5*inverse*inverse)*logR-1.25+.5*inverse
        small=jseries(np.expm1(np.minimum(logR,.1)))-1.5*logR
        return np.where(logR<.05,small,large)
    def jneg(logQ):
        Q=np.exp(logQ);complement=-np.expm1(logQ)
        with np.errstate(divide='ignore',invalid='ignore'):
            large=Q*(1.5*Q-2)/complement**2*logQ-1.25-.5/complement
        return np.where(logQ>-.05,jseries(np.expm1(logQ)),large)
    return li2,jpos,jneg


@lru_cache(maxsize=48)
def _components(shape,vectors,log_ratio,order,backend,power):
    """Exact A integration plus bounded residual, with no large ratio formed."""
    nodes=_phase_nodes(order,power)
    outer=_outer_data(shape,vectors,order,power)
    linear,ilog,difference=_shape_coefficient(shape,power)
    L=log_ratio;mu=exp(-L) if L<746 else 0.
    p1x,p1y=vectors[:2]
    result=[]
    if backend=='numpy':
        import numpy as np
        data=np.asarray(outer);W,kx,ky,weights=(data[:,i] for i in range(4))
        logW=np.log(W);logmuW=np.log(mu+W)
        li2,jpos,jneg=_numpy_functions(np)
        for u,uc,wu in nodes:
            H=u*uc;logH=log(H)
            innerconstant=_li2_unit(1-H)-_j_negative(logH)
            for z,zc,wz in nodes:
                T=uc*z*zc*shape[0]
                Z=H*((kx-zc*p1x)**2+(ky-zc*p1y)**2)
                A=mu+T;q1=-(1-H)*mu+H*W+Z-T
                disc=np.hypot(q1,2*np.sqrt(A*Z))
                V=.5*((1+H)*mu+H*W+T+Z)+.5*disc
                logV=np.log(V);logA=np.logaddexp(-L,log(T) if T else -np.inf)
                logR=np.maximum(0.,logV-logA)
                logQ=np.minimum(0.,logH+logmuW-logV)
                wlog=logH+logmuW
                residual=(.5*wlog*wlog+L*np.log1p(mu/W)+1.5*logV+innerconstant
                          +li2(-np.expm1(-logR))-li2(-np.expm1(logQ))-.5*logQ*logQ
                          +jpos(logR)+jneg(logQ))
                result.append(wu*uc*wz*float(np.dot(weights,residual)))
    else:
        for u,uc,wu in nodes:
            H=u*uc;logH=log(H)
            innerconstant=_li2_unit(1-H)-_j_negative(logH)
            for z,zc,wz in nodes:
                T=uc*z*zc*shape[0]
                logA=_logaddexp(-L,log(T)) if T else -L
                terms=[]
                for W,kx,ky,weight in outer:
                    Z=H*((kx-zc*p1x)**2+(ky-zc*p1y)**2)
                    A=mu+T;q1=-(1-H)*mu+H*W+Z-T
                    disc=hypot(q1,2*sqrt(A*Z))
                    V=.5*((1+H)*mu+H*W+T+Z)+.5*disc
                    logV=log(V);logmuW=log(mu+W)
                    logR=max(0.,logV-logA);logQ=min(0.,logH+logmuW-logV)
                    wlog=logH+logmuW
                    residual=(.5*wlog*wlog+L*log1p(mu/W)+1.5*logV+innerconstant
                              +_li2_unit(-expm1(-logR))-_li2_unit(-expm1(logQ))-.5*logQ*logQ
                              +_j_positive_correction(logR)+_j_negative(logQ))
                    terms.append(weight*residual)
                result.append(wu*uc*wz*fsum(terms))
    residual=fsum(result)
    return {'value':fsum((.125*L*L,linear*L,residual)),
            'log_ratio':L,'quadratic_log_coefficient':.125,'linear_log_coefficient':linear,
            'massless_outer_log_moment':ilog,'coefficient_refinement_difference':difference,
            'residual':residual,'order':order,'endpoint_power':power,'backend':backend,
            'integration_dimension':3 if shape[2]==0 else 4,
            'representation':'Exact analytic A integral; real dilogarithms; explicit scale polynomial plus bounded residual',
            'mass_ratio_underflow_in_bounded_coefficients':mu==0,
            'asymptotic_finite_mass_remainder_upper':130*exp(-L/4),
            'asymptotic_bound_scope':'Whole-phase finite-mass remainder relative to the exact log polynomial and massless residual; uniform over physical normalized shapes. Not a quadrature-error bound.'}


def exact_a_components(invariants,m2,*,order=32,backend='python',power=4):
    """Internal high-scale entry point; inputs have passed shared validation."""
    scale=max(invariants)
    shape=tuple(value/scale for value in invariants)
    vectors=normalized_gram_vectors(invariants)
    L=_scale_log_ratio(scale,m2)
    return dict(_components(shape,vectors,L,order,backend,power))
