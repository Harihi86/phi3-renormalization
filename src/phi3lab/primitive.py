"""Finite BPH integral of the primitive K(2,3) two-loop three-point graph.

All six edges have positive common mass squared M. With loop measure
d^6k/pi^3 d^6q/pi^3, the unweighted scalar integral is logarithmically
divergent and has no UV-divergent proper subgraphs. Thus R I=I(p)-I(0).
The subtracted integral is computed on six smooth five-dimensional charts.
No diagram weight, coupling or (4*pi)^(-6) factor is included.

This module covers this topology for arbitrary momentum-conserving Euclidean
three-point invariants, including degenerate Gram matrices and p=0. It is not
a generic two-loop forest integrator. See docs/primitive-analysis.html.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import exp, fsum, isfinite, log, log1p
from decimal import Decimal, ROUND_CEILING, localcontext

from .finite import _gauss_legendre
from .kinematics import validate_three_point, normalized_gram_vectors, kinematics_status
from .overlap_bounds import _upper_float


def _nonnegative(value, name):
    if isinstance(value,bool):
        raise ValueError(f'{name} must be a real number, not a boolean')
    try:value=float(value)
    except (TypeError,ValueError,OverflowError) as exc:
        raise ValueError(f'{name} must be a finite real number') from exc
    if not isfinite(value) or value<0:
        raise ValueError(f'{name} must be finite and nonnegative')
    return value


def _inputs(invariants,m2):
    invariants=validate_three_point(invariants)
    mass=_nonnegative(m2,'m2')
    if mass==0:raise ValueError('m2 must be positive')
    return invariants,mass


def _log_scaled(value,scale,mass):
    if value<=0 or scale==0:return 0.0
    ratio=scale/mass
    if isfinite(ratio):return log1p(ratio*value)
    log_ratio=log(scale)-log(mass)+log(value)
    if log_ratio>36:return log_ratio+log1p(exp(-log_ratio))
    return log1p(exp(log_ratio))


def _gram_square(a,b,vectors):
    """Nonnegative |a*p1+b*p2|^2 for the shared normalized Gram vectors."""
    p1x,p1y,p2x,p2y=vectors
    return fsum((a*p1x,b*p2x))**2+fsum((a*p1y,b*p2y))**2


def primitive_effective_scale(pair_sums,fractions,invariants):
    """Exact-arithmetic-compatible Gaussian effective external scale w.

    pair_sums=(A,B,C)>0, fractions=(u,v,w) in [0,1], invariants=(s1,s2,s3).
    A+B+C need not be one. Scalars may be Fraction to crosscheck U/F exactly.
    The formula uses p3=-p1-p2 and is homogeneous degree one in pair_sums.
    """
    A,B,C=pair_sums;u,v,w=fractions;s1,s2,s3=invariants
    if min(A,B,C)<=0:raise ValueError('Three positive branch-pair sums are required')
    if any(x<0 or x>1 for x in (u,v,w)):raise ValueError('Branch fractions must lie in [0,1]')
    U=A*B+A*C+B*C
    square=s1*(u-w)**2+s2*(v-w)**2+(s3-s1-s2)*(u-w)*(v-w)
    return A*u*(1-u)*s1+B*v*(1-v)*s2+C*w*(1-w)*s3+A*B*C/U*square


def primitive_chart_density(r,t,u,v,w,invariants=(1,1,1),m2=1):
    """One ordered branch chart, including its parameter Jacobian.

    Pair sums=(r*t,r,1)/(1+r+r*t). To obtain the full integral sum the six
    permutations of the branch invariants. All chart faces have finite limits.
    """
    invariants,mass=_inputs(invariants,m2)
    coordinates=tuple(_nonnegative(x,'chart coordinate') for x in (r,t,u,v,w))
    if max(coordinates)>1:raise ValueError('Chart coordinates must lie in [0,1]')
    r,t,u,v,w=coordinates
    scale=max(invariants)
    if scale==0 or t==0:return 0.0
    s1,s2,s3=(s/scale for s in invariants)
    denominator=1+t+r*t
    total=1+r+r*t
    square=_gram_square(u-w,v-w,normalized_gram_vectors(invariants))
    effective=(r*t*u*(1-u)*s1+r*v*(1-v)*s2+w*(1-w)*s3+r*t/denominator*square)/total
    return -t/denominator**3*_log_scaled(effective,scale,mass)


@lru_cache(maxsize=16)
def _smooth_nodes(order):
    # Both endpoint logarithms acquire a quadratic approach. Normalize the
    # positive weights to the exact uniform interval mass, so a common large
    # logarithm cannot amplify quadrature error in that known constant.
    transformed=[]
    for x,weight in _gauss_legendre(order):
        denominator=x*x+(1-x)*(1-x)
        transformed.append((x*x/denominator,weight*2*x*(1-x)/denominator**2))
    mass=fsum(weight for _,weight in transformed)
    return tuple((x,weight/mass) for x,weight in transformed)


@lru_cache(maxsize=8)
def _phase_nodes(order,shape,smoothed=False):
    nodes=_smooth_nodes(order) if smoothed else _gauss_legendre(order)
    scale=max(shape)
    s1,s2,s3=(s/scale for s in shape)
    vectors=normalized_gram_vectors(shape)
    result=[]
    for u,wu in nodes:
        for v,wv in nodes:
            for w,ww in nodes:
                result.append((u*(1-u)*s1,v*(1-v)*s2,w*(1-w)*s3,
                               _gram_square(u-w,v-w,vectors),wu*wv*ww))
    return tuple(result)


@lru_cache(maxsize=16)
def _discrete_period(order,smoothed=False):
    nodes=_smooth_nodes(order) if smoothed else _gauss_legendre(order)
    return 6*fsum(wr*wt*t/(1+t+r*t)**3 for r,wr in nodes for t,wt in nodes)


def primitive_bph(invariants=(1,1,1),m2=1,*,order=12):
    """Unweighted finite (1-t^0)J for the primitive two-loop K(2,3) graph.

    The deterministic tensor order is per coordinate, not the total point
    count. Six smooth charts each have five integration coordinates. Equal
    branch-invariant permutations are combined before quadrature.
    """
    invariants,mass=_inputs(invariants,m2)
    if type(order) is not int or not 2<=order<=32:
        raise ValueError('Five-dimensional quadrature order must be an integer from 2 to 32')
    scale=max(invariants)
    if scale==0:return 0.0
    smoothed=log(scale)-log(mass)>log(100)
    nodes=_smooth_nodes(order) if smoothed else _gauss_legendre(order)
    # Canonical external relabelling makes permutation symmetry exact also
    # for finite quadrature orders and permits repeated branch charts to merge.
    # Retain the original binary64 scales until the exact Gram reconstruction;
    # normalizing a tiny leg first can underflow its invariant to zero.
    normalized=tuple(sorted(invariants))
    total=[]
    for shape,multiplicity in sorted(Counter(permutations(normalized)).items()):
        phases=_phase_nodes(order,shape,smoothed)
        chart=[]
        for r,wr in nodes:
            for t,wt in nodes:
                denominator=1+t+r*t
                inverse_total=1/(1+r+r*t)
                a=r*t*inverse_total;b=r*inverse_total;c=inverse_total
                correction=r*t/denominator*inverse_total
                value=fsum(weight*_log_scaled(a*hu+b*hv+c*hw+correction*square,scale,mass)
                           for hu,hv,hw,square,weight in phases)
                chart.append(-wr*wt*t/denominator**3*value)
        total.append(multiplicity*fsum(chart))
    # The projective period is analytically1/2. Normalize the positive tensor
    # weights to that exact mass. Otherwise a small period quadrature error is
    # multiplied by a very large common log(smax/m2), and can even violate the
    # proven magnitude bound at coarse order. This is a convergent ratio rule,
    # not multiplication by the diagram symmetry weight, which stays excluded.
    return .5*fsum(total)/_discrete_period(order,smoothed)


def primitive_bound(invariants=(1,1,1),m2=1):
    """Explicit amplitude/scale bounds, distinct from numerical error bars.

    The exact positive projective period is1/2 and w<=max(s_i)/4. The sharp
    logarithmic magnitude bound is .5*log(1+max(s_i)/(4*m2)). Derivatives
    refer to a ray lambda->(lambda*s1,lambda*s2,lambda*s3), at lambda=1.
    """
    kinematics=kinematics_status(invariants)
    invariants,mass=_inputs(kinematics['effective_invariants'],m2)
    scale=max(invariants)
    sd,md=Decimal.from_float(scale),Decimal.from_float(mass)
    with localcontext() as context:
        context.prec=max(100,len(sd.as_tuple().digits)+4,len(md.as_tuple().digits)+4)
        context.rounding=ROUND_CEILING
        ratio=sd/md
        log_bound=(1+sd/(4*md)).ln().next_plus()/2 if scale else Decimal(0)
        small=ratio/8
        first=min(small,Decimal('.5'))
        second=2*first*first
        global_gradient=1/(2*(2*md).sqrt().next_minus())
        compact_gradient=sd.sqrt().next_plus()/(Decimal(2).sqrt().next_minus()*md) if scale else Decimal(0)
        gradient=min(global_gradient,compact_gradient)
        hessian=Decimal(3)/(4*md)
        bounds={'absolute_integral_bound':min(log_bound,small),'small_scale_bound':small,
                'ray_scale_first_derivative_bound':first,'ray_scale_second_derivative_bound':second,
                'independent_momentum_gradient_norm_bound':gradient,
                'independent_momentum_hessian_norm_bound':hessian}
        floats={name:_upper_float(value) for name,value in bounds.items()}
        decimals={name:str(value) for name,value in bounds.items()}
    return {'period_exact':'1/2',**floats,'decimal_upper_bounds':decimals,
            'kinematics':kinematics,
            'formula':'|J_R| <= 0.5*log(1+smax/(4*m2)) <= smax/(8*m2)',
            'scope':'Primitive two-loop K(2,3); no UV proper subgraphs; arbitrary legal Euclidean three-point invariants',
            'derivative_coordinates':'Ray derivatives: lambda*(s1,s2,s3) at lambda=1. Momentum norm derivatives: independent vectors(p1,p2), p3=-p1-p2.',
            'momentum_absolute_convergence':'Proved separately by hard/soft branch decomposition in docs/primitive-analysis.html#momentum-proof',
            'normalization':'d^6k/pi^3 d^6q/pi^3; no 1/S, coupling, or (4*pi)^(-6)',
            'not_an_error_bound':True,
            'numeric_rounding':'Upward Decimal enclosures, then outward binary64 conversion; enclosures apply to converted binary64 inputs.'}


def primitive_benchmark(orders=(8,14,20),invariants=(1,1,1),m2=1):
    """Executed quadrature refinements and analytic bounds for this topology."""
    kinematics=kinematics_status(invariants)
    invariants,mass=_inputs(kinematics['effective_invariants'],m2)
    orders=tuple(orders)
    if len(orders)<2 or any(a>=b for a,b in zip(orders,orders[1:])):
        raise ValueError('Use at least two strictly increasing quadrature orders')
    samples=[{'order':n,'value':primitive_bph(invariants,mass,order=n),
              'period_before_exact_normalization':_discrete_period(n)} for n in orders]
    bound=primitive_bound(kinematics['requested_invariants'],mass)
    value=samples[-1]['value']
    difference=abs(value-samples[-2]['value'])
    return {'id':'L2E3G1','name':'primitive_K2_3','dimension':6,
            'invariants':list(invariants),'m2':mass,'subtraction':'R=(1-t_G^0); no proper UV subgraphs',
            'kinematics':kinematics,
            'value':value,'quadrature':samples,
            'last_refinement_difference':difference,
            'difference_interpretation':'Quadrature refinement diagnostic, not a rigorous numerical-error enclosure',
            'bounds':bound,'bare_DR_simple_pole_coefficient':'1/4 (d=6-2*eps; unweighted)',
            'checks':{'nonpositive':value<=0,'within_analytic_magnitude_bound':abs(value)<=bound['absolute_integral_bound']+1e-14,
                      'quadrature_refinement_stable':difference<=1e-10*max(1,abs(value)),
                      'zero_subtraction_point':primitive_bph((0,0,0),mass,order=2)==0},
            'integration_coordinates':'Six ordered branch-pair charts on[0,1]^5; positive tensor Gauss-Legendre weights normalized to the exact period1/2; identical invariant permutations merged',
            'limits':'Complete finite integral for this primitive topology; not a generic forest integrator or general all-loop BPH proof.'}
