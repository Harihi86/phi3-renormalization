"""Finite BPH vertex correction inserted into the two-loop three-point graph.

The proper triangle and overall graph both have Taylor degree0. The four
recursive terms use the same six original Schwinger labels and positive mass.
Gaussian and radial integration produce a paired logarithmic remainder. Four
blow-up charts resolve its single logarithmic UV face and adjacent corner.
No separate divergent term is integrated or assigned a finite MS value.

J has measure d6k/pi3 d6q/pi3, no1/S or coupling. s1 is the physical external
leg attached to the inserted triangle; s2,s3 are interchangeable outer legs.
"""
from decimal import Decimal, ROUND_CEILING, localcontext
from functools import lru_cache
from math import exp, expm1, fsum, hypot, inf, isfinite, log, log1p, nextafter, pi, sin, sqrt
from time import monotonic
from .finite import _gauss_legendre
from .kinematics import normalized_gram_vectors, validate_three_point
from .nested_triangle_extreme import exact_a_components, _scale_log_ratio

GRAPH_LEG_MAP={'L2E3G2':(0,1,2),'L2E3G3':(1,0,2),'L2E3G5':(2,0,1)}


def _number(value,name,positive=False):
    if isinstance(value,bool):raise ValueError(f'{name} must be a real number, not a boolean')
    try:value=float(value)
    except (TypeError,ValueError,OverflowError) as exc:raise ValueError(f'{name} must be finite and real') from exc
    if not isfinite(value) or value<0 or (positive and not value):
        raise ValueError(f'{name} must be finite and '+('positive' if positive else 'nonnegative'))
    return value


def _invariants(s1,s2,s3):
    return validate_three_point((s1,s2,s3))


def _ell(z):
    """log1p(z)/z with its continuous value1 at zero."""
    if z<1e-5:
        return 1+z*(-.5+z*(1/3+z*(-.25+z*.2)))
    return log1p(z)/z


def _outer_shapes(invariants,nodes):
    s1,s2,s3=invariants;scale=max(invariants)
    if not scale:return ()
    p1x,p1y,p2x,p2y=normalized_gram_vectors(invariants)
    result=[]
    for x,wx in nodes:
        for y,wy in nodes:
            a=x;b=(1-x)*y;c=(1-x)*(1-y)
            W=a*b*s2+a*c*s1+b*c*s3
            norm0=scale*((c*p1x-b*p2x)**2+(c*p1y-b*p2y)**2)
            norm1=scale*(((c-1)*p1x-b*p2x)**2+((c-1)*p1y-b*p2y)**2)
            result.append((W,norm0,norm1,wx*wy*(1-x)))
    return tuple(result)


def _sine_nodes(nodes):
    """Endpoint smoothing x=sin(pi*z/2)^2, with its exact Jacobian."""
    return tuple((sin(pi*x/2)**2,weight*pi*sin(pi*x)/2) for x,weight in nodes)


@lru_cache(maxsize=64)
def _python_normalized(s1,s2,s3,order,smooth_a=False):
    nodes=_gauss_legendre(order)
    # At large external scale, an unresolved A~1/max(si) layer gives finite
    # but inaccurate plain tensor estimates. Integrate A on a logarithmic
    # scale. In the A<=r chart integrate r in [A,1] conditionally, so the
    # transition A*W~1 depends on one integration coordinate rather than a
    # diagonal product of two logarithmic coordinates.
    logarithmic=max(s1,s2,s3)>100
    shapes=_outer_shapes((s1,s2,s3),_sine_nodes(nodes) if logarithmic else nodes)
    a_nodes=_sine_nodes(nodes) if smooth_a else nodes
    if logarithmic:
        length=log1p(max(s1,s2,s3))
        denominator=-expm1(-length)
        radial_nodes=tuple((exp(length*(x-1))*(-expm1(-length*x))/denominator,
                            weight*length*exp(length*(x-1))/denominator)
                           for x,weight in a_nodes)
    else:
        radial_nodes=a_nodes
    outer=[]
    for half in (0,1):
      for chart in (0,1):
        for a,wa in radial_nodes:
            for t,wt in nodes:
                if chart==0:
                    if logarithmic:
                        A=a;length=-log(A)
                        radial=exp(-(1-t)*length)
                        t=exp(-t*length)
                        wt*=length
                    else:
                        radial=a;A=radial*t
                else:
                    A=a;radial=A*t
                C=1-A
                u=radial/2 if half==0 else 1-radial/2
                # In the upper half u rounds to 1 for tiny radial. The exact
                # complementary factor is still representable as radial/2.
                vbase=1-radial/2 if half==0 else radial/2
                h0=.5*(1-radial/2)
                if chart==0:
                    d=t+C*h0
                    first=.5*radial*vbase*h0*(d*d+t*d+t*t)/d**3
                    second=.5*vbase*t*t/d**4
                else:
                    d=1+C*t*h0
                    first=.5*A*vbase*t*h0*(d*d+d+1)/d**3
                    second=.5*vbase/d**4
                for z,wz in nodes:
                    v=vbase*z;w=vbase*(1-z)
                    uv=u*v;uw=u*w;vw=v*w
                    b0=(.5*vbase*vbase if half==0 else u*radial/4)*z*(1-z)
                    extra=C*b0*s1
                    inner=[]
                    for W,norm0,norm1,weight in shapes:
                        G=uv*norm0+uw*norm1+vw*s1
                        N=t*G+extra if chart==0 else G+t*extra
                        delta=C*N/d
                        base=1+A*W
                        value=first*W*_ell(A*W)-second*N/base*_ell(delta/base)
                        inner.append(weight*value)
                    outer.append(wa*wt*wz*fsum(inner))
    return fsum(outer)


@lru_cache(maxsize=64)
def _numpy_normalized(s1,s2,s3,order,smooth_a=False):
    """Vectorize the same paired5D kernel over outer fractions and inner z."""
    import numpy as np
    nodes=_gauss_legendre(order)
    logarithmic=max(s1,s2,s3)>100
    shapes=np.asarray(_outer_shapes((s1,s2,s3),_sine_nodes(nodes) if logarithmic else nodes))
    W,norm0,norm1,shape_weight=(shapes[:,i][None,:] for i in range(4))
    z=np.asarray([x for x,_ in nodes])[:,None]
    wz=np.asarray([weight for _,weight in nodes])[:,None]
    a_nodes=_sine_nodes(nodes) if smooth_a else nodes
    if logarithmic:
        length=log1p(max(s1,s2,s3));denominator=-expm1(-length)
        radial_nodes=tuple((exp(length*(x-1))*(-expm1(-length*x))/denominator,
                           weight*length*exp(length*(x-1))/denominator)
                          for x,weight in a_nodes)
    else:radial_nodes=a_nodes
    def ell(value):
        return np.divide(np.log1p(value),value,out=np.ones_like(value),where=value!=0)
    totals=[]
    for half in (0,1):
      for chart in (0,1):
        chart_terms=[]
        for a,wa in radial_nodes:
          for t,wt in nodes:
            if chart==0:
                if logarithmic:
                    A=a;length=-log(A)
                    radial=exp(-(1-t)*length);t=exp(-t*length);wt*=length
                else:radial=a;A=radial*t
            else:A=a;radial=A*t
            C=1-A;u=radial/2 if half==0 else 1-radial/2
            vbase=1-radial/2 if half==0 else radial/2
            h0=.5*(1-radial/2)
            if chart==0:
                d=t+C*h0
                first=.5*radial*vbase*h0*(d*d+t*d+t*t)/d**3
                second=.5*vbase*t*t/d**4
            else:
                d=1+C*t*h0
                first=.5*A*vbase*t*h0*(d*d+d+1)/d**3
                second=.5*vbase/d**4
            v=vbase*z;w=vbase*(1-z)
            b0=(.5*vbase*vbase if half==0 else u*radial/4)*z*(1-z)
            G=u*v*norm0+u*w*norm1+v*w*s1
            extra=C*b0*s1
            N=t*G+extra if chart==0 else G+t*extra
            delta=C*N/d;base=1+A*W
            value=first*W*ell(A*W)-second*N/base*ell(delta/base)
            chart_terms.append(wa*wt*float(np.sum(wz*shape_weight*value)))
        totals.append(fsum(chart_terms))
    return fsum(totals)


def _backend(backend):
    if backend not in ('auto','python','numpy'):
        raise ValueError("backend must be 'auto', 'python' or 'numpy'")
    if backend!='python':
        try:import numpy
        except ImportError:
            if backend=='numpy':raise ImportError('The requested NumPy backend is not installed') from None
        else:return 'numpy'
    return 'python'


def _log_quadratic(v0,q1,q2,np=None):
    """Integral_0^1 log(1+v0+q1*y+q2*y*y)dy for a positive concave q.

    The roots factor q/q(0)=(1+plus*y)(1+minus*y). Opposite signs permit
    a product formula for the small root, avoiding subtraction cancellation.
    """
    b=q1/(1+v0);c=q2/(1+v0)
    if np is None:
        disc=hypot(b,2*sqrt(max(0.,-c)))
        if b>=0:
            plus=.5*b+.5*disc;minus=c/plus if plus else 0.
        else:
            minus=.5*b-.5*disc;plus=c/minus if minus else 0.
        if minus < -1-2e-14:raise ArithmeticError('Positive quadratic lost its legal root domain')
        minus=max(-1.,minus)
        def primitive(z):
            if abs(z)<1e-4:return z*(.5+z*(-1/6+z*(1/12+z*(-1/20+z/30))))
            return -1. if z==-1 else (1+1/z)*log1p(z)-1
        return fsum((log1p(v0),primitive(plus),primitive(minus)))
    disc=np.hypot(b,2*np.sqrt(np.maximum(0.,-c)))
    positive=.5*b+.5*disc;negative=.5*b-.5*disc
    plus=np.where(b>=0,positive,np.divide(c,negative,out=np.zeros_like(c),where=negative!=0))
    minus=np.where(b<0,negative,np.divide(c,positive,out=np.zeros_like(c),where=positive!=0))
    if np.any(minus < -1-2e-14):raise ArithmeticError('Positive quadratic lost its legal root domain')
    minus=np.maximum(-1.,minus)
    def primitive(z):
        small=z*(.5+z*(-1/6+z*(1/12+z*(-1/20+z/30))))
        with np.errstate(divide='ignore',invalid='ignore'):
            value=(1+1/z)*np.log1p(z)-1
        return np.where(abs(z)<1e-4,small,np.where(z==-1,-1.,value))
    return np.log1p(v0)+primitive(plus)+primitive(minus)


def _paired_y_average(invariants,nodes,A,C,u,vbase,radial,t,d,chart,half,a,z,np=None):
    """Stable fallback for an analytically integrated logarithm difference.

    This evaluates its already-paired ell kernel, avoiding lograw-logct at
    tiny relative corrections. The caller has a fixed outer fraction a and z.
    """
    s1,s2,s3=invariants;scale=max(invariants)
    p1x,p1y,p2x,p2y=normalized_gram_vectors(invariants)
    if np is not None:
        y=np.asarray([x for x,_ in nodes])[:,None]
        weights=np.asarray([w for _,w in nodes])[:,None]
        a=np.asarray(a)[None,:];z=np.asarray(z)[None,:]
    terms=[]
    for y,wy in ([(y,weights)] if np is not None else nodes):
        b=(1-a)*y;c=(1-a)*(1-y)
        W=a*b*s2+a*c*s1+b*c*s3
        norm0=scale*((c*p1x-b*p2x)**2+(c*p1y-b*p2y)**2)
        norm1=scale*(((c-1)*p1x-b*p2x)**2+((c-1)*p1y-b*p2y)**2)
        v=vbase*z;w=vbase*(1-z)
        G=u*v*norm0+u*w*norm1+v*w*s1
        b0=(.5*vbase*vbase if half==0 else u*radial/4)*z*(1-z)
        N=t*G+C*b0*s1 if chart==0 else G+t*C*b0*s1
        base=1+A*W;argument=C*N/(d*base)
        if np is None:terms.append(wy*N/base*_ell(argument))
        else:
            quotient=np.divide(np.log1p(argument),argument,out=np.ones_like(argument),where=argument!=0)
            return np.sum(wy*N/base*quotient,axis=0)
    return fsum(terms)


def _zero_outer_pair(scale,nodes,A,C,u,vbase,radial,t,d,chart,half,z,np=None):
    """Paired fallback after the exact s3=0 outer-fraction symmetry."""
    if np is not None:
        a=np.asarray([x for x,_ in nodes])[:,None]
        weights=np.asarray([weight for _,weight in nodes])[:,None]
        z=np.asarray(z)[None,:]
    terms=[]
    for a,weight in ([(a,weights)] if np is not None else nodes):
        v=vbase*z;w=vbase*(1-z)
        W=scale*a*(1-a)
        G=scale*(u*v*(1-a)**2+u*w*a*a+v*w)
        b0=(.5*vbase*vbase if half==0 else u*radial/4)*z*(1-z)
        N=t*G+C*b0*scale if chart==0 else G+t*C*b0*scale
        base=1+A*W;argument=C*N/(d*base)
        if np is None:terms.append(weight*N/base*_ell(argument))
        else:
            quotient=np.divide(np.log1p(argument),argument,out=np.ones_like(argument),where=argument!=0)
            return np.sum(weight*N/base*quotient,axis=0)
    return fsum(terms)


@lru_cache(maxsize=64)
def _analytic_normalized(s1,s2,s3,order,smooth_a,backend):
    """Four-dimensional high-scale rule: integrate outer Duffy y exactly."""
    np=None
    if backend=='numpy':import numpy as np
    nodes=_gauss_legendre(order);phase=_sine_nodes(nodes)
    zero_outer=s3==0
    z_nodes=phase if smooth_a==2 else nodes
    length=log1p(max(s1,s2,s3));denominator=-expm1(-length)
    radial_nodes=tuple((exp(length*(x-1))*(-expm1(-length*x))/denominator,
                       weight*length*exp(length*(x-1))/denominator)
                      for x,weight in (phase if smooth_a else nodes))
    if np is not None:
        a=np.asarray([x for x,_ in phase])[None,:]
        phase_weight=np.asarray([w for _,w in phase])[None,:]*(1-a)
        z=np.asarray([x for x,_ in z_nodes])[:,None]
        zweight=np.asarray([w for _,w in z_nodes])[:,None]
    result=[]
    for half in (0,1):
      for chart in (0,1):
       for A,wa in radial_nodes:
        C=1-A
        shapes=[]
        if zero_outer:
            shapes=[(None,None,.5,0.,A*s1,_log_quadratic(0.,A*s1,-A*s1,np))]
        for aa,wx in ([] if zero_outer else [(a,phase_weight)] if np is not None else phase):
            outerd=1-aa
            ct0=A*aa*outerd*s1
            ct1=A*(aa*outerd*(s2-s1)+outerd**2*s3)
            ctlog=_log_quadratic(ct0,ct1,-A*outerd**2*s3,np)
            shapes.append((aa,outerd,wx if np is not None else wx*outerd,ct0,ct1,ctlog))
        for y,wy in nodes:
            if chart==0:
                length=-log(A);radial=exp(-(1-y)*length);t=exp(-y*length);wt=wy*length
            else:radial=A*y;t=y;wt=wy
            u=radial/2 if half==0 else 1-radial/2
            vbase=1-radial/2 if half==0 else radial/2;h0=.5*(1-radial/2)
            if chart==0:
                d=t+C*h0
                first=.5*radial*vbase*h0*(d*d+t*d+t*t)/d**3
                second=.5*vbase*t*t/d**4
            else:
                d=1+C*t*h0
                first=.5*A*vbase*t*h0*(d*d+d+1)/d**3
                second=.5*vbase/d**4
            terms=[]
            for zz,wz in ([(z,zweight)] if np is not None else z_nodes):
             v=vbase*zz;w=vbase*(1-zz)
             for aa,outerd,weight,ct0,ct1,ctlog in shapes:
                if zero_outer:
                    g0=s1*(u*v+v*w);g1=-2*s1*u*v
                else:
                    g0=s1*(u*v*outerd**2+u*w*aa*aa+v*w)
                    g1=(s1+s3-s2)*u*outerd*(w*aa-v*outerd)
                b0=(.5*vbase*vbase if half==0 else u*radial/4)*zz*(1-zz)
                N0=t*g0+C*b0*s1 if chart==0 else g0+t*C*b0*s1
                factor=t if chart==0 else 1.
                rawlog=_log_quadratic(ct0+C*N0/d,ct1+C*factor/d*g1,
                                      -A*factor/d*(s1 if zero_outer else outerd**2*s3),np)
                difference=rawlog-ctlog
                if np is None:
                    if abs(difference)<1e-6*max(1.,abs(rawlog),abs(ctlog)):
                        paired=(_zero_outer_pair(s1,phase,A,C,u,vbase,radial,t,d,chart,half,zz)
                                if zero_outer else _paired_y_average((s1,s2,s3),nodes,A,C,u,vbase,radial,t,d,chart,half,aa,zz))
                    else:paired=d/C*difference
                    terms.append(wz*weight*(first*(ctlog/A)-second*paired))
                else:
                    paired=d/C*difference
                    mask=abs(difference)<1e-6*np.maximum(1.,np.maximum(abs(rawlog),abs(ctlog)))
                    if np.any(mask):
                        if zero_outer:
                            paired[mask]=_zero_outer_pair(s1,phase,A,C,u,vbase,radial,t,d,chart,half,zz[mask],np)
                        else:
                            aa_full,zz_full=np.broadcast_arrays(aa,zz)
                            paired[mask]=_paired_y_average((s1,s2,s3),nodes,A,C,u,vbase,radial,t,d,chart,half,
                                                           aa_full[mask],zz_full[mask],np)
                    terms.append(float(np.sum(wz*weight*(first*(ctlog/A)-second*paired))))
            result.append(wa*wt*fsum(terms))
    return fsum(result)


def _point(s1,s2,s3,m2):
    values=_invariants(s1,s2,s3);mass=_number(m2,'m2',True)
    normalized=tuple(s/mass for s in values)
    if any(not isfinite(v) for v in normalized):
        raise ValueError('si/m2 exceeds supported floating-point range')
    normalized=validate_three_point(normalized)
    # The two physical outer legs are interchangeable. Canonicalize them so
    # exchange symmetry is exact also at finite quadrature order, and either
    # zero outer leg reaches the exact three-dimensional specialization.
    normalized=(normalized[0],*sorted(normalized[1:],reverse=True))
    return values,mass,normalized


def _scale_point(s1,s2,s3,m2):
    values=_invariants(s1,s2,s3);mass=_number(m2,'m2',True)
    canonical=(values[0],*sorted(values[1:],reverse=True))
    return values,mass,canonical,_scale_log_ratio(max(values),mass)


def _order(order):
    if type(order)is not int or not 2<=order<=64:
        raise ValueError('order must be an integer between2and64')
    return order


def _evaluate(normalized,order,backend,smooth_a=None):
    if not any(normalized):return 0.
    if 100<max(normalized)<=1e12:
        if smooth_a is None:smooth_a=2
        result=_analytic_normalized(*normalized,order,smooth_a,backend)
    else:
        if smooth_a is None:smooth_a=False
        kernel=_numpy_normalized if backend=='numpy' else _python_normalized
        result=kernel(*normalized,order,smooth_a)
    if not isfinite(result):raise ValueError('The finite value exceeds floating-point range')
    return result


def nested_triangle_bph(s1=1,s2=1,s3=1,m2=1,*,order=12,backend='auto'):
    """Finite recursive BPH amplitude; s1 belongs to the inserted triangle.

    The per-coordinate maximum is64. Above scale100, the complete paired A
    integral is evaluated analytically with real dilogarithms. Its exact
    quadratic/linear log-scale polynomial is separated from a bounded
    residual, so ratios that overflow binary64 remain accessible. Smaller
    scales use four five-dimensional charts. Optional NumPy vectorization
    uses the identical kernel. With backend='auto'
    it is selected if installed, otherwise the standard-library path runs.
    High orders are substantially slower without NumPy. Compare refinements
    and use nested_triangle_diagnostics for high-scale coordinate crosschecks.
    This integral is dimensionless and depends on si/m2 only.
    """
    values,mass,canonical,log_ratio=_scale_point(s1,s2,s3,m2)
    order=_order(order);selected=_backend(backend)
    if log_ratio>log(100.):
        return exact_a_components(canonical,mass,order=order,backend=selected,power=4)['value']
    _,_,normalized=_point(*values,mass)
    return _evaluate(normalized,order,selected)


def nested_triangle_for_graph(graph_id,invariants=(1,1,1),m2=1,*,order=12,backend='auto'):
    """Map fixed-label atlas invariants to the standard physical routing."""
    if graph_id not in GRAPH_LEG_MAP:raise ValueError('Use L2E3G2, L2E3G3 or L2E3G5')
    values=_invariants(*invariants)
    return nested_triangle_bph(*(values[i] for i in GRAPH_LEG_MAP[graph_id]),m2,order=order,backend=backend)


def _exact_a_diagnostics(values,mass,canonical,orders,backend):
    routes={};components={}
    for name,power,trend in (('exact_A_cubic',3,6),('exact_A_quartic',4,8)):
        samples=[]
        for order in orders:
            start=monotonic()
            last=exact_a_components(canonical,mass,order=order,backend=backend,power=power)
            samples.append({'order':order,'value':last['value'],'residual':last['residual'],
                            'seconds_including_cache_reuse':monotonic()-start})
        components[name]=last
        extrapolations=[]
        for low,high in zip(samples,samples[1:]):
            ratio=(high['order']/low['order'])**trend
            extrapolations.append({'orders':[low['order'],high['order']],
                                   'value':high['value']+(high['value']-low['value'])/(ratio-1)})
        routes[name]={'quadrature':samples,'endpoint_power':power,
                      'last_refinement_difference':abs(samples[-1]['value']-samples[-2]['value']),
                      'richardson_power':trend,'richardson_estimates':extrapolations,
                      'last_extrapolated_refinement_difference':(abs(extrapolations[-1]['value']-extrapolations[-2]['value'])
                                                               if len(extrapolations)>1 else None)}
    cubic=routes['exact_A_cubic'];quartic=routes['exact_A_quartic']
    difference=abs(cubic['quadrature'][-1]['value']-quartic['quadrature'][-1]['value'])
    extrapolated_difference=abs(cubic['richardson_estimates'][-1]['value']-quartic['richardson_estimates'][-1]['value'])
    version=None
    if backend=='numpy':
        import numpy as np
        version=np.__version__
    production=components['exact_A_quartic'];scale=max(canonical)
    return {'invariants':list(values),'m2':mass,'orders':list(orders),'backend':backend,'numpy_version':version,
            'canonical_physical_invariants':list(canonical),'scaled_shape_invariants':[s/scale for s in canonical],
            'normalized_evaluation_invariants':None,'log_scale_ratio':production['log_ratio'],
            'integration_dimension':production['integration_dimension'],
            'analytic_reduction':'Complete paired A integral evaluated with real dilogarithms; exact log-scale polynomial separated before phase quadrature',
            'production_route':'exact_A_quartic','production_value':production['value'],
            'scale_decomposition':production,'routes':routes,
            'coordinate_route_difference':difference,'extrapolated_route_difference':extrapolated_difference,
            'coordinate_routes_agree_at_1e7':difference<=1e-7,
            'extrapolated_routes_agree_at_1e7':extrapolated_difference<=1e-7,
            'extrapolated_pair_average':(cubic['richardson_estimates'][-1]['value']+quartic['richardson_estimates'][-1]['value'])/2,
            'normalization':'d6k/pi3 d6q/pi3; no1/S or coupling',
            'scope':'Exact finite recursive BPH A integration and cubic/quartic endpoint maps. Refinement and extrapolation are numerical diagnostics, not rigorous quadrature error bounds. The separately proven finite-mass asymptotic remainder has a different scope.'}


def nested_triangle_diagnostics(s1=1,s2=1,s3=1,m2=1,*,orders=(24,32,40),backend='auto'):
    """Compare exact coordinate maps and observed refinement trends.

    The joint-smoothed route reuses nested_triangle_bph's high-scale cache.
    Raw-A and sine-square A smoothing give different node distributions for
    the same integral. A third route smooths A and inner-z jointly. Richardson
    extrapolation is a measured convergence diagnostic, not an error theorem.
    A minimum of two distinct increasing orders is required. In a stdlib-only
    environment use manageable orders; high orders are expensive even with
    the four-dimensional analytic reduction.
    """
    values,mass,canonical,log_ratio=_scale_point(s1,s2,s3,m2)
    try:orders=tuple(_order(n) for n in orders)
    except TypeError as exc:raise ValueError('Supply at least two increasing orders') from exc
    if len(orders)<2 or any(a>=b for a,b in zip(orders,orders[1:])):
        raise ValueError('Supply at least two distinct increasing orders')
    selected=_backend(backend)
    if log_ratio>log(100.):
        return _exact_a_diagnostics(values,mass,canonical,orders,selected)
    _,_,normalized=_point(*values,mass)
    routes={}
    highscale=100<max(normalized)<=1e12
    route_parameters=[('raw_A',False),('smoothed_A',True)]
    if highscale:route_parameters.append(('joint_smoothed',2))
    for name,smooth in route_parameters:
        samples=[]
        for order in orders:
            start=monotonic();value=_evaluate(normalized,order,selected,smooth)
            samples.append({'order':order,'value':value,'seconds_including_cache_reuse':monotonic()-start})
        extrapolations=[]
        for low,high in zip(samples,samples[1:]):
            ratio=(high['order']/low['order'])**4
            extrapolations.append({'orders':[low['order'],high['order']],
                                   'value':high['value']+(high['value']-low['value'])/(ratio-1)})
        routes[name]={'quadrature':samples,'last_refinement_difference':abs(samples[-1]['value']-samples[-2]['value']),
                      'richardson_power':4,'richardson_estimates':extrapolations,
                      'last_extrapolated_refinement_difference':(abs(extrapolations[-1]['value']-extrapolations[-2]['value'])
                                                               if len(extrapolations)>1 else None)}
    raw=routes['raw_A'];smooth=routes['smoothed_A']
    difference=abs(raw['quadrature'][-1]['value']-smooth['quadrature'][-1]['value'])
    extrapolated_difference=abs(raw['richardson_estimates'][-1]['value']-smooth['richardson_estimates'][-1]['value'])
    version=None
    if selected=='numpy':
        import numpy as np
        version=np.__version__
    return {'invariants':list(values),'m2':mass,'orders':list(orders),'backend':selected,'numpy_version':version,
            'normalized_evaluation_invariants':list(normalized),
            'integration_dimension':(3 if normalized[2]==0 else 4) if highscale else 5,
            'analytic_reduction':'Outer Duffy fraction integrated by positive concave log-quadratic primitive; tiny log differences use the paired y kernel' if 100<max(normalized)<=1e12 else None,
            'production_route':'joint_smoothed' if highscale else 'raw_A',
            'production_value':routes['joint_smoothed' if highscale else 'raw_A']['quadrature'][-1]['value'],
            'routes':routes,'coordinate_route_difference':difference,
            'extrapolated_route_difference':extrapolated_difference,
            'coordinate_routes_agree_at_1e7':difference<=1e-7,
            'extrapolated_routes_agree_at_1e7':extrapolated_difference<=1e-7,
            'extrapolated_pair_average':(raw['richardson_estimates'][-1]['value']+smooth['richardson_estimates'][-1]['value'])/2,
            'normalization':'d6k/pi3 d6q/pi3; no1/S or coupling',
            'scope':'Two coordinate maps of the finite paired BPH integral; n^-4 is an observed trend. Neither route agreement nor extrapolation is a rigorous quadrature-error enclosure.'}


def _upper_float(value):
    result=float(value)
    if not isfinite(result):return None
    return nextafter(result,inf) if Decimal.from_float(result)<value else result


def nested_triangle_bounds(s1=1,s2=1,s3=1,m2=1):
    """Proven parameter L1 / amplitude bounds with outward numeric enclosures.

    S=max(si), C=103/192+7*log(2)/288, L=log(1+S/(3*m2)).
    Integral|paired parameter density| <= min(C*S/m2,9L/4+5L²/8).
    The separate momentum-space absolute-convergence argument in the linked
    derivation identifies this value with the original recursive BPH integral.
    This bound is not a numerical quadrature error estimate.
    """
    values=_invariants(s1,s2,s3);m2=_number(m2,'m2',True)
    S=Decimal.from_float(max(values));M=Decimal.from_float(m2)
    with localcontext() as context:
        context.prec=max(100,len(S.as_tuple().digits)+5,len(M.as_tuple().digits)+5)
        context.rounding=ROUND_CEILING
        C=Decimal(103)/192+7*Decimal(2).ln().next_plus()/288
        L=(1+S/(3*M)).ln().next_plus() if S else Decimal(0)
        linear=C*S/M;global_bound=9*L/4+5*L*L/8
        bound=min(linear,global_bound)
        raw={'linear_constant':C,'log_scale':L,'linear_bound':linear,
             'global_log_bound':global_bound,'absolute_integral_bound':bound}
        result={k:_upper_float(v) for k,v in raw.items()}
        result['decimal_upper_bounds']={k:str(v) for k,v in raw.items()}
    return {'s1':values[0],'s2':values[1],'s3':values[2],'m2':m2,**result,
            'constant_exact':'103/192+7*log(2)/288',
            'formulas':{'linear':'C*max(si)/m2','global':'9*L/4+5*L^2/8, L=log(1+max(si)/(3*m2))'},
            'scope':'This nested triangle graph, paired parameter L1 and amplitude magnitude; not momentum L1 with the same constant',
            'not_claimed':['quadrature error bound','arbitrary-loop theorem','different renormalization scheme'],
            'numeric_rounding':'Upward Decimal arithmetic and outward binary64 conversion for the converted inputs'}


def nested_triangle_benchmark():
    """One standard routing and all three fixed-label placements at an asymmetric point."""
    records=[]
    for values in ((1.,1.,1.),(1.,2.,2.5)):
        for graph_id in GRAPH_LEG_MAP:
            mapping=GRAPH_LEG_MAP[graph_id]
            mapped=tuple(values[i] for i in mapping)
            quadrature=[{'order':n,'value':nested_triangle_bph(*mapped,order=n)} for n in (6,8,12)]
            difference=abs(quadrature[-1]['value']-quadrature[-2]['value'])
            records.append({'id':graph_id,'invariants':list(values),'standard_invariants':list(mapped),
                            'inserted_physical_leg':mapping[0]+1,'value':quadrature[-1]['value'],
                            'quadrature':quadrature,'last_refinement_difference':difference,
                            'bounds':nested_triangle_bounds(*mapped),
                            'checks':{'zero_subtraction_point':nested_triangle_bph(0,0,0)==0.,
                                      'quadrature_refinement_stable':difference<2e-7},
                            'normalization':'d6k/pi3 d6q/pi3; no1/S or coupling',
                            'limits':'Actual finite BPH value for this vertex-insertion graph. Refinement is diagnostic, not certified quadrature error.'})
    return records
