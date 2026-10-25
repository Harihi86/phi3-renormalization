"""Exact domain and error accounting, conditional on independently joined tiles.

This module never issues a graph certificate. The caller must also validate
original graph identities, all saved native records, analytic proofs, source
hashes and actual successful exits. Partial domains cannot enter close().
"""
from fractions import Fraction as Q
from math import comb

CAP=1000000
XL,XR,YH=-84,37,72
ABSOLUTE_BUDGET=Q(1,10**8)
SIDES={'aligned':'L4E3G10','opposite':'L4E3G36'}


def bounded(value):
    assert type(value) in (int,Q),'Only exact integer or Fraction inputs are allowed'
    value=Q(value)
    assert max(value.numerator.bit_length(),value.denominator.bit_length())<=CAP
    return value


def policy(cfg):
    side=cfg['side'];assert side in SIDES and cfg['graph']==SIDES[side]
    assert cfg['order']==20,'This closure is bound to the actual N=20 campaigns'
    assert tuple(cfg['x_log2_domain'])==(XL,XR)
    assert tuple(map(tuple,cfg['w_charts']))==((0,1),(0,YH))
    assert tuple(cfg['original_radial_tail_domain'])==(-10,72)
    assert cfg['total_coordinate_area']==(XR-XL)*(YH+1)==8833
    assert cfg['absolute_finite_integral_budget']==ABSOLUTE_BUDGET
    assert cfg['bits']==512 and cfg['dyadic_depth_cap']==24 and cfg['native_call_cap']==500000
    assert tuple(cfg['gauss_orders'])==tuple(range(8,65,4)) and cfg['gauss_order']==64
    assert cfg['finite_domain_complete'] is True and cfg['audit_tiles'] is None
    assert cfg['whole_graph_certificates_added']==0 and cfg['requires_final_independent_audit'] is True
    return side


def keys(cfg,observed):
    policy(cfg);seen=set()
    for item in observed:
        assert len(item)==3 and all(type(x) is int for x in item)
        chart,x,y=item
        assert chart in (0,1) and XL<=x<XR and 0<=y<(1 if chart==0 else YH)
        key=tuple(item);assert key not in seen,'Duplicated finite-domain tile'
        seen.add(key)
    # Independent chart cardinalities plus coordinate bounds prove the full
    # Cartesian grid without trusting a producer's total counter or ordering.
    assert sum(c==0 for c,x,y in seen)==XR-XL,'Incomplete first chart'
    assert sum(c==1 for c,x,y in seen)==(XR-XL)*YH,'Incomplete second chart'
    assert len(seen)==8833
    return seen


def errors():
    n=20;epsilon=Q(1,2**10);radius=Q(2**72);split=radius/4;log=Q(71)
    count=sum(comb(2*j+7,7) for j in range(1,n+1));major=Q(count,4)
    lower=Q(101,5000)*major*epsilon**5
    p2=log**2+2*log+2;p3=log**3+3*log**2+6*log+6
    upper=2*major*((p2+4*p3)/(72*split)
        +16/radius*(Q(1,50400)+log**3/216)
        +40/radius**2*(Q(1,50400)+split*log**2/72))
    series=(Q(41,9)+Q(71,3628800))*(Q(1,4*(n+1))+Q(10,3)+5*n
        +Q(16,3)*comb(n,2)+Q(8,3)*comb(n,3))/4**n
    # The ordered rectangle contains the original [epsilon,R]^2 after
    # exchange symmetry. Thus its omitted region is bounded by the already
    # proved absolute radial tails; no second tail subtraction is permitted.
    assert Q(2)**XL<=2*epsilon/(radius+1)
    assert (Q(2)**XR)**2>=4*radius
    assert Q(2)**YH>=radius
    return {'lower_radial_error':bounded(lower),'upper_radial_error':bounded(upper),
        'radial_tail_error':bounded(lower+upper),'external_series_error':bounded(series),
        'all_nonquadrature_errors':bounded(lower+upper+series),
        'ordered_rectangle_contains_original_radial_square':True}


def enclosure(finite,terms):
    expected=errors();assert terms==expected,'Changed or omitted analytic remainder'
    lo,hi=map(bounded,finite);assert lo<=hi
    assert (hi-lo)/2<=ABSOLUTE_BUDGET,'Finite-domain allowance exceeded'
    extra=expected['all_nonquadrature_errors']
    interval=(bounded(lo-extra),bounded(hi+extra));mid=bounded(sum(interval)/2)
    rad=bounded((interval[1]-interval[0])/2)
    target=bounded((Q(1,10**8)+Q(3,100)*abs(mid))/16)
    assert rad<target,'The original strict whole-integral gate is not met'
    return {'finite_interval':(lo,hi),'certified_interval_if_all_evidence_verified':interval,
        'radius':rad,'strict_target':target,'strict_gate_passed':True,**expected,
        'whole_graph_certificate':False,'conditional_arithmetic_only':True}


def close(cfg,joined):
    side=policy(cfg);keys(cfg,[r['tile'] for r in joined]);lo=hi=Q();leaves=nodes=regions=boxes=0
    for stem in ('production','central','geometry'):
        assert len({r[stem] for r in joined})==len(joined),'A proof file was reused for distinct tiles'
    for row in joined:
        assert row['side']==side and row['graph']==cfg['graph']
        assert row['all_saved_batches_and_geometry_partitions_replayed'] is True
        assert row['finite_domain_complete'] is False and row['whole_graph_certificate'] is False
        a,b=map(bounded,row['interval']);assert a<=b
        assert (b-a)/2<=ABSOLUTE_BUDGET/cfg['total_coordinate_area']
        assert type(row['leaves']) is int and row['leaves']>0
        assert all(type(row[k]) is int and row[k]>=0 for k in ('independent_central_nodes','geometry_regions','geometry_subboxes'))
        for stem in ('production','central','geometry'):
            assert isinstance(row[stem],str) and row[stem]
            assert len(row[stem+'_sha256'])==64 and all(c in '0123456789abcdef' for c in row[stem+'_sha256'])
        lo=bounded(lo+a);hi=bounded(hi+b);leaves+=row['leaves'];nodes+=row['independent_central_nodes']
        regions+=row['geometry_regions'];boxes+=row['geometry_subboxes']
    return {'graph':cfg['graph'],'side':side,'tiles':len(joined),'leaves':leaves,
        'independent_central_nodes':nodes,'geometry_regions':regions,'geometry_subboxes':boxes,
        **enclosure((lo,hi),errors())}
