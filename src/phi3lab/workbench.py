"""User-facing finite overlap calculations with separate numerical/analytic data."""
from math import isfinite
from time import perf_counter
from . import __version__
from .overlap import overlap_bph
from .overlap_bounds import overlap_bound
from .provenance import source_fingerprint


def evaluate_overlap(s=1, m2=1, order=16, *, refine=True):
    if type(order) is not int or not 8 <= order <= 32:
        raise ValueError('积分阶数须为 8 到 32 的整数')
    try:s,m2=float(s),float(m2)
    except (TypeError,ValueError,OverflowError) as exc:
        raise ValueError('请输入有效的外动量平方与质量平方') from exc
    if not isfinite(s) or s<0:raise ValueError('外动量平方须为非负有限数')
    if not isfinite(m2) or m2<=0:raise ValueError('共同质量平方须为正的有限数')
    start=perf_counter()
    value=overlap_bph(s,m2,order=order)
    bounds=overlap_bound(s,m2)
    previous_order=max(4,order-4)
    previous=overlap_bph(s,m2,order=previous_order) if refine else None
    difference=abs(value-previous) if refine else None
    tolerance=1e-12+1e-9*abs(value)
    return {
        'software_version':__version__,'source_fingerprint':source_fingerprint(),
        'graph':'L2E2G1','s':float(s),'m2':float(m2),'order':order,
        'value':value,'symmetry_weight':'1/2','weighted_value':value/2,
        'previous_order':previous_order if refine else None,
        'previous_value':previous,'refinement_difference':difference,
        'refinement_tolerance':tolerance,
        'refinement_stable':difference<=tolerance if refine else None,
        'seconds':perf_counter()-start,
        'analytic_bounds':bounds,
        'normalization':'d^6k/pi^3 d^6q/pi^3; no coupling or (4*pi)^(-6); 1/S is applied only in weighted_value',
        'BPH_recursion':'(1-t_p^2)[I-D(k)^3 D(q)D(p-q)-D(q)^3 D(k)D(p-k)]',
        'numerical_status':'Order differences diagnose stability; no rigorous quadrature error enclosure is claimed.',
        'proof_scope':'This particular massive Euclidean two-loop graph with zero-momentum Taylor conditions.',
    }
