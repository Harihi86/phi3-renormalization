"""All ten two-loop graphs in one explicit massive zero-momentum BPH scheme."""
from functools import lru_cache
from decimal import Decimal
from math import isfinite, log
from time import perf_counter
from . import __version__
from .overlap import overlap_bph
from .overlap_bounds import overlap_bound
from .provenance import source_fingerprint
from .two_loop_registry import SCHEME, one_point_bph, two_loop_catalog
from .kinematics import kinematics_status, validate_three_point


@lru_cache(maxsize=1)
def _registry():
    return {row['id']: row for row in two_loop_catalog()}


def _number(value, name, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(name + ' 须为有限数')
    try:
        value = float(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(name + ' 须为有限数') from exc
    if not isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(name + (' 须为正的有限数' if positive else ' 须为非负有限数'))
    return value


def _route(row, invariants, mass):
    family = row['family']
    if family == 'overlap':
        return lambda n: overlap_bph(invariants[0], mass, order=n), overlap_bound(invariants[0], mass)
    if family == 'nested_bubble_two_point':
        from .nested_bubble import nested_bubble_two_point, nested_bubble_two_point_bounds
        return (lambda n: nested_bubble_two_point(invariants[0], mass, order=n),
                nested_bubble_two_point_bounds(invariants[0], mass))
    if family == 'primitive_three_point':
        from .primitive import primitive_bph, primitive_bound
        return lambda n: primitive_bph(invariants, mass, order=n), primitive_bound(invariants, mass)
    if family == 'nested_triangle':
        from .nested_triangle import nested_triangle_bph, nested_triangle_bounds
        focus = row['focus_label'] - 1
        other = [j for j in range(3) if j != focus]
        routed = (invariants[focus], invariants[other[0]], invariants[other[1]])
        return (lambda n: nested_triangle_bph(*routed, m2=mass, order=n),
                nested_triangle_bounds(*routed, m2=mass))
    if family == 'nested_bubble_three_point':
        from .nested_bubble import nested_bubble_three_point, nested_bubble_three_point_bounds
        focus = row['focus_label'] - 1
        other = [j for j in range(3) if j != focus]
        routed = (invariants[other[0]], invariants[other[1]], invariants[focus])
        return (lambda n: nested_bubble_three_point(*routed, m2=mass, order=n),
                nested_bubble_three_point_bounds(*routed, m2=mass))
    raise ValueError('尚无此图的数值路线')


def evaluate_two_loop(graph='L2E2G1', invariants=None, m2=1, order=None, *, refine=True):
    """Return finite amplitude, separate weight, analytic bounds and diagnostics.

    Two-point invariants are (p_E^2,); three-point invariants are always
    (p1_E^2,p2_E^2,p3_E^2), in the graph's fixed physical external label order.
    Family-specific integration orders are reported explicitly. No physical
    vertex sign, coupling, ordinary (2*pi)^(-6) measure, or MS conversion is
    applied implicitly.
    """
    if not isinstance(graph, str) or graph not in _registry():
        raise ValueError('请选择图集中的有效两圈图编号')
    row = _registry()[graph]
    mass = _number(m2, '共同质量平方', positive=True)
    family = row['family']
    legs = len(row['external_labels'])
    expected = 0 if legs == 1 else 1 if legs == 2 else 3
    if invariants is None:
        invariants = (1,) * expected
    if isinstance(invariants, (str, bytes)):
        raise ValueError('请按固定外腿标签输入不变量序列')
    try:
        invariants = tuple(_number(x, '外动量不变量') for x in invariants)
    except TypeError as exc:
        raise ValueError('请按固定外腿标签输入不变量序列') from exc
    if len(invariants) != expected:
        raise ValueError(f'此图需要 {expected} 个外动量不变量')
    requested_invariants = invariants
    gram_status = kinematics_status(invariants) if legs == 3 else None
    if legs == 3:
        invariants = validate_three_point(invariants)
    is_bubble = family.startswith('nested_bubble')
    default_order = 64 if is_bubble else 20 if family == 'nested_triangle' else 16
    selected_order = default_order if order is None else order
    max_order = 256 if is_bubble else 64 if family == 'nested_triangle' else 32
    if type(selected_order) is not int or not 8 <= selected_order <= max_order:
        raise ValueError(f'此图积分阶数须为 8 到 {max_order} 的整数')
    start = perf_counter()
    if family == 'one_point_local':
        exact = one_point_bph(mass)
        value, previous, bounds = 0.0, 0.0, {'absolute_integral_bound': 0.0,
                                                       'scope': exact['reason']}
        previous_order = selected_order
        method = 'Exact R-operation identity before integration'
    else:
        evaluator, bounds = _route(row, invariants, mass)
        value = evaluator(selected_order)
        previous_order = max(4, selected_order // 2 if is_bubble else selected_order - 4)
        previous = evaluator(previous_order) if refine else None
        method = family
    if not isfinite(value):
        raise ValueError('计算值超出浮点数支持范围，请调整尺度或精度')
    difference = abs(value - previous) if previous is not None else None
    tolerance = 1e-12 + 1e-9 * abs(value)
    precision_diagnostics = None
    external_scale=max(invariants,default=0.)
    log_scale_ratio=log(external_scale)-log(mass) if external_scale else None
    if family == 'nested_triangle' and refine and selected_order >= 24 and log_scale_ratio is not None and log_scale_ratio > log(100):
        from .nested_triangle import nested_triangle_diagnostics
        focus = row['focus_label'] - 1
        other = [j for j in range(3) if j != focus]
        routed = (invariants[focus], invariants[other[0]], invariants[other[1]])
        precision_diagnostics = nested_triangle_diagnostics(*routed, m2=mass,
                                   orders=(max(8, selected_order-8), selected_order))
    precision_scope = None
    extreme_scale_method = None
    if log_scale_ratio is not None and log_scale_ratio > log(100):
        extreme_scale_method = {'nested_triangle':'Exact analytic A integration; large logarithms separated before quadrature',
            'overlap':'Analytic counterterm moments and exact finite log-domain corrections',
            'nested_bubble_two_point':'Positive log-radius BPH integration resolving both physical scales',
            'nested_bubble_three_point':'Log-radius and endpoint-smoothed angular BPH integration',
            'primitive_three_point':'Endpoint smoothing with exact positive period normalization'}.get(family)
    enclosing_bound = bounds.get('absolute_integral_bound')
    if enclosing_bound is not None:
        consistent_bound = abs(value) <= enclosing_bound
    else:
        decimal_bound = bounds.get('decimal_upper_bounds', {}).get('absolute_integral_bound')
        consistent_bound = decimal_bound is not None and Decimal.from_float(abs(value)) <= Decimal(decimal_bound)
    return {'software_version': __version__, 'source_fingerprint': source_fingerprint(),
            'graph': graph, 'family': family, 'focus_label': row['focus_label'],
            'external_invariants': list(invariants), 'm2': mass,
            'requested_external_invariants': list(requested_invariants),
            'kinematics_status': gram_status,
            'order': selected_order, 'value': value,
            'order_limits': {'min': 8, 'max': max_order, 'default': default_order},
            'symmetry_weight': '1' if row['symmetry'] == 1 else f"1/{row['symmetry']}",
            'weighted_value': value / row['symmetry'],
            'vertex_factor': f'(-g)^{legs+2}, kept separate',
            'mass_dimension': 6-2*legs,
            'previous_order': previous_order if refine else None,
            'previous_value': previous if refine else None,
            'refinement_difference': difference if refine else None,
            'refinement_tolerance': tolerance,
            'refinement_stable': difference <= tolerance if refine else None,
            'checks': {'finite_value': True, 'within_analytic_amplitude_bound': consistent_bound,
                       'quadrature_refinement_stable': difference <= tolerance if refine else None},
            'precision_diagnostics': precision_diagnostics, 'precision_scope': precision_scope,
            'extreme_scale_method': extreme_scale_method, 'log_scale_ratio':log_scale_ratio,
            'seconds': perf_counter() - start,
            'analytic_bounds': bounds, 'method': method,
            'scheme': SCHEME,
            'normalization': 'd^6k/pi^3 d^6q/pi^3; no coupling, vertex sign or (4*pi)^(-6); graph 1/S appears only in weighted_value',
            'numerical_status': 'Quadrature order differences are stability diagnostics, not rigorous error enclosures.',
            'one_point_note': one_point_bph(mass)['not_claimed'] if family == 'one_point_local' else None}


def two_loop_benchmarks(graph_ids=None):
    selected = set(_registry()) if graph_ids is None else set(graph_ids)
    if selected - set(_registry()):
        raise ValueError('Unknown two-loop benchmark graph')
    return [evaluate_two_loop(ident) for ident in _registry() if ident in selected]
