"""Topology-derived coverage for the complete supported two-loop atlas.

This registry does not turn structural checks into integration claims. The
one-point result is the exact R-operation identity in the declared full
zero-momentum Taylor scheme, not a claim about the bare tadpole.
"""
from .graphs import generate, divergent_subgraphs, symmetry


SCHEME = {
    'method': 'Recursive BPH',
    'dimension': 6,
    'mass': 'One common positive mass squared',
    'overall_orders': {'1': 4, '2': 2, '3': 0},
    'subgraph_orders': {'bubble': 2, 'triangle': 0},
    'subtraction_point': 'All independent physical/subgraph external momenta zero',
    'measure': 'd^6k/pi^3 d^6q/pi^3',
    'excludes': ['coupling', 'symmetry weight', '(4*pi)^(-6)'],
    'one_point_condition': 'Full Taylor subtraction fixes the renormalized one-point graph to zero; it does not assert a vanishing bare integral or independently prove a vacuum expectation value condition.'
}


def classify_two_loop(g):
    """Classify topology and derive the physical label that fixes its mapping."""
    if g.loops != 2:
        raise ValueError('Exactly two loops are required')
    legs = len(g.external)
    subs = divergent_subgraphs(g, 6)
    if legs == 1:
        return {'family': 'one_point_local', 'focus_label': None,
                'identity': '(1-t_G^4) Rbar I = 0: there is no independent physical external momentum'}
    if legs == 2:
        if len(subs) == 2 and all(s['external_legs'] == 3 for s in subs):
            return {'family': 'overlap', 'focus_label': None}
        if len(subs) == 1 and subs[0]['external_legs'] == 2:
            return {'family': 'nested_bubble_two_point', 'focus_label': None}
    if legs == 3:
        if not subs:
            return {'family': 'primitive_three_point', 'focus_label': None}
        if len(subs) == 1:
            sub = subs[0]
            vertices = set(sub['vertices'])
            if sub['external_legs'] == 3:
                labels = [label for label, v in g.external if v in vertices]
                if len(labels) == 1:
                    return {'family': 'nested_triangle', 'focus_label': labels[0]}
            if sub['external_legs'] == 2:
                neighbors = {v for a, b in g.edges for v, w in ((a, b), (b, a))
                             if v not in vertices and w in vertices}
                labels = [label for label, v in g.external if v not in neighbors]
                if len(labels) == 1:
                    return {'family': 'nested_bubble_three_point', 'focus_label': labels[0]}
    raise ValueError('No supported two-loop topology matches this graph')


def two_loop_catalog():
    rows = []
    for legs in (1, 2, 3):
        for index, graph in enumerate(generate(2, legs), 1):
            rows.append({'id': f'L2E{legs}G{index}', **classify_two_loop(graph),
                         'external_labels': [label for label, _ in graph.external],
                         'symmetry': symmetry(graph),
                         'proper_subtractions': [
                             {'edge_ids': s['edge_ids'], 'order': s['omega']}
                             for s in divergent_subgraphs(graph, 6)]})
    return rows


def one_point_bph(m2=1):
    """Exact zero obtained before integrating any divergent forest term."""
    import math
    if isinstance(m2, bool):
        raise ValueError('Positive finite mass squared required')
    try:
        m2 = float(m2)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError('Positive finite mass squared required') from exc
    if not math.isfinite(m2) or m2 <= 0:
        raise ValueError('Positive finite mass squared required')
    return {'graph': 'L2E1G1', 'm2': m2, 'value': 0.0,
            'absolute_integral_bound': 0.0, 'status': 'exact_before_integration',
            'recursion': '(1-t_G^4) Rbar I = 0',
            'reason': 'Momentum conservation fixes the sole physical external momentum to zero. The prepared integrand has no independent external variable, so t_G^4 is the identity.',
            'local_counterterm': '-Rbar J_eta(0), defined with a common regulator before removal',
            'scheme': SCHEME,
            'not_claimed': 'The bare tadpole is not zero; the graphwise condition is part of this full Taylor subtraction scheme, not an independent assertion about the physical vacuum.'}
