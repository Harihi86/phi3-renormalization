from fractions import Fraction
from itertools import permutations

import pytest

from phi3_public.demo import describe, render_svg
from phi3_public.graphs import (
    Graph,
    canonical,
    degrees,
    divergent_subgraphs,
    generate,
    is_1pi,
    symmetry,
    validate_scope,
    weight,
)


def test_known_bubble_and_triangle_weights():
    bubble = Graph(2, ((0, 1), (0, 1)), ((1, 0), (2, 1)))
    triangle = Graph(3, ((0, 1), (1, 2), (0, 2)), ((1, 0), (2, 1), (3, 2)))
    assert weight(bubble) == Fraction(1, 2)
    assert weight(triangle) == Fraction(1, 1)
    assert describe(bubble)["omega_d6"] == 2
    assert describe(triangle)["omega_d6"] == 0
    assert not divergent_subgraphs(triangle)


@pytest.mark.parametrize("loops,legs", [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)])
def test_enumeration_obeys_independent_graph_identities(loops, legs):
    graphs = generate(loops, legs)
    assert graphs
    assert len({canonical(g)[0] for g in graphs}) == len(graphs)
    for g in graphs:
        validate_scope(g)
        assert degrees(g) == (3,) * g.n
        assert 3 * g.n == 2 * len(g.edges) + legs
        assert g.loops == len(g.edges) - g.n + 1 == loops
        assert g.n == 2 * loops + legs - 2
        assert 6 * g.loops - 2 * len(g.edges) == 6 - 2 * legs
        assert is_1pi(range(g.n), g.edges)


def test_seven_fixed_label_two_loop_three_point_topologies():
    assert len(generate(2, 3)) == 7


def test_all_vertex_relabellings_preserve_identity_and_weight():
    g = generate(2, 3)[0]
    for p in permutations(range(g.n)):
        other = Graph(
            g.n,
            tuple((p[a], p[b]) for a, b in g.edges),
            tuple((label, p[v]) for label, v in g.external),
        )
        assert canonical(g)[0] == canonical(other)[0]
        assert symmetry(g) == symmetry(other)


def test_external_labels_are_preserved():
    multiple_label_orbits = False
    for g in generate(2, 3):
        keys = set()
        for labels in permutations((1, 2, 3)):
            other = Graph(
                g.n, g.edges, tuple((labels[i], v) for i, (_, v) in enumerate(g.external))
            )
            keys.add(canonical(other)[0])
        assert len(keys) * symmetry(g, fixed=False) == 6 * symmetry(g, fixed=True)
        multiple_label_orbits |= len(keys) > 1
    assert multiple_label_orbits


@pytest.mark.parametrize(
    "args",
    [(True, (), ()), (2, ((0, 2),), ()), (2, ((0, 1.2),), ()), (1, ((0, 0),), ((1, 0), (1, 0)))],
)
def test_invalid_discrete_inputs_rejected(args):
    with pytest.raises(ValueError):
        Graph(*args)


def test_tree_is_not_1pi_and_unsupported_generation_rejected():
    assert not is_1pi(range(2), ((0, 1),))
    with pytest.raises(ValueError):
        generate(3, 3)


def test_svg_contains_actual_parallel_edges():
    g = generate(1, 2)[0]
    svg = render_svg(g)
    assert svg.count('<path d="M') == len(g.edges) + len(g.external)
    assert "weight=1/2" in svg


def test_subgraph_boundary_half_edges_counted_separately():
    for g in generate(2, 3):
        for sub in divergent_subgraphs(g):
            assert sub["omega"] == 6 * sub["loops"] - 2 * len(sub["edge_ids"])
            assert len(sub["boundary_ports"]) == sub["external_legs"]
            assert sub["omega"] == 6 - 2 * sub["external_legs"]
