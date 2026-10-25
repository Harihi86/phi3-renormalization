"""Executable topological checks and original SVG output; no amplitude claim."""

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from .graphs import Graph, degrees, divergent_subgraphs, generate, symmetry, validate_scope, weight


def describe(graph: Graph):
    validate_scope(graph)
    internal, vertices, external = len(graph.edges), graph.n, len(graph.external)
    return {
        **graph.to_dict(),
        "loops": graph.loops,
        "degrees": list(degrees(graph)),
        "omega_d6": 6 * graph.loops - 2 * internal,
        "trivalence_identity": 3 * vertices == 2 * internal + external,
        "euler_identity": graph.loops == internal - vertices + 1,
        "symmetry_order_fixed_external_labels": symmetry(graph),
        "combinatorial_weight": str(weight(graph)),
        "proper_divergent_subgraphs_d6": divergent_subgraphs(graph),
    }


def render_svg(graph: Graph):
    validate_scope(graph)
    n = graph.n
    if n == 1:
        positions = [(250.0, 145.0)]
    elif n == 2:
        positions = [(180.0, 145.0), (320.0, 145.0)]
    else:
        positions = [
            (
                250 + 60 * math.cos(-math.pi / 2 + 2 * math.pi * i / n),
                145 + 60 * math.sin(-math.pi / 2 + 2 * math.pi * i / n),
            )
            for i in range(n)
        ]
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="500" height="330" viewBox="0 0 500 330">',
        '<rect width="500" height="330" fill="#fcfaf6"/>',
        '<g fill="none" stroke="#204650" stroke-width="2.3">',
    ]
    counts, seen = Counter(graph.edges), Counter()
    for a, b in graph.edges:
        x1, y1 = positions[a]
        x2, y2 = positions[b]
        index = seen[(a, b)]
        seen[(a, b)] += 1
        if a == b:
            pieces.append(
                f'<path d="M {x1} {y1} C {x1 - 65} {y1 - 90}, {x1 + 65} {y1 - 90}, {x1} {y1}"/>'
            )
        else:
            length = math.hypot(x2 - x1, y2 - y1)
            offset = 55 * (index - (counts[(a, b)] - 1) / 2)
            cx, cy = (
                (x1 + x2) / 2 - offset * (y2 - y1) / length,
                (y1 + y2) / 2 + offset * (x2 - x1) / length,
            )
            pieces.append(f'<path d="M {x1} {y1} Q {cx} {cy}, {x2} {y2}"/>')
    pieces.append("</g>")
    for label, vertex in graph.external:
        x, y = positions[vertex]
        dx, dy = x - 250, y - 145
        length = math.hypot(dx, dy) or 1
        if n == 1:
            dx, dy, length = 0, 1, 1
        ex, ey = x + 55 * dx / length, y + 55 * dy / length
        label_x = ex + 7 if dx >= 0 else ex - 7
        anchor = "start" if dx >= 0 else "end"
        pieces.extend(
            [
                f'<path d="M {x} {y} L {ex} {ey}" stroke="#ba5a34" stroke-width="2"/>',
                f'<text x="{label_x}" y="{ey + 4}" text-anchor="{anchor}" font-family="sans-serif" font-size="13">p{label}</text>',
            ]
        )
    for i, (x, y) in enumerate(positions):
        pieces.extend(
            [
                f'<circle cx="{x}" cy="{y}" r="4.5" fill="#204650"/>',
                f'<text x="{x + 8}" y="{y + 15}" font-family="sans-serif" font-size="11">v{i}</text>',
            ]
        )
    pieces.append(
        f'<text x="250" y="303" text-anchor="middle" font-family="sans-serif" font-size="14">'
        f"L={graph.loops}, E={len(graph.external)}, weight={weight(graph)} (fixed external labels)</text></svg>"
    )
    return "\n".join(pieces)


def main():
    parser = argparse.ArgumentParser(description="Bounded bare phi-cubed multigraph demonstration")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input", type=Path, help="Optional graph JSON: vertices, edges, external")
    args = parser.parse_args()
    if args.input:
        graphs = [Graph.from_dict(json.loads(args.input.read_text(encoding="utf-8")))]
    else:
        graphs = [generate(1, 2)[0], generate(1, 3)[0]]
    reports = [describe(graph) for graph in graphs]
    args.output.mkdir(parents=True, exist_ok=False)
    for index, graph in enumerate(graphs):
        (args.output / f"graph-{index + 1}.svg").write_text(render_svg(graph), encoding="utf-8")
    report = {
        "scope": "bare amputated scalar graphs; fixed labelled external legs; d=6 power counting",
        "amplitude_calculation": "NOT_RUN: this bare-graph demo does not evaluate amplitudes; use phi3-ct for counterterms",
        "examples": reports,
        "bounded_two_loop_three_point_count": len(generate(2, 3)),
    }
    (args.output / "graphs.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "examples": len(reports),
                "two_loop_three_point_graphs": len(generate(2, 3)),
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
