"""Small exact multigraph engine. No third-party dependencies or floating weights.

Edges have stable list indices: subgraphs are EDGE sets, never induced vertex sets.
External half-edges are labelled and amputated. Parallel edges and self-loops remain.
The exhaustive algorithms deliberately target <= 5 vertices (two-loop E <= 3).
"""

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, prod


@dataclass(frozen=True)
class Graph:
    n: int
    edges: tuple[tuple[int, int], ...]
    external: tuple[tuple[int, int], ...]  # (fixed leg label, attached vertex)
    kinds: tuple[str, ...] = ()

    def __post_init__(self):
        # A graph input is exact discrete data. In particular, bool is an int
        # subclass and int(1.9) silently truncates: neither is an accepted ID.
        if type(self.n) is not int or self.n < 1:
            raise ValueError("vertices must be a positive integer (not a boolean or decimal)")

        def integer_pairs(value, name):
            if not isinstance(value, (list, tuple)):
                raise ValueError(f"{name} must be an array of integer pairs")
            result = []
            for pair in value:
                if (
                    not isinstance(pair, (list, tuple))
                    or len(pair) != 2
                    or any(type(v) is not int for v in pair)
                ):
                    raise ValueError(f"Each {name} entry must contain exactly two integers")
                result.append(tuple(pair))
            return tuple(result)

        edges = integer_pairs(self.edges, "edges")
        external = integer_pairs(self.external, "external")
        if any(not (0 <= a < self.n and 0 <= b < self.n) for a, b in edges):
            raise ValueError("Edge endpoints must be vertex indices from 0 to vertices-1")
        if any(not 0 <= v < self.n for _, v in external):
            raise ValueError("External attachments must be vertex indices from 0 to vertices-1")
        if any(label <= 0 for label, _ in external):
            raise ValueError("External labels must be positive integers")
        if len({label for label, _ in external}) != len(external):
            raise ValueError("External labels must be unique")
        if not isinstance(self.kinds, (list, tuple)) or any(
            not isinstance(kind, str) or not kind.strip() for kind in self.kinds
        ):
            raise ValueError("vertex_kinds must be an array of nonempty strings")
        if self.kinds and len(self.kinds) != self.n:
            raise ValueError("One kind is required per vertex")
        object.__setattr__(self, "edges", tuple(tuple(sorted(e)) for e in edges))
        object.__setattr__(self, "external", tuple(sorted(external)))
        object.__setattr__(self, "kinds", tuple(self.kinds))

    @property
    def vertex_kinds(self):
        return self.kinds or ("phi3",) * self.n

    @property
    def loops(self):
        return len(self.edges) - self.n + len(components(range(self.n), self.edges))

    def to_dict(self):
        return {
            "vertices": self.n,
            "edges": [list(e) for e in self.edges],
            "external": [list(e) for e in self.external],
            "vertex_kinds": list(self.vertex_kinds),
        }

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, Mapping):
            raise ValueError("Graph JSON must be an object")
        missing = [key for key in ("vertices", "edges", "external") if key not in data]
        if missing:
            raise ValueError("Missing graph field(s): " + ", ".join(missing))
        return cls(data["vertices"], data["edges"], data["external"], data.get("vertex_kinds", ()))


def components(vertices, edges):
    todo = set(vertices)
    result = []
    adjacency = {v: set() for v in todo}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    while todo:
        seed = min(todo)
        seen = {seed}
        stack = [seed]
        while stack:
            for b in adjacency[stack.pop()] - seen:
                seen.add(b)
                stack.append(b)
        result.append(frozenset(seen))
        todo -= seen
    return tuple(result)


def connected(vertices, edges):
    vertices = tuple(vertices)
    return bool(vertices) and len(components(vertices, edges)) == 1


def is_1pi(vertices, edges):
    vertices = tuple(vertices)
    edges = tuple(edges)
    return connected(vertices, edges) and all(
        connected(vertices, edges[:i] + edges[i + 1 :]) for i in range(len(edges))
    )


def degrees(g):
    c = Counter(v for e in g.edges for v in e)
    c.update(v for _, v in g.external)
    return tuple(c[i] for i in range(g.n))


def validate_scope(g):
    if g.n > 5 or not (1 <= len(g.external) <= 3) or g.loops not in (1, 2):
        raise ValueError("Supported scope: 1-2 loops, 1-3 external legs, <=5 vertices")
    if degrees(g) != (3,) * g.n or g.vertex_kinds != ("phi3",) * g.n:
        raise ValueError("Input must contain bare trivalent scalar vertices")
    if not is_1pi(range(g.n), g.edges):
        raise ValueError("Input must be connected and 1PI")


def _encoding(g, perm, fixed):
    counts = Counter(tuple(sorted((perm[a], perm[b]))) for a, b in g.edges)
    kinds = [""] * g.n
    for v in range(g.n):
        kinds[perm[v]] = g.vertex_kinds[v]
    if fixed:
        ext = tuple((label, perm[v]) for label, v in g.external)
    else:
        ext = tuple(sorted(Counter(perm[v] for _, v in g.external).items()))
    return (tuple(kinds), ext, tuple(counts[a, b] for a in range(g.n) for b in range(a, g.n)))


@lru_cache(maxsize=None)
def canonical(g, fixed=True):
    best = None
    chosen = None
    for p in permutations(range(g.n)):
        key = _encoding(g, p, fixed)
        if best is None or key < best:
            best = key
            chosen = p
    kinds = [""] * g.n
    for v in range(g.n):
        kinds[chosen[v]] = g.vertex_kinds[v]
    result = Graph(
        g.n,
        tuple(sorted(tuple(sorted((chosen[a], chosen[b]))) for a, b in g.edges)),
        tuple(sorted((label, chosen[v]) for label, v in g.external)),
        tuple(kinds),
    )
    return best, result


@lru_cache(maxsize=None)
def symmetry(g, fixed=True):
    identity = _encoding(g, tuple(range(g.n)), fixed)
    vertex_aut = sum(_encoding(g, p, fixed) == identity for p in permutations(range(g.n)))
    multiplicities = Counter(g.edges)
    edge_aut = prod(
        factorial(count) * (2**count if a == b else 1) for (a, b), count in multiplicities.items()
    )
    # For this scope external legs attach at distinct vertices. In the general
    # unfixed case external half-edges at the same vertex can also permute.
    external_aut = (
        1 if fixed else prod(factorial(v) for v in Counter(v for _, v in g.external).values())
    )
    return vertex_aut * edge_aut * external_aut


def weight(g, fixed=True):
    return Fraction(1, symmetry(g, fixed))


@lru_cache(maxsize=None)
def generate(loops, legs):
    if loops not in (1, 2) or legs not in (1, 2, 3):
        raise ValueError("Generation supports loops=1,2 and legs=1,2,3")
    n = 2 * loops + legs - 2
    pairs = tuple((a, b) for a in range(n) for b in range(a, n))
    target = [2 if v < legs else 3 for v in range(n)]
    seen = {}

    def visit(i, remaining, chosen):
        if i == len(pairs):
            if any(remaining):
                return
            edges = tuple(pair for pair, m in zip(pairs, chosen) for _ in range(m))
            if not is_1pi(range(n), edges):
                return
            g = Graph(n, edges, tuple((label + 1, label) for label in range(legs)))
            key, cg = canonical(g)
            seen[key] = cg
            return
        a, b = pairs[i]
        maximum = remaining[a] // 2 if a == b else min(remaining[a], remaining[b])
        for count in range(maximum + 1):
            r = remaining.copy()
            r[a] -= count
            r[b] -= count
            if any(r[v] and not any(v in p for p in pairs[i + 1 :]) for v in range(n)):
                continue
            visit(i + 1, r, chosen + [count])

    visit(0, target, [])
    return tuple(seen[k] for k in sorted(seen))


def subgraph(g, edge_ids):
    """Retain every omitted incident half-edge as a DISTINCT external port."""
    ids = frozenset(edge_ids)
    vertices = sorted({v for i in ids for v in g.edges[i]})
    mapping = {v: i for i, v in enumerate(vertices)}
    boundary = []
    for label, v in g.external:
        if v in mapping:
            boundary.append(("external", label, 0, v))
    for i, (a, b) in enumerate(g.edges):
        if i in ids:
            continue
        if a in mapping:
            boundary.append(("cut", i, 0, a))
        if b in mapping:
            boundary.append(("cut", i, 1, b))
    ext = tuple((i + 1, mapping[x[-1]]) for i, x in enumerate(boundary))
    sg = Graph(
        len(vertices),
        tuple((mapping[g.edges[i][0]], mapping[g.edges[i][1]]) for i in sorted(ids)),
        ext,
        tuple(g.vertex_kinds[v] for v in vertices),
    )
    return sg, tuple(boundary)


def divergent_subgraphs(g, dimension=6):
    result = []
    for size in range(1, len(g.edges)):
        for ids in combinations(range(len(g.edges)), size):
            sg, boundary = subgraph(g, ids)
            if sg.loops < 1 or not is_1pi(range(sg.n), sg.edges):
                continue
            omega = dimension * sg.loops - 2 * len(sg.edges)
            if omega >= 0:
                result.append(
                    {
                        "edge_ids": list(ids),
                        "vertices": sorted({v for i in ids for v in g.edges[i]}),
                        "loops": sg.loops,
                        "external_legs": len(boundary),
                        "omega": omega,
                        "symmetry": symmetry(sg),
                        "boundary_ports": [list(x) for x in boundary],
                    }
                )
    return result


def henges(g, dimension=6):
    result = []
    for deleted in range(len(g.edges)):
        rest = tuple(i for i in range(len(g.edges)) if i != deleted)
        bridges = tuple(
            i for i in rest if not connected(range(g.n), tuple(g.edges[j] for j in rest if j != i))
        )
        retained = tuple(i for i in rest if i not in bridges)
        blocks = components(range(g.n), tuple(g.edges[i] for i in retained))
        block_records = []
        mapping = {}
        for j, vertices in enumerate(blocks):
            for v in vertices:
                mapping[v] = j
            ids = [i for i in retained if g.edges[i][0] in vertices]
            loops = len(ids) - len(vertices) + 1
            block_records.append(
                {
                    "vertices": sorted(vertices),
                    "edge_ids": ids,
                    "loops": loops,
                    "omega": dimension * loops - 2 * len(ids),
                    "needs_uv_subtraction": loops > 0 and dimension * loops - 2 * len(ids) >= 0,
                }
            )
        skeleton = tuple(
            (mapping[g.edges[i][0]], mapping[g.edges[i][1]]) for i in (*bridges, deleted)
        )
        result.append(
            {
                "minimum_edge": deleted,
                "blocks": block_records,
                "bridges": list(bridges),
                "skeleton_edges": [list(e) for e in skeleton],
                "skeleton_loops": len(skeleton) - len(blocks) + 1,
            }
        )
    return result


def contract(g, subgraphs):
    """Contract a vertex-disjoint spinney; KEEP omitted edges, even new self-loops."""
    used = set()
    mapping = {}
    kinds = []
    removed = set()
    for record in subgraphs:
        vertices = set(record["vertices"])
        if used & vertices:
            raise ValueError("Contract only a vertex-disjoint spinney")
        used |= vertices
        new = len(kinds)
        kinds.append("C" + str(record["external_legs"]))
        mapping.update({v: new for v in vertices})
        removed.update(record["edge_ids"])
    for v in range(g.n):
        if v not in mapping:
            mapping[v] = len(kinds)
            kinds.append(g.vertex_kinds[v])
    return Graph(
        len(kinds),
        tuple((mapping[a], mapping[b]) for i, (a, b) in enumerate(g.edges) if i not in removed),
        tuple((label, mapping[v]) for label, v in g.external),
        tuple(kinds),
    )
