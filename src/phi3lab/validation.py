"""Independent labelled-half-edge Wick enumeration and exact weight comparison."""
from collections import Counter
from fractions import Fraction
from math import factorial
from .graphs import Graph, canonical, generate, is_1pi, weight


def pairings(items):
    if not items:
        yield ();return
    a=items[0]
    for i in range(1,len(items)):
        b=items[i]
        for rest in pairings(items[1:i]+items[i+1:]):yield ((a,b),)+rest


def wick_check(loops,legs):
    graphs=generate(loops,legs);n=2*loops+legs-2
    stubs=tuple(v for v in range(n) for _ in range(2 if v<legs else 3))
    # Unlike generation by adjacency matrices, this visits every labelled pairing.
    counts=Counter();visited=0
    for matching in pairings(tuple(range(len(stubs)))):
        visited+=1
        edges=tuple(tuple(sorted((stubs[a],stubs[b]))) for a,b in matching)
        if not is_1pi(range(n),edges):continue
        key,_=canonical(Graph(n,tuple(sorted(edges)),tuple((i+1,i) for i in range(legs))))
        counts[key]+=1
    normalization=Fraction(factorial(n)//factorial(n-legs)*3**legs,
                           factorial(n)*factorial(3)**n)
    expected={canonical(g)[0]:weight(g) for g in graphs}
    actual={k:count*normalization for k,count in counts.items()}
    return {"loops":loops,"legs":legs,"pairings_visited":visited,
            "fixed_attachment_pairings":[counts[canonical(g)[0]] for g in graphs],
            "normalization":str(normalization),"weights":[str(actual.get(canonical(g)[0],0)) for g in graphs],
            "ok":actual==expected,"method":"Labelled half-edge pairings; factor (V)_E 3^E / [V! (3!)^V]."}


def zero_leg_insertion_check(loops):
    """Bare combinatorial check only: split every line of a two-point graph.

    Label the new external leg 3. An inserted bare -g and the mass derivative
    sign agree for Feynman amplitudes; here we compare positive 1/S weights.
    Different Taylor orders invalidate a naive renormalized identity; no such
    identity is asserted by this topology/weight test.
    """
    counts=Counter();embeddings=0
    for g in generate(loops,2):
        for i,(a,b) in enumerate(g.edges):
            inserted=Graph(g.n+1,g.edges[:i]+g.edges[i+1:]+((a,g.n),(g.n,b)),g.external+((3,g.n),))
            key,_=canonical(inserted);counts[key]+=weight(g);embeddings+=1
    expected={canonical(g)[0]:weight(g) for g in generate(loops,3)}
    return {"loops":loops,"insertion_positions":embeddings,"three_point_graphs":len(counts),
            "weights":[str(counts[canonical(g)[0]]) for g in generate(loops,3)],
            "ok":dict(counts)==expected,"qualification":"Bare topology and weight identity only; no renormalized mass-derivative identity claimed."}
