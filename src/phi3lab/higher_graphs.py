"""Exact colored multigraph generation through four phi-cubed loops.

This is a graph enumerator, not a catalogue of topologies or an integral
compiler.  It returns the same Graph objects as graphs.py without extending
that module's numerical coverage.  No optional dependency is required.

Completeness: in a connected 1PI trivalent graph with at most three external
legs the external attachments are distinct, except for the one-vertex
one-leg tadpole (which is included).  Fill the upper adjacency matrix one
row at a time with the prescribed remaining degrees.  Unprocessed vertices
with identical colors and processed columns are interchangeable; sorting
the entries just assigned inside each such class loses no isomorphism
class.  Final individualization/refinement canonicalization removes the
remaining duplicates exactly.  Multiple edges and loops are discrete data.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial, prod
from time import monotonic

from .graphs import Graph, degrees, is_1pi


def _matrix(g):
    matrix=[[0]*g.n for _ in range(g.n)]
    for a,b in g.edges:
        matrix[a][b]+=1
        if a!=b:matrix[b][a]+=1
    return tuple(tuple(row) for row in matrix)


def _colors(g,fixed):
    attachments=[[] for _ in range(g.n)]
    for label,v in g.external:attachments[v].append(label)
    if fixed:
        return tuple((g.vertex_kinds[v],tuple(attachments[v])) for v in range(g.n))
    return tuple((g.vertex_kinds[v],len(attachments[v])) for v in range(g.n))


def _refine(matrix,partition):
    """Deterministic equitable refinement; loops have two incident ends."""
    while True:
        refined=[]
        for cell in partition:
            groups={}
            for v in cell:
                signature=(matrix[v][v],tuple(sum(matrix[v][w]*(2 if v==w else 1)
                                                 for w in block) for block in partition))
                groups.setdefault(signature,[]).append(v)
            refined.extend(tuple(groups[key]) for key in sorted(groups))
        refined=tuple(refined)
        if len(refined)==len(partition):return refined
        partition=refined


def _label_data(g,fixed):
    matrix=_matrix(g);colors=_colors(g,fixed);groups={}
    for v,color in enumerate(colors):groups.setdefault(color,[]).append(v)
    initial=tuple(tuple(groups[color]) for color in sorted(groups))
    best=None;chosen=None;automorphisms=0;leaves=0

    def search(partition):
        nonlocal best,chosen,automorphisms,leaves
        partition=_refine(matrix,partition)
        candidates=[(len(cell),i) for i,cell in enumerate(partition) if len(cell)>1]
        if candidates:
            _,index=min(candidates);cell=partition[index]
            for v in cell:
                search(partition[:index]+((v,),tuple(w for w in cell if w!=v))+partition[index+1:])
            return
        order=tuple(cell[0] for cell in partition)
        key=(tuple(colors[v] for v in order),
             tuple(matrix[order[a]][order[b]] for a in range(g.n) for b in range(a,g.n)))
        leaves+=1
        if best is None or key<best:best=key;chosen=order;automorphisms=1
        elif key==best:automorphisms+=1
    search(initial)
    mapping={v:i for i,v in enumerate(chosen)}
    canonical_graph=Graph(g.n,tuple(sorted(tuple(sorted((mapping[a],mapping[b]))) for a,b in g.edges)),
                          tuple((label,mapping[v]) for label,v in g.external),
                          tuple(g.vertex_kinds[v] for v in chosen))
    return best,canonical_graph,automorphisms,leaves


@lru_cache(maxsize=32768)
def _canonical_data(g,fixed=True):
    if not isinstance(g,Graph):raise ValueError('Use a Graph with exact integer vertices and edges')
    if type(fixed)is not bool:raise ValueError('fixed must be a boolean')
    return _label_data(g,fixed)


def canonical(g,fixed=True):
    """Exact colored canonical key and representative; no factorial scan."""
    if not isinstance(g,Graph):raise ValueError('Use a Graph with exact integer vertices and edges')
    if type(fixed)is not bool:raise ValueError('fixed must be a boolean')
    key,representative,_,_=_canonical_data(g,fixed)
    return key,representative


def symmetry(g,fixed=True):
    """Full half-edge symmetry: vertex automorphisms, parallels and loop flips."""
    if not isinstance(g,Graph):raise ValueError('Use a Graph with exact integer vertices and edges')
    if type(fixed)is not bool:raise ValueError('fixed must be a boolean')
    vertex_aut=_canonical_data(g,fixed)[2]
    edge_aut=prod(factorial(count)*(2**count if a==b else 1)
                  for (a,b),count in Counter(g.edges).items())
    external_aut=1 if fixed else prod(factorial(count) for count in Counter(v for _,v in g.external).values())
    return vertex_aut*edge_aut*external_aut


def weight(g,fixed=True):
    return Fraction(1,symmetry(g,fixed))


def validate_scope(g):
    if not isinstance(g,Graph):raise ValueError('Use a Graph input')
    if not 1<=len(g.external)<=3 or not 1<=g.loops<=4 or g.n>9:
        raise ValueError('Higher graph scope: 1–4 loops, 1–3 external legs, at most9 vertices')
    if degrees(g)!=(3,)*g.n or g.vertex_kinds!=('phi3',)*g.n:
        raise ValueError('Use bare trivalent scalar vertices')
    if not is_1pi(range(g.n),g.edges):raise ValueError('Use a connected 1PI graph')


def _arguments(loops,legs):
    if type(loops)is not int or type(legs)is not int or not 1<=loops<=4 or not 1<=legs<=3:
        raise ValueError('Generation requires integer loops=1–4 and legs=1–3')


@lru_cache(maxsize=12)
def _generate_data(loops,legs):
    _arguments(loops,legs)
    start=monotonic();n=2*loops+legs-2
    target=[2 if v<legs else 3 for v in range(n)]
    matrix=[[0]*n for _ in range(n)]
    remaining=target.copy();seen={}
    stats={'rows_visited':0,'degree_solutions':0,'one_pi_labeled_candidates':0,
           'canonical_search_leaves':0}
    colors=tuple((v+1,) if v<legs else () for v in range(n))

    def row(i):
        stats['rows_visited']+=1
        if i==n:
            if any(remaining):return
            stats['degree_solutions']+=1
            edges=tuple((a,b) for a in range(n) for b in range(a,n) for _ in range(matrix[a][b]))
            if not is_1pi(range(n),edges):return
            stats['one_pi_labeled_candidates']+=1
            g=Graph(n,edges,tuple((label+1,label) for label in range(legs)))
            key,cg,_,leaves=_canonical_data(g)
            stats['canonical_search_leaves']+=leaves;seen[key]=cg
            return
        # A loop in a multi-vertex trivalent graph leaves at most one internal
        # connection. It is then disconnected or has a bridge, so cannot be1PI.
        loop_counts=range(remaining[i]//2+1) if n==1 else (0,)
        for loops_here in loop_counts:
            matrix[i][i]=loops_here
            need=remaining[i]-2*loops_here
            classes={}
            for v in range(i+1,n):
                history=(colors[v],tuple(matrix[a][v] for a in range(i)))
                classes.setdefault(history,[]).append(v)
            previous={}
            for cell in classes.values():
                for a,b in zip(cell,cell[1:]):previous[b]=a

            def fill(j,left):
                if j==n:
                    if left:return
                    remaining[i]=0
                    # No proper processed component may be already closed.
                    if i+1<n:
                        open_vertices={v for v in range(n) if remaining[v]}
                        reachable=set(open_vertices);stack=list(open_vertices)
                        while stack:
                            v=stack.pop()
                            for w in range(n):
                                if matrix[v][w] and w not in reachable:
                                    reachable.add(w);stack.append(w)
                        if open_vertices and len(reachable)<n:
                            remaining[i]=need+2*loops_here;return
                    row(i+1)
                    remaining[i]=need+2*loops_here
                    return
                maximum=min(left,remaining[j])
                if j in previous:maximum=min(maximum,matrix[i][previous[j]])
                future=sum(remaining[k] for k in range(j+1,n))
                for count in range(max(0,left-future),maximum+1):
                    matrix[i][j]=matrix[j][i]=count;remaining[j]-=count
                    fill(j+1,left-count)
                    remaining[j]+=count
                matrix[i][j]=matrix[j][i]=0
            fill(i+1,need)
            matrix[i][i]=0
    row(0)
    graphs=tuple(seen[key] for key in sorted(seen))
    stats.update({'loops':loops,'external_legs':legs,'vertices':n,'internal_edges':3*loops+legs-3,
                  'fixed_label_graph_count':len(graphs),'generation_seconds':monotonic()-start,
                  'algorithm':'Degree-constrained adjacency rows, future-vertex orbit ordering, exact colored individualization/refinement'})
    return graphs,stats


def generate(loops,legs):
    """All connected 1PI phi3 multigraphs with fixed external labels1..legs."""
    _arguments(loops,legs)
    return _generate_data(loops,legs)[0]


def generation_report(loops,legs):
    """Actual enumerated counts, exact weights and diagnostic work counts."""
    _arguments(loops,legs)
    graphs,stats=_generate_data(loops,legs)
    unlabeled={canonical(g,False)[0] for g in graphs}
    total=sum((weight(g) for g in graphs),Fraction())
    return {**stats,'unfixed_external_graph_count':len(unlabeled),
            'fixed_label_weight_sum':str(total),'fixed_label_weight_numerator':total.numerator,
            'fixed_label_weight_denominator':total.denominator,
            'external_convention':'External labels are fixed pointwise; internal vertices and parallel half-edges may permute',
            'integral_coverage':'Graph enumeration only; no numerical integral claimed by this module'}
