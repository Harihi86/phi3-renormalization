"""Cycle response polynomials with independently decorated branches."""
from itertools import combinations
from fractions import Fraction as Q
from radial_tree_polynomials import add,multiply,product,edge_sum

def compile_tree(tree,variables,edges):
    children=tree['children'];compiled=[compile_tree(c,variables,edges) for c in children]
    representatives=[min(v for e in c['edge_ids'] for v in edges[e]) for c in children]
    assert len(representatives)==len(set(representatives))
    branches=[];used=[]
    for arc in tree['arcs']:
        indices=[representatives.index(v) for v in arc['child_representatives']];used.extend(indices)
        factor=product(compiled[i][0] for i in indices)
        response=multiply(edge_sum(arc['edges']),factor)
        for i in indices:
            index=children[i]['source_part_index'];assert index is not None
            slot=1<<variables.index(index) if index in variables else 0
            response=add(response,product([compiled[i][1],*(compiled[j][0] for j in indices if j!=i),{(0,slot):Q(1)}]))
        branches.append((factor,response))
    assert sorted(used)==list(range(len(children)))
    assert len(branches)==tree['external_legs'] in (2,3)
    U=add(*(product([response,*(p for j,(p,q) in enumerate(branches) if i!=j)]) for i,(factor,response) in enumerate(branches)))
    # The target data have every external squared channel equal to one.
    F=add(*(product([branches[i][1],branches[j][1],*(p for h,(p,q) in enumerate(branches) if h not in (i,j))]) for i,j in combinations(range(len(branches)),2)))
    return U,F
