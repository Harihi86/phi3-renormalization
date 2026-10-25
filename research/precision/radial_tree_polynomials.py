"""Recursive graph-cycle U/F algebra, independent of Gaussian routing."""
from collections import defaultdict
from fractions import Fraction as Q

def add(*polys):
    out=defaultdict(Q)
    for poly in polys:
        for key,value in poly.items():out[key]+=value
    return {key:value for key,value in out.items() if value}

def multiply(a,b):
    out=defaultdict(Q)
    for (am,aslot),av in a.items():
        for (bm,bslot),bv in b.items():
            assert not am&bm and not aslot&bslot, 'Disjoint subtrees and independent Taylor variables'
            out[(am|bm,aslot|bslot)]+=av*bv
    return dict(out)

def product(polys):
    value={(0,0):Q(1)}
    for poly in polys:value=multiply(value,poly)
    return value

def edge_sum(edges):return {(1<<e,0):Q(1) for e in edges}

def compile_tree(tree,variables):
    children=tree['children'];compiled=[compile_tree(c,variables) for c in children]
    arcs=tree['arcs'];decorated=max(arcs,key=lambda a:len(a['edges']))
    ordinary=[a for a in arcs if a is not decorated]
    assert all(len(a['edges'])==1 for a in ordinary)
    prod=product(pair[0] for pair in compiled)
    response=multiply(edge_sum(decorated['edges']),prod)
    for j,child in enumerate(children):
        index=child['source_part_index'];assert index is not None
        slot=1<<variables.index(index) if index in variables else 0
        contribution=product([compiled[j][1],*(pair[0] for i,pair in enumerate(compiled) if i!=j),{(0,slot):Q(1)}])
        response=add(response,contribution)
    simple_edges=[a['edges'][0] for a in ordinary]
    U=add(response,multiply(edge_sum(simple_edges),prod))
    if tree['external_legs']==2:
        assert len(simple_edges)==1
        F=multiply(edge_sum(simple_edges),response)
    else:
        assert tree['external_legs']==3 and len(simple_edges)==2
        F=add(multiply(edge_sum(simple_edges),response),multiply({(sum(1<<e for e in simple_edges),0):Q(1)},prod))
    return U,F

def original_positive(state,family,invariants):
    out=defaultdict(Q)
    for row in state[family]:
        mask,slot,channel=row[:3];value=Q(row[3] if len(row)==4 else 1)
        if family=='F':value*=Q(invariants[channel])
        out[(mask,slot)]+=value
    return {key:value for key,value in out.items() if value}

def gaussian_full(descriptor,term,gram,target_variables):
    # Independent Cauchy--Binet construction using every original edge.
    from itertools import combinations
    import triple_symbolic_binding as G
    L=descriptor['loops'];P=descriptor['external_legs']-1;k=len(term['jet_subgraphs']);route=term['scaled_routing']
    assert set(term['jet_subgraphs'])==set(target_variables)
    results=[]
    for family,degree in ((0,L),(1,L+1)):
        out=defaultdict(Q)
        for edges in combinations(range(descriptor['edges']),degree):
            if family==0:
                det=G.minor(route,edges,tuple(range(L)),k);poly=G.multiply(det,det,k)
            else:
                dets=[G.minor(route,edges,(*range(L),L+i),k) for i in range(P)];poly=defaultdict(Q)
                for i in range(P):
                    for j in range(P):G.add(poly,G.multiply(dets[i],dets[j],k),Q(gram[i][j]))
            for powers,value in poly.items():
                if not value:continue
                assert all(p%2==0 for p in powers), 'Odd Taylor coefficient remains'
                slot=sum(1<<target_variables.index(v) for v,p in zip(term['jet_subgraphs'],powers) if p==2)
                out[(sum(1<<e for e in edges),slot)]+=value
        results.append({key:value for key,value in out.items() if value})
    return tuple(results)
