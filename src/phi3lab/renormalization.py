"""BPH recursive subtraction expansion with auxiliary forest checks.

These are combinatorial statements. They are NOT an analytic convergence proof.
Taylor operators retain a typed polynomial space, not an untyped contracted point.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from .graphs import canonical, contract, symmetry


def compatible(a,b):
    ea,eb=set(a["edge_ids"]),set(b["edge_ids"])
    return ea<=eb or eb<=ea or set(a["vertices"]).isdisjoint(b["vertices"])


def forests(subs):
    result=[]
    for size in range(len(subs)+1):
        for chosen in combinations(range(len(subs)),size):
            if all(compatible(subs[a],subs[b]) for a,b in combinations(chosen,2)):
                result.append(chosen)
    return tuple(result)


def spinneys(subs):
    return tuple(chosen for size in range(len(subs)+1)
                 for chosen in combinations(range(len(subs)),size)
                 if all(set(subs[a]["vertices"]).isdisjoint(subs[b]["vertices"])
                        for a,b in combinations(chosen,2)))


def bogoliubov_forests(subs):
    """Independently expand recursive disjoint spinneys into embedded forests."""
    def recurse(available):
        out={frozenset()}
        for size in range(1,len(available)+1):
            for top in combinations(available,size):
                if not all(set(subs[a]["vertices"]).isdisjoint(subs[b]["vertices"])
                           for a,b in combinations(top,2)):continue
                choices=[]
                for a in top:
                    proper=[b for b in available if set(subs[b]["edge_ids"])<set(subs[a]["edge_ids"])]
                    choices.append(recurse(proper))
                for children in product(*choices):
                    out.add(frozenset(top).union(*children))
        return out
    return recurse(tuple(range(len(subs))))


def local_basis(legs,omega):
    if legs==0:return ["vacuum constant"]
    if legs==1:return ["linear phi; no independent external momentum"]
    if legs==2:return ["mass: phi^2"]+(["kinetic: (partial phi)^2"] if omega>=2 else [])
    if legs==3 and omega==0:return ["coupling: phi^3"]
    return [f"local external-momentum polynomial of degree <= {omega}"]


def operator_terms(subs,omega):
    """Expand recursive Bogoliubov subtractions; forest enumeration is an oracle.

    This is a bounded combinatorial expression, not a general proof of the
    forest theorem or a numerical implementation of all Taylor operators.
    """
    result=[]
    chosen_sets=sorted(bogoliubov_forests(subs),key=lambda f:(len(f),tuple(sorted(f))))
    for chosen in chosen_sets:
        ordered=sorted(chosen,key=lambda j:(len(subs[j]["edge_ids"]),j))
        operators=[{"subgraph":f"gamma{j+1}","edge_ids":subs[j]["edge_ids"],
                    "degree":subs[j]["omega"],"basis":local_basis(subs[j]["external_legs"],subs[j]["omega"])}
                   for j in ordered]
        result.append({"sign":(-1)**len(chosen),"operators_inner_first":operators})
        if omega>=0:
            result.append({"sign":(-1)**(len(chosen)+1),
                           "operators_inner_first":operators+[{"subgraph":"G","degree":omega}]})
    return result


def structural_checks(g,subs,hs):
    direct={frozenset(f) for f in forests(subs)}
    recursive=bogoliubov_forests(subs)
    coverage=[]
    for f in spinneys(subs):
        internal=set().union(*(set(subs[j]["edge_ids"]) for j in f))
        expected=set(range(len(g.edges)))-internal
        actual=set()
        for h in hs:
            if all(any(set(subs[j]["edge_ids"])<=set(b["edge_ids"]) for b in h["blocks"]) for j in f):
                actual.add(h["minimum_edge"])
        coverage.append({"spinney":[j+1 for j in f],"allowed_minimum_edges":sorted(actual),
                         "cograph_edges":sorted(expected),"ok":expected==actual})
    return {"forest_equals_recursive_bogoliubov":direct==recursive,
            "one_loop_henge_skeletons":all(h["skeleton_loops"]==1 for h in hs),
            "spinney_coverage":coverage,
            "all_spinney_coverage_ok":all(c["ok"] for c in coverage),
            "qualification":"Combinatorial support/coverage only; no numerical sector identity or analytic convergence asserted."}


def insertion_columns(graphs,all_subgraphs):
    """One-loop counterterms inserted into one-loop cographs at two-loop order.

    A local countervertex includes its own 1/S_gamma. Thus a parent embedding
    contributes S_gamma/S_parent to the corresponding countergraph column.
    Multiple embeddings are summed EXACTLY, without multiplying by Henge counts.
    """
    columns={}
    for graph_id,g in graphs.items():
        for sub in all_subgraphs[graph_id]:
            if sub["loops"]!=1 or g.loops!=2:continue
            cg=contract(g,[sub]);key,cg=canonical(cg)
            if key not in columns:
                columns[key]={"cograph":cg.to_dict(),"counterterm_type":f'C{sub["external_legs"]}',
                              "basis":local_basis(sub["external_legs"],sub["omega"]),
                              "expected_coefficient":str(Fraction(1,symmetry(cg))),"entries":[],"total":Fraction()}
            ratio=Fraction(sub["symmetry"],symmetry(g))
            columns[key]["entries"].append({"parent":graph_id,"edge_ids":sub["edge_ids"],"coefficient":str(ratio)})
            columns[key]["total"]+=ratio
    result=[]
    for k in sorted(columns):
        c=columns[k];c["total"]=str(c["total"])
        c["coefficient_matches_cograph_symmetry"]=Fraction(c["total"])==Fraction(c["expected_coefficient"])
        result.append(c)
    return result
