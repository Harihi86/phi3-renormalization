"""Compose exact graph analysis into a versioned, portable JSON data set."""
from fractions import Fraction
from . import __version__
from .graphs import generate,divergent_subgraphs,henges,symmetry,validate_scope
from .renormalization import forests,operator_terms,structural_checks,insertion_columns,local_basis
from .routing import momentum_routing
from .validation import wick_check,zero_leg_insertion_check
from .examples import demo
from .parametric import symanzik
from .normalization import one_loop_pole_benchmarks
from .finite import finite_benchmark
from .overlap import overlap_benchmark
from .overlap_bounds import overlap_bound
from .two_loop import two_loop_benchmarks
from .two_loop_registry import SCHEME, two_loop_catalog


def analyze(g,graph_id,dimension=6):
    validate_scope(g)
    if dimension not in (4,6):raise ValueError("Use dimension 4 or 6")
    subs=divergent_subgraphs(g,dimension);hs=henges(g,dimension)
    omega=dimension*g.loops-2*len(g.edges)
    return {"id":graph_id,**g.to_dict(),"loops":g.loops,"external_legs":len(g.external),
            "omega":omega,"symmetry":symmetry(g),"weight":str(Fraction(1,symmetry(g))),
            "local_basis":local_basis(len(g.external),omega) if omega>=0 else [],
            "proper_divergent_subgraphs":subs,"henges":hs,
            "proper_forests":[list(f) for f in forests(subs)],"forest_terms":operator_terms(subs,omega),
            "routing":momentum_routing(g),"checks":structural_checks(g,subs,hs),
            "symanzik":symanzik(g)}


def build_dataset(loop_orders=(1,2),leg_counts=(1,2,3),dimension=6,wick=False):
    graphs={};records=[];validation=[];summaries=[]
    for loops in loop_orders:
        for legs in leg_counts:
            gs=generate(loops,legs)
            summaries.append({"loops":loops,"legs":legs,"count":len(gs),
                              "weight_sum":str(sum((Fraction(1,symmetry(g)) for g in gs),Fraction()))})
            for i,g in enumerate(gs,1):
                ident=f"L{loops}E{legs}G{i}";graphs[ident]=g
                records.append(analyze(g,ident,dimension))
            if wick:validation.append(wick_check(loops,legs))
    columns=insertion_columns(graphs,{r["id"]:r["proper_divergent_subgraphs"] for r in records})
    overlap=overlap_benchmark() if dimension==6 and 2 in loop_orders and 2 in leg_counts else None
    if overlap:overlap['analytic_bounds']=overlap_bound(1)
    requested_two_loop=[r['id'] for r in records if r['loops']==2]
    finite_two_loop=two_loop_benchmarks(requested_two_loop) if dimension==6 else []
    return {"schema_version":4,"software_version":__version__,"dimension":dimension,
            "method_status":{"primary":"BPH: Henge/Bogoliubov recursive R operation",
                             "forest_role":"Auxiliary explicit expansion and finite-graph combinatorial crosscheck",
                             "general_proof":"Not established by software tests; analytic scale bounds remain separate",
                             "terminology_note":"The supervisor's BPH naming remark concerns historical attribution; it does not state that Zimmermann's forest formula is unproved."},
            "scope":"Massive Euclidean single-real-scalar phi^3; connected amputated 1PI graphs; fixed external labels; 1-2 loops, E=1-3; tadpoles retained.",
            "convention":"Interaction +g phi^3/3!, bare Feynman vertex -g; positive symmetry denominator S, diagram weight 1/S. This is not the overall effective-action Gamma sign convention.",
            "limits":["Forest output is an ordered Taylor-operator expression, not a general integrated amplitude.",
                      "Three representative rational subtracted integrands are implemented separately with explicitly declared routings.",
                      "No general analytic convergence/uniform-momentum bound is proved by these finite checks.",
                      "The full zero-momentum Taylor scheme fixes the graphwise renormalized one-point amplitude to zero; no independent statement about the bare tadpole or physical vacuum is inferred.",
                      "Henge regions are distinct integration regions: repeated spinneys are not extra multiplicities.",
                      "Exact Euclidean Symanzik polynomials and pySecDec exports describe bare dimensionally regulated integrals, not BPH-renormalized amplitudes.",
                      "Finite integration covers the supported six-dimensional two-loop atlas in the declared massive Euclidean Taylor scheme. No beta function, arbitrary-loop integrator or general all-loop theorem is claimed.",
                      "Analytic bounds and original momentum-space convergence proofs are recorded separately by topology. Numerical order differences are not rigorous quadrature error bounds."],
            "summaries":summaries,"graphs":records,"counterterm_columns":columns,
            "wick_validation":validation,"zero_leg_insertion_validation":[zero_leg_insertion_check(l) for l in loop_orders],
            "representative_integrands":demo() if dimension==6 else {},
            "one_loop_pole_benchmarks":one_loop_pole_benchmarks(dimension),
            "finite_one_loop_benchmarks":finite_benchmark() if dimension==6 else [],
            "two_loop_overlap":overlap,
            "two_loop_finite":finite_two_loop,
            "two_loop_scheme":SCHEME if dimension==6 else None,
            "two_loop_catalog":[row for row in two_loop_catalog() if row['id'] in requested_two_loop]}
