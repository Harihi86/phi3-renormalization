"""Reproducible two-loop three-point counterterm closure build (no integrals)."""
import csv
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from fractions import Fraction
from .ct_bph import build_path_a
from .ct_action import build_path_b
from .ct_verify import compare_paths, check_henge
from .ct_wick import enumerate_wick
from .ct_graph import identity
from .graphs import Graph, generate, symmetry

VERSION='4.3.0-ct1'
CONVENTIONS={
 'dimension':6,'signature':'Euclidean','mass_squared':'m^2 > 0',
 'action':'S = integral [ (partial phi)^2/2 + m^2 phi^2/2 + g phi^3/3! + a phi + delta_m2 phi^2/2! + delta_Z (partial phi)^2/2 + delta_g phi^3/3! ]',
 'functional':'Z[J] = exp(-sum lambda_r/r! integral (delta/delta J)^r) exp(J D J/2); normalized by Z[0]',
 'connected':'W[J] = log(Z[J]/Z[0]); differentiate first, set J=0 afterwards',
 'amputation':'Strip one external free propagator per differentiated source; does not remove internal bridges',
 'one_pi':'Connected and no internal propagator bridge; external labels 1,2,3 fixed pointwise',
 'bare_vertex':'-g','counterterm_vertex':'negative of its coefficient in the action',
 'positive_weight':'1/S; NOT the signed amplitude coefficient',
 'bare_two_loop_three_point':'(-g)^5/S times I_G = -g^5 I_G/S',
 'loop_measure':'d^6q/pi^3 per loop, matching retained numerical core',
 'two_point_taylor':'total degree <=2 in the independent vector p; not degree 2 in p^2',
 'three_point_taylor':'degree 0 in two independent vectors p1,p2; p3=-p1-p2',
 'subgraph_momenta':'Taylor scale applies only inside the selected subgraph. Its loop momenta remain fixed; outer loops acting as boundary momenta are scaled only in the subgraph factors.',
 'tadpoles':'Retained in Z and in enumeration. No normal ordering or zero-VEV assumption. Linear finite counterterm remains unspecified; C1 hosts cannot enter this amputated 1PI three-point table because of an internal bridge.',
 'local_symbols':{
  'C1':'Formal one-point local vertex; finite linear condition unspecified. It is structurally absent from the target 1PI family.',
  'C2m':'vertex -g^2 K_m/2; K_m = integral D(q)^2',
  'C2z':'vertex -g^2 K_z p^2/2; K_z = integral [-D(q)^3 + (4/6) q^2 D(q)^4]',
  'C3':'vertex +g^3 K_3; K_3 = integral D(q)^3',
  'C3L2':'vertex -delta_g^(2); fixed as -T0 of path B bare + lower-CT action diagrams'},
 'tensor_reduction':'For two-point local coefficients, odd linear terms integrate to zero under q -> -q; q_mu q_nu integrates to delta_mu_nu q^2/6 in a rotationally invariant common UV regulator.',
 'symbolic_integrals':'K_m,K_z,K_3 and two-loop projected integrals remain regulated symbolic objects; this build computes exact combinatorics, not their numerical integrals',
 'external_reference':'pending_external_graph_by_graph_comparison; no missing reference image has been inferred',
 'precision_162':'Research ledger: 24 strict analytic, 48 empirical dual-kernel, 90 unresolved. This exact build performs no precision acceptance or large integration',
}


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def column_name(c):
    g=Graph.from_dict(c['host']);ct=next(v for v,k in enumerate(g.vertex_kinds) if k!='phi3')
    if c['operator_type']=='C3':
        return 'G'+''.join(str(label) for label,v in g.external if v==ct)
    neighbors={b if a==ct else a for a,b in g.edges if ct in (a,b)}
    legs=''.join(str(label) for label,v in g.external if v in neighbors)
    return ('M' if c['operator_type']=='C2m' else 'Z')+legs


def build(out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    def progress(stage):
        write_json(out/'progress.json',{'stage':stage,'time_utc':datetime.now(timezone.utc).isoformat(),'complete':False})
        print(stage,flush=True)
    progress('Path A: recursive BPH and all Henges')
    a=build_path_a()
    progress('Path B: counterterm-action Wick expansion and port gluing')
    b=build_path_b()
    progress('Exact cell comparison, including structural zeros')
    verification=compare_paths(a,b)
    if not verification['all_equal']:raise ArithmeticError('Exact correspondence mismatch')
    for r in a['rows']:
        if not check_henge(Graph.from_dict(r['graph']),r['henges']):raise ArithmeticError('Henge edge-map oracle failed')
    progress('Full generating functional and low-loop count audit')
    counts=[]
    for loops in (1,2):
        for legs in (1,2,3):
            audit=enumerate_wick(('phi3',)*(2*loops+legs-2),tuple(range(1,legs+1)),True)
            accepted=[r for r in audit['records'] if r['accepted']]
            expected={identity(g):str(Fraction(1,symmetry(g))) for g in generate(loops,legs)}
            actual={r['id']:r['positive_weight'] for r in accepted}
            if actual!=expected:raise ArithmeticError('Independent Wick catalogue mismatch')
            write_json(out/'audit'/f'functional-L{loops}E{legs}.json',audit)
            stages={'functional_classes':len(audit['records']),
                    'vacuum_factor_classes':sum(bool(r['vacuum_components']) for r in audit['records']),
                    'normalized_classes':sum(not r['vacuum_components'] for r in audit['records']),
                    'connected_amputated_classes':sum(r['connected'] for r in audit['records']),
                    'one_pi_classes':len(accepted),'unfixed_orbits':len({identity(Graph.from_dict(r['graph']),False) for r in accepted})}
            counts.append({'loops':loops,'legs':legs,**stages,'weight_sum':str(sum((Fraction(r['positive_weight']) for r in accepted),Fraction())),
                           'pairings_visited_in_orbit_representatives':audit['pairings_visited'],
                           'full_labelled_wick_multiplicity':audit['full_wick_multiplicity'],
                           'audit_file':f'audit/functional-L{loops}E{legs}.json','exact_wick_equals_adjacency_and_automorphisms':True})
    b_export={k:v for k,v in b.items() if k!='audits'}
    b_export['audit_files']={}
    for i,(key,audit) in enumerate(b['audits'].items()):
        name=f'audit/action-{i:02d}.json';write_json(out/name,audit);b_export['audit_files'][key]=name
    cols=sorted(a['columns'],key=lambda c:({'C2m':0,'C2z':1,'C3':2}[c['operator_type']],column_name(c)))
    for c in cols:c['name']=column_name(c)
    data={'schema':'phi3lab.counterterm-closure.v1','software_version':VERSION,'conventions':CONVENTIONS,
          'rows':a['rows'],'columns':cols,'verification':verification,'count_audit':counts,
          'action_sectors':b['sectors'],'path_files':{'A':'path-a.json','B':'path-b.json'},
          'environment':{'python':sys.version,'platform':platform.platform()},
          'references':[{'title':'A Simple Proof of the BPH Theorem (Kennedy, 1996)','url':'https://arxiv.org/abs/hep-th/9612113'},
                        {'title':'Automating Renormalization of Quantum Field Theories (Kennedy, Binoth, Rippon, 2007)','url':'https://arxiv.org/abs/0712.1016'}]}
    write_json(out/'path-a.json',a);write_json(out/'path-b.json',b_export)
    write_json(out/'correspondence.json',data);write_json(out/'verification.json',verification)
    write_json(out/'conventions.json',CONVENTIONS)
    names={c['id']:c['name'] for c in cols};names.update(bare='I',C3L2='C3L2')
    old={r['id']:r['old_id'] for r in a['rows']}
    with (out/'table.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.writer(f);writer.writerow(['graph','stable_id','column','exact_A','exact_B','equal','zero_reason'])
        for cell in verification['cells']:
            writer.writerow([old[cell['row']],cell['row'],names[cell['column']],json.dumps(cell['A']),json.dumps(cell['B']),cell['equal'],cell['zero_reason'] or ''])
    progress('Bilingual workbench, SVG and editable TikZ/LaTeX')
    from .ct_report import write_reports
    write_reports(data,a,b_export,out)
    write_json(out/'completion.json',{'exit_code':0,'complete':True,'version':VERSION,
              'rows':verification['rows'],'cells':verification['checked_cells'],'exact_equal':True,
              'time_utc':datetime.now(timezone.utc).isoformat()})
    write_json(out/'progress.json',{'stage':'complete','complete':True})
    print(json.dumps({'rows':verification['rows'],'columns':verification['columns'],'cells':verification['checked_cells'],
                      'nonzero_cells':verification['nonzero_cells'],'symbolic_terms':verification['symbolic_monomials'],'out':str(out)}),flush=True)
    return data
