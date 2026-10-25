"""Run: python3 -m phi3lab build --wick; python3 -m phi3lab validate."""
import argparse,json,sys
from pathlib import Path
from .analysis import build_dataset,analyze
from .graphs import Graph
from .reporting import write_html,write_latex,svg
from .validation import wick_check
from .examples import demo
from .interop import write_pysecdec_export
from .provenance import write_manifest,verify_manifest
from .workbench import evaluate_overlap
from .two_loop import evaluate_two_loop


def run():
    parser=argparse.ArgumentParser(description="Auditable massive phi^3 recursive BPH graph compiler through four loops")
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('counterterms',help='Exact independent two-path counterterm correspondence; no integration')
    p.add_argument('--out',type=Path,default=Path('results/counterterms'))
    p=commands.add_parser('build',help='Generate the complete supported graph atlas and offline report')
    p.add_argument('--loops',nargs='+',type=int,choices=[1,2],default=[1,2])
    p.add_argument('--legs',nargs='+',type=int,choices=[1,2,3],default=[1,2,3])
    p.add_argument('--dimension',type=int,choices=[4,6],default=6)
    p.add_argument('--wick',action='store_true',help='Independently verify labelled half-edge Wick counts')
    p.add_argument('--out',type=Path,default=Path('results'))
    p=commands.add_parser('validate',help='Run all six independent Wick counting comparisons')
    p.add_argument('--out',type=Path)
    p=commands.add_parser('verify',help='Check current source and every recorded build artifact against the manifest')
    p.add_argument('--out',type=Path,default=Path('results'))
    p=commands.add_parser('analyze',help='Analyze one bare trivalent 1PI graph JSON')
    p.add_argument('input',type=Path);p.add_argument('--dimension',type=int,choices=[4,6],default=6)
    p.add_argument('--out',type=Path)
    commands.add_parser('examples',help='Evaluate exact rational sample points for three subtracted integrands')
    p=commands.add_parser('overlap',help='Evaluate the full finite two-loop BPH overlapping graph and its analytic bounds')
    p.add_argument('--s',type=float,default=1,help='Euclidean external momentum squared')
    p.add_argument('--m2',type=float,default=1,help='Positive common mass squared')
    p.add_argument('--order',type=int,default=16)
    p.add_argument('--out',type=Path)
    p=commands.add_parser('two-loop',help='Evaluate any supported two-loop graph in the same complete BPH Taylor scheme')
    p.add_argument('graph',help='Graph ID, e.g. L2E3G4')
    p.add_argument('--invariants',nargs='*',type=float,default=None,help='p^2 for two-point or p1^2 p2^2 p3^2 for three-point graphs')
    p.add_argument('--m2',type=float,default=1)
    p.add_argument('--order',type=int,default=None)
    p.add_argument('--out',type=Path)
    p=commands.add_parser('serve',help='Open the local interactive BPH workbench')
    p.add_argument('--port',type=int,default=8767)
    p.add_argument('--no-open',action='store_true')
    p=commands.add_parser('build-automatic',help='Compile all698 generated1--4-loop graphs to executable recursive finite densities')
    p.add_argument('--out',type=Path,default=Path('results/automatic'))
    p.add_argument('--language',choices=['cn','en'])
    p=commands.add_parser('automatic',help='Construct a finite BPH Gaussian density directly from a generated or input graph')
    p.add_argument('--loops',type=int,choices=[1,2,3,4],default=2)
    p.add_argument('--legs',type=int,choices=[1,2,3],default=3)
    p.add_argument('--index',type=int,default=1)
    p.add_argument('--input',type=Path)
    p.add_argument('--invariants',nargs='*',type=float,default=None)
    p.add_argument('--m2',type=float,default=1)
    p.add_argument('--integrate',action='store_true')
    p.add_argument('--power',type=int,default=14)
    p.add_argument('--replicates',type=int,default=8)
    p.add_argument('--beta',type=float,default=.35)
    p.add_argument('--sampler',choices=['auto','hepp','dirichlet','recursive','mixture'],default='auto')
    p.add_argument('--out',type=Path)
    args=parser.parse_args()
    try:
        if args.command=='counterterms':
            from .ct_closure import build
            build(args.out)
        elif args.command=='build-automatic':
            from .automatic_catalog import build_automatic
            data=build_automatic(args.out,language=args.language)
            print(json.dumps({'compiled_graphs':data['compiled_graphs'],'groups':data['groups'],'source_unchanged':data['source_unchanged'],'report':str((args.out/'index.html').resolve())},indent=2))
        elif args.command=='automatic':
            from .automatic_catalog import evaluate_automatic
            graph=Graph.from_dict(json.loads(args.input.read_text())) if args.input else None
            result=evaluate_automatic(args.loops,args.legs,args.index,graph=graph,invariants=args.invariants,m2=args.m2,integrate=args.integrate,power=args.power,replicates=args.replicates,beta=args.beta,sampler=args.sampler)
            output=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)
            if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(output)
            else:print(output)
        elif args.command=='build':
            data=build_dataset(tuple(sorted(set(args.loops))),tuple(sorted(set(args.legs))),args.dimension,args.wick)
            checks=all(r['checks']['forest_equals_recursive_bogoliubov'] and r['checks']['one_loop_henge_skeletons']
                       and r['checks']['all_spinney_coverage_ok'] and r['routing']['conservation_ok'] for r in data['graphs'])
            checks=checks and all(c['coefficient_matches_cograph_symmetry'] for c in data['counterterm_columns'])
            checks=checks and all(v['ok'] for v in data['wick_validation'])
            checks=checks and all(v['ok'] for v in data['zero_leg_insertion_validation'])
            checks=checks and all(all(r['symanzik']['checks'].values()) for r in data['graphs'])
            checks=checks and all(all(r['checks'].values()) for r in data['one_loop_pole_benchmarks'])
            checks=checks and all(all(r['checks'].values()) for r in data['finite_one_loop_benchmarks'])
            if data['two_loop_overlap']:
                checks=checks and all(data['two_loop_overlap']['checks'].values())
            if data['two_loop_finite']:
                checks=checks and all(all(r['checks'].values()) for r in data['two_loop_finite'])
                checks=checks and {r['graph'] for r in data['two_loop_finite']}=={r['id'] for r in data['graphs'] if r['loops']==2}
            if not checks:raise RuntimeError('Validation failed; refusing to publish a passing report')
            args.out.mkdir(parents=True,exist_ok=True)
            if args.dimension==6 and set(args.loops)=={1,2} and set(args.legs)=={1,2,3}:
                from .automatic_catalog import build_automatic
                automatic=build_automatic(args.out/'automatic')
                data['automatic_bph']={'compiled_graphs':automatic['compiled_graphs'],'groups':automatic['groups'],'source_fingerprint':automatic['source_fingerprint']}
            (args.out/'analysis.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
            write_latex(data,args.out/'forests.tex')
            artifacts=[args.out/'analysis.json',args.out/'forests.tex']
            artifacts.extend(write_pysecdec_export(data,args.out/'pysecdec'))
            (args.out/'figures').mkdir(exist_ok=True)
            (args.out/'inputs').mkdir(exist_ok=True)
            for record in data['graphs']:
                input_path=args.out/'inputs'/(record['id']+'.json')
                input_path.write_text(json.dumps(Graph.from_dict(record).to_dict(),indent=2),encoding='utf-8')
                artifacts.append(input_path)
                path=args.out/'figures'/(record['id']+'.svg')
                path.write_text(svg(Graph.from_dict(record)),encoding='utf-8');artifacts.append(path)
            write_html(data,args.out/'index.html');artifacts.append(args.out/'index.html')
            if (args.out/'automatic').exists():artifacts.extend(p for p in (args.out/'automatic').rglob('*') if p.is_file())
            for name in ('external_comparison.json','pysecdec_comparison.json','overlap_comparison.json','nested_triangle_comparison.json','native_primitive_comparison.json','test_verification.json','test-results.txt','extreme_triangle_comparison.json','extreme_native_comparison.json','extreme_bubble_comparison.json','extreme_overlap_comparison.json','higher_graph_feyngen_comparison.json','automatic_bph_comparison.json','automatic_integral_comparison.json','automatic_higher_integrals.json'):
                if (args.out/name).exists():artifacts.append(args.out/name)
            write_manifest(args.out,artifacts,data)
            print(f"Generated {len(data['graphs'])} graphs; structural checks passed. Report: {(args.out/'index.html').resolve()}")
        elif args.command=='validate':
            data={'checks':[wick_check(l,e) for l in (1,2) for e in (1,2,3)]}
            data['ok']=all(v['ok'] for v in data['checks'])
            s=json.dumps(data,indent=2)
            if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(s,encoding='utf-8')
            print(s)
            if not data['ok']:return 1
        elif args.command=='verify':
            data=verify_manifest(args.out)
            print(json.dumps(data,ensure_ascii=False,indent=2))
            if not data['ok']:return 1
        elif args.command=='overlap':
            result=evaluate_overlap(args.s,args.m2,args.order)
            output=json.dumps(result,ensure_ascii=False,indent=2)
            if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(output,encoding='utf-8')
            else:print(output)
        elif args.command=='two-loop':
            result=evaluate_two_loop(args.graph,args.invariants,args.m2,args.order)
            output=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)
            if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(output,encoding='utf-8')
            else:print(output)
        elif args.command=='serve':
            from .server import serve
            serve(Path(__file__).resolve().parents[1],args.port,not args.no_open)
        elif args.command=='analyze':
            data=analyze(Graph.from_dict(json.loads(args.input.read_text(encoding='utf-8'))),'input',args.dimension)
            s=json.dumps(data,ensure_ascii=False,indent=2)
            if args.out:args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(s,encoding='utf-8')
            else:print(s)
        else:print(json.dumps(demo(),indent=2))
    except (ValueError,KeyError,OSError,RuntimeError) as exc:
        print(f"Error: {exc}",file=sys.stderr);return 2
    return 0


if __name__=='__main__':raise SystemExit(run())
