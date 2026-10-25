"""Freshness-aware evidence and reproducible build manifests."""
import hashlib,json,platform
from datetime import datetime,timezone
from pathlib import Path
from . import __version__

ROOT = Path(__file__).resolve().parents[1]


def source_hashes():
    paths=list((ROOT/'phi3lab').glob('*.py'))+list((ROOT/'native').glob('*.cpp'))
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def source_fingerprint():
    return hashlib.sha256(json.dumps(source_hashes(),sort_keys=True).encode()).hexdigest()


def supporting_hashes():
    paths=[ROOT/'README.md',ROOT/'NOTICE.md',ROOT/'启动分析.command',ROOT/'启动分析.bat',ROOT/'requirements-numerics.txt']
    for pattern in ('docs/*','tools/*.py','tests/*.py','native/*.md','native/*.json','native/build/*.json','native/build/*.dylib','native/build/*.so','analysis-project/*.tex','analysis-project/*.pdf','full-analysis-project/*.tex','full-analysis-project/*.pdf','extreme-analysis-project/*.tex','extreme-analysis-project/*.pdf','automatic-analysis-project/*.tex','automatic-analysis-project/*.pdf'):
        paths.extend(ROOT.glob(pattern))
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(set(paths)) if p.is_file()}


def external_evidence(folder):
    result=[]
    for name in ('external_comparison.json','pysecdec_comparison.json','test_verification.json','overlap_comparison.json','nested_triangle_comparison.json','native_primitive_comparison.json','extreme_triangle_comparison.json','extreme_native_comparison.json','extreme_bubble_comparison.json','extreme_overlap_comparison.json','higher_graph_feyngen_comparison.json','automatic_bph_comparison.json','automatic_integral_comparison.json','automatic_higher_integrals.json'):
        path=Path(folder)/name
        if not path.exists():
            result.append({'file':name,'status':'not_run','description':'尚未运行'})
            continue
        try:
            evidence=json.loads(path.read_text(encoding='utf-8'))
            fresh=evidence.get('source_fingerprint')==source_fingerprint()
            ok=evidence.get('all_ok') is True
            status='passed' if fresh and ok else 'failed' if fresh else 'stale'
            desc={'passed':'当前源码：已通过','failed':'当前源码：检查失败','stale':'旧版本证据：需重新核查'}[status]
            result.append({'file':name,'status':status,'description':desc,
                           'version':evidence.get('version',evidence.get('commit',''))})
        except (ValueError,OSError,TypeError):
            result.append({'file':name,'status':'invalid','description':'无法读取证据文件'})
    return result


def write_manifest(folder, files, data):
    folder=Path(folder).resolve()
    hashes={str(Path(p).resolve().relative_to(folder)):hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in sorted(map(Path,files))}
    manifest={'schema_version':2,'software_version':__version__,
              'generated_utc':datetime.now(timezone.utc).isoformat(),
              'python':platform.python_version(),'machine':platform.machine(),
              'source_fingerprint':source_fingerprint(),'source_hashes':source_hashes(),
              'supporting_hashes':supporting_hashes(),
              'artifact_hashes':hashes,
              'graph_count':len(data['graphs']),
              'dimension':data['dimension'],
              'wick_cases':len(data['wick_validation']),
              'external_evidence':external_evidence(folder),
              'limits':'Hashes establish file consistency, not mathematical correctness. Unit tests and numerical integration are not implied by build.'}
    (folder/'provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    return manifest


def verify_manifest(folder):
    folder=Path(folder).resolve()
    data=json.loads((folder/'provenance.json').read_text(encoding='utf-8'))
    problems=[]
    if data.get('source_fingerprint')!=source_fingerprint():problems.append('Source files changed since this build')
    for name,expected in data.get('supporting_hashes',{}).items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT):problems.append('Unsafe supporting path: '+name);continue
        if not path.is_file():problems.append('Missing supporting file: '+name)
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=expected:problems.append('Changed supporting file: '+name)
    for name,expected in data.get('artifact_hashes',{}).items():
        path=(folder/name).resolve()
        if not path.is_relative_to(folder):
            problems.append('Unsafe artifact path: '+name);continue
        if not path.is_file():problems.append('Missing: '+name)
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=expected:problems.append('Changed: '+name)
    if not data.get('artifact_hashes'):problems.append('No build artifact hashes recorded')
    return {'ok':not problems,'software_version':data.get('software_version'),
            'checked_artifacts':len(data.get('artifact_hashes',{})),'problems':problems,
            'external_evidence':external_evidence(folder)}
