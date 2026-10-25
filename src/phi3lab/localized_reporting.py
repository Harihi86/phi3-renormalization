"""Single-language workbenches and expression captions for Chinese/English releases.

The numerical graph programs are unchanged. Locale selection affects only
human-facing text, navigation, and generated TeX captions.
"""
from html import escape
import hashlib
import json
from pathlib import Path
import re

STRINGS = {
    'cn': {
        'lang': 'zh-CN', 'title': '从图构造 BPH 有限积分',
        'eyebrow': '图生成 · BPH 递归 · 六维有质量标量场',
        'intro': '由边、顶点和固定外腿标签识别发散子图，按博戈柳博夫递归构造相容动量路由、泰勒系数和高斯积分，再解析积分整体尺度。生产计算不按图号选择预写的拓扑公式。',
        'manual': '完整推导与验证', 'coverage': '本轮积分及独立核对记录', 'catalog': '全部表达式目录',
        'compiled': '已编译的图', 'integrated': '本轮实际完整求积的代表图', 'independent': '本轮有独立核对的代表图',
        'mapped': '由已求积代表及映射证书覆盖的图', 'zero': '由完整一点条件严格为零的图', 'representation_scope': '代表图实际执行完整积分；其余图通过顶点、内边、外腿标签与运动学保持的精确映射传播结果。映射覆盖数包含代表图，不能视为逐张独立求积。',
        'scope': '表达式覆盖、完整求积、独立方法核对与严格解析证书分开计数。随机组标准误和采样细化差是数值诊断，不能替代严格总误差界。未核验的项目内 Z 部分不作前提。',
        'development': '开发记录：本轮新结果尚未完成发布核查。', 'verified': '发布记录：请以逐图记录中的实际输入和核查状态为准。',
        'choose': '选择生成图', 'loops': '圈数', 'legs': '外腿数', 'graph': '本组图',
        'json': '下载可执行表达式', 'tex': '下载表达式排版源文件', 'construction': '查看递归构造',
        'construction_help': '表达式保存每项的整数路由、边界缩放幂次、泰勒次数和递归符号。完整递归和在同一原始参数点合并后求积。',
        'facts': '{id} · {vertices} 个顶点 · {edges} 条内边 · {parts} 个真发散子图 · {terms} 个递归项 · {dimension} 维参数积分',
        'plan': '已实际生成 {terms} 个递归项，最多 {jets} 个独立泰勒系数变量。等参数处的完整有限密度为 {density}。',
        'calculate': '在本机计算', 'mass': '共同质量平方', 'power': '每组样本数', 'reps': '独立随机组数', 'sampler': '采样方法',
        'auto': '自动选择', 'hepp': '按边尺度排序采样', 'dirichlet': '狄利克雷采样', 'recursive': '按递归结构采样', 'program': '按递归程序采样', 'mixture': '多通道混合采样',
        'compile': '重新构造并检查表达式', 'integrate': '运行完整递归积分', 'idle': '尚未运行当前选择。',
        'ready': '当前显示的是已生成的图和表达式记录。', 'compiling': '正在从图重新构造表达式……',
        'integrating': '正在计算完整递归和，请保持页面开启……', 'compiled_ok': '图导出表达式及密度检查已完成。',
        'integrated_ok': '完整递归积分已执行；报告误差为数值诊断。', 'file_error': '请通过启动文件开启本地服务后计算。',
        'error': '计算未完成。请核对输入并查看原始诊断记录。', 'diagnostic': '下载原始诊断记录',
        'basis': 'BPH 递归 · 由输入图构造', 'parts': '真发散子图', 'terms': '递归项', 'density': '等参数有限密度',
        'value': '有限积分', 'se': '独立随机组标准误', 'refinement': '半样本细化差', 'samples': '实际采样点数',
        'precision': '多精度复算点数', 'seconds': '计算用时（秒）', 'kernel': '本入口未调用预写的拓扑积分公式。',
        'precision_help': '程序监测子图尺度及递归系数抵消；需要时以多精度重新计算完整递归和。可选原生执行器缺失时使用标准库高精度十进制路径。库源码和二进制身份均需匹配。',
        'progress': '本轮完整积分覆盖', 'category': '子发散结构', 'total': '该类图总数', 'executed': '已完整求积', 'checked': '已独立核对',
        'nested': '嵌套子发散', 'overlap': '重叠子发散', 'disjoint': '不交子发散', 'primitive': '无真子发散', 'mixed': '同时含嵌套与重叠', 'single': '单个真子发散',
        'pending': '等待实际执行记录', 'graph_image': '当前图的内边、顶点及固定外腿',
    },
    'en': {
        'lang': 'en', 'title': 'Construct finite BPH integrals from graphs',
        'eyebrow': 'GRAPH GENERATION · BPH RECURSION · MASSIVE SIX-DIMENSIONAL SCALAR THEORY',
        'intro': 'Edges, vertices and fixed external-leg labels determine divergent subgraphs, compatible momentum routings, Taylor coefficients and Gaussian integrals through the Bogoliubov recursion. The common overall scale is integrated analytically. Production execution does not select a stored topology formula by graph number.',
        'manual': 'Derivation and validation', 'coverage': 'Current complete-integral and independent checks', 'catalog': 'Complete expression catalogue',
        'compiled': 'Compiled graphs', 'integrated': 'Representatives actually integrated in this release', 'independent': 'Representatives independently checked in this release',
        'mapped': 'Graphs covered by integrated representatives and mapping certificates', 'zero': 'Graphs identically zero by the full one-point condition', 'representation_scope': 'Representatives undergo actual complete integration. Results for other graphs are propagated through exact vertex, edge and external-label maps preserving the kinematics. Mapped coverage includes the representatives and does not count separate integrations of every graph.',
        'scope': 'Expression coverage, complete integration, independent methods and rigorous analytic certificates are counted separately. Between-scramble standard errors and refinement differences are numerical diagnostics, not rigorous total error bounds. The unverified Z portion of this project is not a premise.',
        'development': 'Development record: new results have not yet completed release verification.', 'verified': 'Release record: consult each graph record for its actual input and verification status.',
        'choose': 'Select a generated graph', 'loops': 'Loops', 'legs': 'External legs', 'graph': 'Graph in this group',
        'json': 'Download executable expression', 'tex': 'Download typesetting source', 'construction': 'Inspect the recursive construction',
        'construction_help': 'The expression retains integer routing rows, boundary-scaling powers, Taylor orders and recursive signs for every term. The complete recursive sum is combined at each common original-edge parameter point before integration.',
        'facts': '{id} · {vertices} vertices · {edges} internal edges · {parts} proper divergent subgraphs · {terms} recursive terms · {dimension} parameter dimensions',
        'plan': 'The program generated {terms} recursive terms with at most {jets} independent Taylor-coefficient variables. The complete finite density at equal parameters is {density}.',
        'calculate': 'Calculate locally', 'mass': 'Common squared mass', 'power': 'Samples per scramble', 'reps': 'Independent scrambles', 'sampler': 'Sampling method',
        'auto': 'Automatic selection', 'hepp': 'Edge-ordering sectors', 'dirichlet': 'Dirichlet', 'recursive': 'Recursive structure', 'program': 'Recursive program', 'mixture': 'Multiple-channel mixture',
        'compile': 'Reconstruct and check the expression', 'integrate': 'Integrate the complete recursive sum', 'idle': 'The current selection has not been executed.',
        'ready': 'The generated graph and expression records are shown.', 'compiling': 'Reconstructing the expression from the graph...',
        'integrating': 'Evaluating the complete recursive sum; keep this page open...', 'compiled_ok': 'Graph-derived expression and density checks completed.',
        'integrated_ok': 'Complete recursive integration executed; reported errors are numerical diagnostics.', 'file_error': 'Start the local service with the launcher before calculating.',
        'error': 'The calculation did not complete. Check the input and inspect the raw diagnostic record.', 'diagnostic': 'Download the raw diagnostic record',
        'basis': 'BPH recursion · constructed from the input graph', 'parts': 'Proper divergent subgraphs', 'terms': 'Recursive terms', 'density': 'Equal-parameter finite density',
        'value': 'Finite integral', 'se': 'Between-scramble standard error', 'refinement': 'Half-power refinement difference', 'samples': 'Actual sample count',
        'precision': 'Multiprecision recomputation points', 'seconds': 'Elapsed seconds', 'kernel': 'This entry point did not call a stored topology-specific integral formula.',
        'precision_help': 'The program monitors subgraph scales and coefficient cancellation and recomputes the complete recursive sum in multiprecision when needed. An independent standard-library Decimal path remains available if the optional native executor is unavailable. Source and library identities must match.',
        'progress': 'Current complete-integration coverage', 'category': 'Subdivergence structure', 'total': 'Graphs in this class', 'executed': 'Completely integrated', 'checked': 'Independently checked',
        'nested': 'Nested subdivergences', 'overlap': 'Overlapping subdivergences', 'disjoint': 'Disjoint subdivergences', 'primitive': 'Primitive', 'mixed': 'Nested and overlapping', 'single': 'Single proper subdivergence',
        'pending': 'Awaiting actual execution records', 'graph_image': 'Internal edges, vertices and fixed external legs of the selected graph',
    },
}


def _locale(language):
    if language in ('cn', 'zh', 'zh-CN'):
        return 'cn'
    if language == 'en':
        return 'en'
    raise ValueError('Use cn or en')


def localize_expression_tex(text, language='cn'):
    """Translate only human captions; mathematical program text is preserved."""
    if _locale(language) == 'en':
        return text
    for old, new in (
        (': graph-derived recursive BPH density', '：图导出的递归 BPH 密度'),
        ('Recursive word ', '递归项 '), (', sign ', '，符号 '),
        ('Inner-first embedded parts: ', '由内到外的嵌入子图：'),
        ('. Jet orders: ', '。泰勒系数最高次数：'),
        ('by the declared full one-point Taylor condition', '由所声明的完整一点泰勒条件'),
    ):
        text = text.replace(old, new)
    return text


def localize_expression_tree(folder, language='cn'):
    """Localize copied TeX assets and return fresh byte hashes for the manifest."""
    folder = Path(folder)
    result = {}
    for path in sorted(folder.glob('*.tex')):
        path.write_text(localize_expression_tex(path.read_text(encoding='utf-8'), language), encoding='utf-8')
        result[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def write_localized_html(data, path, *, language='cn', publication=None,
                         manual_href='../../manual/main.pdf',
                         coverage_href='../recursive_integral_comparison.json',
                         catalog_href='catalog.json'):
    """Write a single-language workbench; publication counts are actual run data.

    publication: version, verified, representative_integrated_count (or legacy
    complete_graph_count), propagated_graph_count, exact_zero_count,
    independent_graph_count, quality_rows[{loops, category, total, executed,
    cross_checked}], samplers. Mapped coverage includes its representatives;
    it never means every mapped graph was separately integrated.
    Missing new-run counts are displayed as zero,
    never inferred from the number of compiled expressions.
    """
    language = _locale(language)
    t = STRINGS[language]
    publication = dict(publication or {})
    samplers = publication.get('samplers', ['auto', 'recursive', 'mixture', 'hepp', 'dirichlet'])
    payload = json.dumps({'catalog': data, 'publication': publication, 'text': t}, ensure_ascii=False).replace('</', '<\\/')
    def tr(key): return escape(t[key])
    def metric(value, label): return f'<div class="metric"><b>{value}</b>{tr(label)}</div>'
    options = ''.join(f'<option value="{escape(value,quote=True)}">{tr(value) if value in t else escape(value)}</option>' for value in samplers)
    quality = publication.get('quality_rows', [])
    rows = ''.join('<tr>'+''.join(f'<td>{escape(str(v))}</td>' for v in (r['loops'], t.get(r['category'], r['category']), r['total'], r['executed'], r['cross_checked']))+'</tr>' for r in quality)
    if not rows: rows=f'<tr><td colspan="5">{tr("pending")}</td></tr>'
    h = r'''<!doctype html><html lang="LANG"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TITLE</title>
<style>*{box-sizing:border-box}body{margin:0;background:#f2f6f3;color:#183d38;font:16px/1.75 system-ui}main{max-width:1100px;margin:32px auto;padding:0 22px 64px}h1{font-size:32px;line-height:1.35}h2{font-size:23px}a{color:#126f60}.card{background:#fff;border:1px solid #d4e2da;border-radius:14px;padding:24px;margin:20px 0}.eyebrow,.small{font-size:13px;color:#567369}.grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}.metric{background:#edf4ee;padding:14px;border-radius:8px}.metric b{font-size:25px;display:block}label{display:block;font-size:14px}input,select,button{width:100%;padding:10px;border:1px solid #b6cec1;border-radius:7px;background:#fff;color:inherit;font:inherit}button{background:#1b6655;color:#fff;cursor:pointer;line-height:1.4}button.secondary{background:#eef4ef;color:#195a49}button:disabled{opacity:.5;cursor:wait}.buttons{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:18px}pre{font:14px/1.9 ui-monospace,monospace;white-space:pre-wrap;overflow-wrap:anywhere;background:#edf4ee;border-radius:8px;padding:18px}summary{cursor:pointer}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:9px;border-bottom:1px solid #d9e5dd}.scroll{overflow:auto}svg{display:block;margin:auto;width:100%;max-width:460px;height:260px}#status{padding-top:16px}nav{line-height:2.1}@media(max-width:700px){main{padding:0 15px 45px;margin:22px auto}.grid,.buttons{grid-template-columns:1fr}.card{padding:17px}h1{font-size:26px}}
</style><main><div class="eyebrow">PHI3LAB VERSION · EYEBROW</div><h1>TITLE</h1><nav><a href="MANUAL_HREF">MANUAL</a> · <a href="COVERAGE_HREF">COVERAGE</a> · <a href="CATALOG_HREF">CATALOG</a></nav>
<section class="card"><p>INTRO</p><div class="grid">METRICS</div><p class="small">STAGE SCOPE</p><p class="small">REPRESENTATION_SCOPE</p></section>
<section class="card"><h2>CHOOSE</h2><div class="grid"><label>LOOPS<select id="loops"><option>1</option><option>2</option><option selected>3</option><option>4</option></select></label><label>LEGS<select id="legs"><option>1</option><option>2</option><option selected>3</option></select></label><label>GRAPH<select id="graph"></select></label></div><div id="diagram" role="img" aria-label="GRAPH_IMAGE"></div><p id="facts"></p><p><a id="expression">JSON_LABEL</a> · <a id="tex">TEX_LABEL</a></p><details><summary>CONSTRUCTION</summary><p id="plan"></p><p class="small">CONSTRUCTION_HELP</p></details></section>
<section class="card"><h2>CALCULATE</h2><div class="grid"><label>MASS<input id="m2" type="number" value="1" min="0" step="any"></label><label>p₁²<input id="s1" type="number" value="1" min="0" step="any"></label><label>p₂²<input id="s2" type="number" value="1" min="0" step="any"></label><label>p₃²<input id="s3" type="number" value="1" min="0" step="any"></label><label>POWER<select id="power"><option value="12">4096</option><option value="14" selected>16384</option><option value="16">65536</option><option value="18">262144</option></select></label><label>REPS<select id="reps"><option>4</option><option selected>8</option><option>16</option></select></label><label>SAMPLER<select id="sampler">SAMPLER_OPTIONS</select></label></div><div class="buttons"><button id="compile" class="secondary">COMPILE</button><button id="integrate">INTEGRATE</button></div><p id="status" role="status">READY</p><pre id="numeric">IDLE</pre><p><a id="diagnostic" hidden>DIAGNOSTIC</a></p><p class="small">PRECISION_HELP</p></section>
<section class="card"><h2>PROGRESS</h2><div class="scroll"><table><thead><tr><th>LOOPS</th><th>CATEGORY</th><th>TOTAL</th><th>EXECUTED</th><th>CHECKED</th></tr></thead><tbody>QUALITY_ROWS</tbody></table></div><p class="small">SCOPE</p></section></main>
<script type="application/json" id="payload">PAYLOAD</script><script>
const bundle=JSON.parse(document.getElementById('payload').textContent),data=bundle.catalog,t=bundle.text,$=id=>document.getElementById(id);let row,serial=0,diagnosticURL=null;
const format=(key,values)=>t[key].replace(/\{(\w+)\}/g,(_,k)=>values[k]);
function update(){const rows=data.graphs.filter(g=>g.loops===Number($('loops').value)&&g.legs===Number($('legs').value));$('graph').replaceChildren(...rows.map(g=>{const o=document.createElement('option');o.value=g.id;o.textContent=g.id;return o}));select();}
function select(){serial++;row=data.graphs.find(g=>g.id===$('graph').value);$('facts').textContent=format('facts',{id:row.id,vertices:row.vertices,edges:row.edges,parts:row.proper_uv_parts,terms:row.recursive_terms,dimension:row.parameter_dimension});$('expression').href=row.expression;$('tex').href=row.tex;$('plan').textContent=format('plan',{terms:row.recursive_terms,jets:row.maximum_jet_dimension,density:row.equal_parameter_density});$('s1').disabled=row.legs===1;$('s2').disabled=$('s3').disabled=row.legs!==3;$('numeric').textContent=t.idle;$('status').textContent=t.ready;$('diagnostic').hidden=true;draw(row.graph);}
function draw(g){const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','0 0 400 290');const p=Array.from({length:g.vertices},(_,i)=>[200+96*Math.cos(2*Math.PI*i/g.vertices),145+96*Math.sin(2*Math.PI*i/g.vertices)]),seen={};for(const[a,b]of g.edges){const key=[a,b].join(','),k=seen[key]||0;seen[key]=k+1;const e=document.createElementNS(ns,'path'),x=p[a],y=p[b],mx=(x[0]+y[0])/2,my=(x[1]+y[1])/2,dx=y[0]-x[0],dy=y[1]-x[1],r=Math.hypot(dx,dy)||1,offset=k?16:-8;e.setAttribute('d',a===b?`M${x[0]},${x[1]}c-35,-55 35,-55 0,0`:`M${x}Q${mx-dy/r*offset},${my+dx/r*offset} ${y}`);e.setAttribute('stroke','#34776a');e.setAttribute('stroke-width','2');e.setAttribute('fill','none');svg.append(e);}p.forEach((v,i)=>{const c=document.createElementNS(ns,'circle');c.setAttribute('cx',v[0]);c.setAttribute('cy',v[1]);c.setAttribute('r','5');c.setAttribute('fill','#164e43');svg.append(c);const text=document.createElementNS(ns,'text');text.setAttribute('x',v[0]+8);text.setAttribute('y',v[1]-8);text.textContent=i;svg.append(text);});for(const[label,v]of g.external){const x=p[v],dx=x[0]-200,dy=x[1]-145,r=Math.hypot(dx,dy)||1,text=document.createElementNS(ns,'text');text.setAttribute('x',x[0]+dx/r*30-4);text.setAttribute('y',x[1]+dy/r*30+4);text.textContent='p'+label;svg.append(text);}$('diagram').replaceChildren(svg);}
function saveDiagnostic(record){if(diagnosticURL)URL.revokeObjectURL(diagnosticURL);diagnosticURL=URL.createObjectURL(new Blob([JSON.stringify(record,null,2)],{type:'application/json'}));$('diagnostic').href=diagnosticURL;$('diagnostic').download='diagnostic.json';$('diagnostic').hidden=false;}
async function run(integrate){const token=++serial,chosen=row,query=new URLSearchParams({loops:chosen.loops,legs:chosen.legs,index:chosen.index,m2:$('m2').value,integrate:integrate?'1':'0',power:$('power').value,replicates:$('reps').value,sampler:$('sampler').value});for(let j=1;j<=(chosen.legs===1?0:chosen.legs===2?1:3);j++)query.set('s'+j,$('s'+j).value);$('compile').disabled=$('integrate').disabled=true;$('status').textContent=integrate?t.integrating:t.compiling;$('diagnostic').hidden=true;try{if(location.protocol==='file:')throw{localFile:true};const response=await fetch('/api/automatic?'+query),r=await response.json();if(!response.ok)throw{record:r};if(token!==serial)return;const e=r.expression;let lines=[r.graph_id,t.basis,`${t.parts} = ${e.proper_subgraphs.length}`,`${t.terms} = ${e.terms.length}`,`${t.density} = ${r.equal_parameter_density}`];if(r.integral){const x=r.integral,p=x.precision_control||{};lines.push('',`${t.value} = ${x.value}`,`${t.se} = ${x.standard_error}`,`${t.refinement} = ${x.power_refinement_difference}`,`${t.samples} = ${x.total_samples}`,`${t.precision} = ${p.precision_recomputed_points??p.decimal_recomputed_points??0}`,`${t.seconds} = ${Number(x.seconds).toFixed(2)}`);}lines.push('',t.kernel);$('numeric').textContent=lines.join('\n');$('status').textContent=integrate?t.integrated_ok:t.compiled_ok;saveDiagnostic(r);}catch(error){if(token===serial){$('status').textContent=error.localFile?t.file_error:t.error;saveDiagnostic(error.record||{message:String(error)});}}finally{$('compile').disabled=$('integrate').disabled=false;}}
$('loops').onchange=$('legs').onchange=update;$('graph').onchange=select;$('compile').onclick=()=>run(false);$('integrate').onclick=()=>run(true);update();
</script></html>'''
    replacements = {
        'LANG': t['lang'], 'VERSION': str(publication.get('version', '4.1')), 'TITLE': t['title'],
        'MANUAL_HREF': manual_href, 'COVERAGE_HREF': coverage_href, 'CATALOG_HREF': catalog_href,
        'METRICS': metric(len(data['graphs']), 'compiled')+metric(int(publication.get('representative_integrated_count', publication.get('complete_graph_count', 0))), 'integrated')+metric(int(publication.get('propagated_graph_count', 0)), 'mapped')+metric(int(publication.get('independent_graph_count', 0)), 'independent')+metric(int(publication.get('exact_zero_count', 0)), 'zero'),
        'STAGE': t['verified' if publication.get('verified') is True else 'development'],
        'JSON_LABEL': t['json'], 'TEX_LABEL': t['tex'], 'SAMPLER_OPTIONS': options, 'QUALITY_ROWS': rows,
        'PAYLOAD': payload,
    }
    for key in ('eyebrow','intro','manual','coverage','catalog','scope','choose','loops','legs','graph','construction','construction_help','calculate','mass','power','reps','sampler','compile','integrate','idle','ready','diagnostic','precision_help','progress','category','total','executed','checked','graph_image','representation_scope'):
        replacements[key.upper()] = escape(t[key])
    # Tokens are replaced in one pass so translated text is never rewritten.
    h = re.sub(r'\b('+'|'.join(sorted(map(re.escape,replacements),key=len,reverse=True))+r')\b',lambda match:replacements[match.group()],h)
    path = Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(h,encoding='utf-8')
    return {'language': language, 'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def _tex_escape(value):
    replacements = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    return ''.join(replacements.get(c, c) for c in str(value))


def write_localized_manual_data(publication, path, *, language='cn'):
    """Create the language-specific TeX input from frozen actual release records.

    Shared numeric fields: version, date, verified, source_fingerprint,
    catalogue_count, representative_integrated_count, propagated_graph_count,
    independent_graph_count, exact_zero_count. quality_rows use the HTML schema.
    Optional numerical_rows contain id, value, standard_error, difference,
    method (production|independent|paired), status (pass|pending|failed).
    Caller supplies counts from executed records, never from intended targets.
    """
    language = _locale(language)
    t = STRINGS[language]
    p = dict(publication or {})
    cn = language == 'cn'
    text = {
        'status': (('本轮发布核查已完成；具体精度仍以逐图诊断为准。' if p.get('verified') is True else '本轮结果尚未完成发布核查；下列数值只从实际执行记录导入。') if cn else ('Release verification completed; numerical precision remains subject to each graph record.' if p.get('verified') is True else 'Release verification is incomplete; the values below are imported only from actual execution records.')),
        'missing': '尚未冻结' if cn else 'Not frozen',
        'pending': '尚无已导入的本轮执行记录' if cn else 'No current execution records imported',
        'production': '生产递归核' if cn else 'Production recursive kernel',
        'independent': '独立多项式核' if cn else 'Independent polynomial kernel',
        'paired': '独立方法配对比较' if cn else 'Paired independent comparison',
        'pass': '核查通过' if cn else 'Passed',
        'failed': '需继续修复' if cn else 'Requires further repair',
    }
    macros = {
        'ReleaseVersion': p.get('version', '4.1'),
        'ReleaseDate': p.get('date', text['missing']),
        'ReleaseStatus': text['status'],
        'ReleaseFingerprint': p.get('source_fingerprint', text['missing']),
        'ReleaseCatalogueCount': int(p.get('catalogue_count', 698)),
        'ReleaseActualCount': int(p.get('representative_integrated_count', p.get('complete_graph_count', 0))),
        'ReleaseMappedCount': int(p.get('propagated_graph_count', 0)),
        'ReleaseIndependentCount': int(p.get('independent_graph_count', 0)),
        'ReleaseZeroCount': int(p.get('exact_zero_count', 0)),
    }
    lines = ['% Generated from frozen execution records; do not infer counts from targets.']
    lines += ['\\newcommand{\\'+key+'}{'+_tex_escape(value)+'}' for key,value in macros.items()]
    rows = p.get('quality_rows', [])
    quality = []
    for r in rows:
        values = [int(r['loops']), t[r['category']], int(r['total']), int(r['executed']), int(r['cross_checked'])]
        quality.append(' & '.join(_tex_escape(x) for x in values)+r' \\')
    if not quality:
        quality = [r'\multicolumn{5}{l}{'+text['pending']+r'} \\']
    lines.append(r'\newcommand{\ReleaseQualityRows}{'+'\n'.join(quality)+'}')
    numerical = []
    for r in p.get('numerical_rows', []):
        values = [r['id'], r.get('value','--'), r.get('standard_error','--'), r.get('difference','--'), text[r.get('method','production')], text.get(r.get('status','pending'),text['pending'])]
        numerical.append(' & '.join(_tex_escape(x) for x in values)+r' \\')
    if not numerical:
        numerical = [r'\multicolumn{6}{l}{'+text['pending']+r'} \\']
    lines.append(r'\newcommand{\ReleaseNumericalRows}{'+'\n'.join(numerical)+'}')
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return {'language': language, 'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
