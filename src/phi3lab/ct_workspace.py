"""Unified bilingual workbench over retained computation and CT services."""
from pathlib import Path
from fractions import Fraction
import json
from .graphs import generate as low_generate
from .higher_graphs import generate as high_generate,canonical
from .ct_identity import record
from .ct_layout import geometry,svg,tikz
from .graphs import Graph
from .reporting import svg as legacy_svg


def build(root):
    root=Path(root);out=root/'results/workspace';out.mkdir(parents=True,exist_ok=True)
    fig=out/'figures';fig.mkdir(exist_ok=True);tex=out/'latex';tex.mkdir(exist_ok=True)
    old=json.loads((root/'results/counterterms/correspondence.json').read_text(encoding='utf-8'))
    old_ids={r['old_id']:r['id'] for r in old['rows']};rows=[];families=[]
    registry=root/'results/series/families.json'
    locations=json.loads(registry.read_text(encoding='utf-8')) if registry.exists() else {}
    for loops in range(1,5):
        for legs in (1,2,3):
            graphs=(low_generate if loops<=2 else high_generate)(loops,legs)
            family_dir=locations.get(f'L{loops}E{legs}',f'L{loops}E{legs}')
            report_path=root/f'results/series/{family_dir}/report.json'
            report=json.loads(report_path.read_text(encoding='utf-8')) if report_path.exists() else None
            passed=report is not None and report['verification']['all_equal'] and len(report['rows'])==len(graphs)
            status='verified' if passed or (loops,legs)==(2,3) else 'pending'
            report_link=f'../series/{family_dir}/index-{{lang}}.html' if report else '../counterterms/index-{lang}.html' if (loops,legs)==(2,3) else None
            families.append({'loops':loops,'legs':legs,'graphs':len(graphs),'report':report_link,
                'ct_status':status,
                'numerical_status':'not_newly_accepted','linear_finite_condition':'unspecified' if legs==1 else 'not_applicable_to_1PI_hosts'})
            for index,g in enumerate(graphs,1):
                legacy=f'L{loops}E{legs}G{index}';r=record(g,legacy)
                if legacy in old_ids and r['id']!=old_ids[legacy]:raise ArithmeticError('CT1 identity changed')
                display=Graph.from_dict(r['canonical_graph']);geo=geometry(display)
                (fig/f'{legacy}.svg').write_text(svg(display,geo),encoding='utf-8')
                (fig/f'{legacy}-legacy.svg').write_text(legacy_svg(g),encoding='utf-8')
                (tex/f'{legacy}.tex').write_text(tikz(display,geo),encoding='utf-8')
                rows.append({**r,'loops':loops,'legs':legs,'index':index,'weight':str(Fraction(1,r['symmetry'])),
                    'layout':geo,'ct_status':status,'report':report_link,
                    'numerical_status':'not_newly_accepted',
                    'historical_expression':f'../../legacy/cn/results/automatic/expressions/{legacy}.json'})
    data={'schema':'phi3lab.unified-workspace.v1','rows':rows,'families':families,
          'scope':'fixed-label connected amputated 1PI; loops 1..4, external legs 1..3',
          'coverage_claim':'catalogue enumeration is separate from counterterm correspondence and numerical accuracy',
          'unresolved_precision_cases':90,'strict_analytic_cases':24,'empirical_cases':48}
    (out/'catalog.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    for lang in ('cn','en'):(out/f'index-{lang}.html').write_text(page(data,lang),encoding='utf-8')
    return {'graphs':len(rows),'families':len(families),'output':str(out)}


def page(data,lang):
    cn=lang=='cn';T=lambda a,b:a if cn else b
    labels={'graphs':T('图集与身份','Graphs and identity'),'expression':T('递归表达式','Recursive expression'),
            'counterterms':T('反项对应','Counterterm correspondence'),'numerical':T('数值计算','Numerical calculation'),
            'validation':T('验证记录','Validation records'),'verified':T('对应已通过','Correspondence verified'),
            'pending':T('对应待完成','Correspondence pending'),'all':T('全部状态','All states'),
            'empty':T('此筛选下没有图。','No graphs match this filter.'),
            'wait':T('计算中；本次结果不自动构成精度验收。','Computing; this result is not automatic accuracy acceptance.'),
            'offline':T('离线页面可查看图集；请运行启动程序以使用计算接口。','Offline catalogue is available. Run the launcher to use calculation endpoints.'),
            'pending_detail':T('此族的完整反项对应尚未验收。图集与历史表达式可查看；不隐藏未决图。','Full counterterm correspondence for this family is not accepted yet. The complete catalogue and historical expressions remain available.'),
            'numeric_note':T('保留原积分核。采样结果、精确组合验证和精度验收分别记录；162 组研究总账：严格解析 24、经验验收 48、待解决 90。','Retained integration kernels. Sampling, exact combinatorics and accuracy acceptance are separate; the 162-case research ledger has 24 strict analytic, 48 empirical and 90 unresolved cases.'),
            'identity_note':T('身份取决于结构、类型与固定外腿标签；切换布局不改变图。','Identity depends on structure, types and fixed external labels; changing layout does not change the graph.'),
            'linear_note':T('一点反项的有限线性条件未指定。已有一点数值方案不构成新物理条件的选择。','The finite linear counterterm condition is unspecified. The retained one-point numerical scheme does not select a new physical condition.'),
            'shared':T('共享布局','Shared layout'),'legacy':T('旧版布局','Legacy layout')}
    options=''.join(f'<option value="{k}">{labels[k]}</option>' for k in ('graphs','expression','counterterms','numerical','validation'))
    template='''<!doctype html><html lang="LANG"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>TITLE</title>
<style>body{max-width:1200px;margin:30px auto;padding:0 22px;background:#f3f6f5;color:#163b45;font:16px/1.6 system-ui}h1{font-size:30px}nav,section{background:white;border:1px solid #dbe4e0;border-radius:12px;padding:18px;margin:18px 0}label{display:inline-block;margin:5px 14px 5px 0}select,input,button{font:inherit;padding:6px;border:1px solid #b8c9c3;border-radius:6px}button{cursor:pointer}a{color:#246b5b}#gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:14px}figure{margin:0}figure img{width:100%;max-height:280px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px;max-height:550px;overflow:auto}.badge{background:#e4efe9;padding:4px 12px;border-radius:20px}.small{font-size:13px;color:#536e68}#graph{max-width:250px}</style>
<h1>TITLE</h1><p>SUBTITLE</p><p><a href="index-OTHER.html">LANGSWITCH</a> · <a href="../counterterms/index-LANG.html">CT1</a> · <a href="../../legacy/LANG/results/automatic/index.html">LEGACY</a> · <a href="catalog.json">DATA</a></p>
<nav><label>LOOPS <select id="loops"><option>1</option><option selected>2</option><option>3</option><option>4</option></select></label><label>LEGS <select id="legs"><option>1</option><option>2</option><option selected>3</option></select></label><label>FUNCTION <select id="view">OPTIONS</select></label><label>STATUS <select id="status"><option value="all">ALLSTATUS</option><option value="verified">VERIFIED</option><option value="pending">PENDING</option></select></label><label>GRAPH <select id="graph"></select></label><p id="summary"></p></nav>
<section><div id="gallery"></div><p id="identity" class="small"></p><p id="links"></p><div id="panel"></div></section>
<p class="small">BOUNDARY</p><script>const DATA=PAYLOAD,L=LABELS,LANG='LANG';const $=s=>document.getElementById(s);let selected=[];
function refill(){selected=DATA.rows.filter(r=>r.loops==$('loops').value&&r.legs==$('legs').value&&($('status').value==='all'||r.ct_status===$('status').value));$('graph').replaceChildren(...selected.map(r=>{const o=document.createElement('option');o.value=r.legacy_id;o.textContent=r.legacy_id+' · 1/S='+r.weight;return o}));render()}
function render(){const r=selected.find(r=>r.legacy_id===$('graph').value);if(!r){$('summary').textContent=L.empty;$('gallery').innerHTML='';$('panel').innerHTML='';return}
$('summary').textContent=selected.length+' '+L.graphs+' · '+L[r.ct_status];$('gallery').innerHTML='<figure><img src="figures/'+r.legacy_id+'.svg"><figcaption>'+L.shared+'</figcaption></figure><figure><img src="figures/'+r.legacy_id+'-legacy.svg"><figcaption>'+L.legacy+'</figcaption></figure>';
$('identity').textContent=L.identity_note+' '+r.legacy_id+' → '+r.id;
$('links').innerHTML='<a href="figures/'+r.legacy_id+'.svg">SVG</a> · <a href="latex/'+r.legacy_id+'.tex">TikZ / LaTeX</a> · <a href="../../legacy/'+LANG+'/results/automatic/expressions/'+r.legacy_id+'.json">'+L.expression+'</a>';
const v=$('view').value;$('panel').replaceChildren();if(r.legs===1){const p=document.createElement('p');p.textContent=L.linear_note;$('panel').append(p)}
if(v==='graphs'){const pre=document.createElement('pre');pre.textContent=JSON.stringify({graph:r.graph,weight:r.weight,identity_schema:r.identity_schema,to_canonical:r.to_canonical},null,2);$('panel').append(pre)}
if(v==='counterterms'){if(r.report){const url=r.report.replace('{lang}',LANG);$('panel').innerHTML='<a href="'+url+'?graph='+r.legacy_id+'">'+r.legacy_id+' · '+L.counterterms+'</a>'}else{$('panel').textContent=L.pending_detail+(r.legs===1?' '+L.linear_note:'')}}
if(v==='expression'){$('panel').innerHTML='<button id="evaluate">'+L.expression+'</button><pre id="result"></pre>';$('evaluate').onclick=()=>evaluate(r,false)}
if(v==='numerical'){$('panel').innerHTML='<p>'+L.numeric_note+(r.legs===1?' '+L.linear_note:'')+'</p><label>m² <input id="m2" value="1" size="5"></label>'+Array.from({length:r.legs===3?3:r.legs===2?1:0},(_,i)=>'<label>s'+(i+1)+' <input id="s'+(i+1)+'" value="1" size="5"></label>').join('')+'<label>2^ <input id="power" value="8" size="3"></label><button id="evaluate">'+L.numerical+'</button><pre id="result"></pre>';$('evaluate').onclick=()=>evaluate(r,true)}
if(v==='validation'){$('panel').innerHTML='<p>'+L.numeric_note+'</p><a href="../../evidence/ct-next/baseline-tests/status.json">Baseline regression</a> · <a href="../counterterms/verification.json">CT1: 77 cells</a> · <a href="../../evidence/ct-next/baseline-ssh-comparison.json">Windows / Mac</a>'}}
async function evaluate(r,integrate){if(location.protocol==='file:'){$('result').textContent=L.offline;return}const q=new URLSearchParams({loops:r.loops,legs:r.legs,index:r.index,integrate:integrate?'1':'0'});if(integrate){q.set('m2',$('m2').value);q.set('power',$('power').value);q.set('replicates','2');for(let i=1;i<=3;i++)if($('s'+i))q.set('s'+i,$('s'+i).value)}$('result').textContent=L.wait;try{const response=await fetch('/api/automatic?'+q);$('result').textContent=JSON.stringify(await response.json(),null,2)}catch(e){$('result').textContent=String(e)}}
for(const id of ['loops','legs','status'])$(id).onchange=refill;for(const id of ['graph','view'])$(id).onchange=render;refill();</script></html>'''
    replacements={'TITLE':T('Phi³ Lab · 统一工作台','Phi³ Lab · Unified workbench'),
      'SUBTITLE':T('一至四圈 · 一至三条外腿 · 图生成、递归表达式、反项与数值计算','One to four loops · one to three external legs · graphs, recursive expressions, counterterms and numerics'),
      'CT1':T('两圈三点精确对应','Two-loop three-point correspondence'),'LEGACY':T('已有数值工作台','Retained numerical workbench'),
      'LANGSWITCH':'English' if cn else '中文','OTHER':'en' if cn else 'cn','LOOPS':T('圈数','Loops'),'LEGS':T('外腿数','External legs'),
      'FUNCTION':T('功能','Function'),'STATUS':T('验证状态','Verification status'),'GRAPH':T('图','Graph'),
      'ALLSTATUS':labels['all'],'VERIFIED':labels['verified'],'PENDING':labels['pending'],'OPTIONS':options,
      'BOUNDARY':labels['numeric_note'],'PAYLOAD':json.dumps(data,ensure_ascii=False).replace('</','<\\/'),
      'LABELS':json.dumps(labels,ensure_ascii=False),'DATALINK':T('完整目录 JSON','Complete catalogue JSON')}
    # Replace placeholders once, without interpreting replacement content.
    import re
    replacements['LANG']=lang
    template=template.replace('>DATA</a>','>DATALINK</a>')
    return re.sub(r'\b(?:'+ '|'.join(sorted(replacements,key=len,reverse=True))+r')\b',lambda m:replacements[m[0]],template)
