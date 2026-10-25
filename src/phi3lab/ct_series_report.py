"""Complete bilingual family report, with row-wise evidence and shared geometry.

The browser loads one row at a time. Every matrix cell is exported to CSV;
zero cells are represented explicitly in that export and in verification.json.
The UI's sparse row payload only compresses verified structural zeros.
"""
from pathlib import Path
import csv,html,json
from .graphs import Graph,subgraph
from .ct_identity import record
from .ct_layout import geometry,layout,anchors,svg,tikz


def write_json(path,data):
    path.write_text(json.dumps(data,ensure_ascii=False,separators=(',',':')),encoding='utf-8')


def figures_for_row(row,columns,out):
    """All coordinates descend from the parent, with explicit maps saved."""
    g=Graph.from_dict(row['graph']);pos=layout(g);name=row['old_id'];records={}
    def emit(suffix,graph,positions=None,marked=(),selected=None,bridges=(),provenance=None,external=None):
        stem=name+'-'+suffix;geo=geometry(graph,positions,external_positions=external,provenance=provenance)
        records[suffix]={'svg':'figures/'+stem+'.svg','tex':'latex/'+stem+'.tex','graph':graph.to_dict(),'geometry':geo}
        (out/'figures'/f'{stem}.svg').write_text(svg(graph,geo,marked,selected,bridges),encoding='utf-8')
        (out/'latex'/f'{stem}.tex').write_text(tikz(graph,geo,marked,selected,bridges),encoding='utf-8')
        return suffix
    emit('parent',g,pos)
    for i,s in enumerate(row['recursive_ast']['proper_parts']):
        emit(f'sub{i}',g,pos,s['edge_ids'],provenance={'original_edges':s['edge_ids']})
        sg,ports=subgraph(g,s['edge_ids']);verts=s['vertices'];external={}
        for j,(kind,index,end,v) in enumerate(ports,1):
            if kind=='external':external[j]=anchors(g)[index]
            else:
                other=g.edges[index][1-end]
                external[j]=tuple((pos[v][k]+pos[other][k])/2 for k in (0,1))
        emit(f'source{i}',sg,{j:pos[v] for j,v in enumerate(verts)},external=external,
             provenance={'source_to_parent_vertices':dict(enumerate(verts)),'source_to_parent_edges':dict(enumerate(s['edge_ids'])),'boundary_ports':ports})
    for i,t in enumerate(row['terms']):
        if t['kind']!='insertion':continue
        # Each decorated forest is a term; the top spinney defines its host.
        top=t['top_subgraphs'];mapping=t['contraction_mapping'];groups=[s['vertices'] for s in top]
        used={v for group in groups for v in group};groups.extend([v] for v in range(g.n) if v not in used)
        host=Graph.from_dict(columns[t['column']]['host']);vm={v:mapping['vertices'][j] for j,group in enumerate(groups) for v in group}
        hp={mapping['vertices'][j]:tuple(sum(pos[v][k] for v in group)/len(group) for k in (0,1)) for j,group in enumerate(groups)}
        provenance={'parent':row['id'],'parent_to_host_vertices':vm,'contraction_mapping':mapping,'top_subgraphs':top,'all_projectors':t['nodes']}
        emit(f't{i}-forest',g,pos,{e for n in t['nodes'] for e in n['edges']},provenance=provenance)
        emit(f't{i}-host',host,hp,provenance=provenance)
    for h in row.get('henges') or []:
        index=h['minimum_edge'];vm={int(k):v for k,v in h['original_vertex_to_skeleton'].items()}
        marked={e for block in h['blocks'] if block['loops'] for e in block['edge_ids']}
        emit(f'henge-e{index}',g,pos,marked,index,h['bridges'],{'original_edge_map':h['original_edge_map']})
        sg=Graph(len(h['blocks']),tuple(tuple(e) for e in h['skeleton_edges']),tuple((l,vm[v]) for l,v in g.external),
                 tuple('H' if b['loops'] else 'phi3' for b in h['blocks']))
        hp={j:tuple(sum(pos[v][k] for v in b['vertices'])/len(b['vertices']) for k in (0,1)) for j,b in enumerate(h['blocks'])}
        emit(f'skeleton-e{index}',sg,hp,provenance={'parent_to_skeleton_vertices':vm,'original_edge_map':h['original_edge_map']})
    return records


def write_reports(a,b,verification,out):
    out=Path(out)
    for folder in ('rows','figures','latex','layouts'):(out/folder).mkdir(exist_ok=True)
    columns={c['id']:c for c in a['columns']};ordered=['bare']+sorted(columns)+sorted(a['overall_columns'])
    column_meta=[]
    for i,cid in enumerate(ordered):
        if cid in columns:
            kinds=[k for k in columns[cid]['host']['vertex_kinds'] if k!='phi3'];label=f'H{i}: '+', '.join(kinds)
            host=Graph.from_dict(columns[cid]['host']);geo=geometry(host)
            (out/'figures'/f'column-{i}.svg').write_text(svg(host,geo),encoding='utf-8')
            (out/'latex'/f'column-{i}.tex').write_text(tikz(host,geo),encoding='utf-8')
        else:label=cid
        column_meta.append({'id':cid,'name':label,'figure':f'figures/column-{i}.svg' if cid in columns else None})
    brow={r['id']:r for r in b['rows']};byrow={r['id']:[] for r in a['rows']}
    for cell in verification['cells']:byrow.setdefault(cell['row'],[]).append(cell)
    rows=[]
    with (out/'table.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.writer(stream);writer.writerow(['old_id','graph_id','column','A_exact','B_exact','equal','structural_zero','zero_reason'])
        for row in a['rows']:
            cells=byrow[row['id']];ok=all(c['equal'] for c in cells)
            for c in cells:writer.writerow([row['old_id'],row['id'],c['column'],json.dumps(c['A']),json.dumps(c['B']),c['equal'],c['zero'],c['zero_reason']])
            figs=figures_for_row(row,columns,out);write_json(out/'layouts'/f'{row["old_id"]}.json',figs)
            payload={'A':row,'B':brow.get(row['id']),'nonzero_cells':[c for c in cells if not c['zero']],
                     'checked_cells':len(cells),'all_equal':ok,'zero_rule':'Any remaining catalogue column is a verified structural zero; see full table.csv.',
                     'figures':{k:{field:v[field] for field in ('svg','tex')} for k,v in figs.items()}}
            write_json(out/'rows'/f'{row["old_id"]}.json',payload)
            rows.append({'id':row['id'],'old_id':row['old_id'],'weight':row['positive_weight'],'all_equal':ok,
                         'terms':len(row['terms']),'data':'rows/'+row['old_id']+'.json','layouts':'layouts/'+row['old_id']+'.json',
                         'identity':record(Graph.from_dict(row['graph']),row['old_id'])})
    manifest={'schema':'phi3lab.ct-family-report.v1','loops':a['loops'],'legs':a['legs'],'rows':rows,'columns':column_meta,
              'verification':{k:v for k,v in verification.items() if k not in ('cells','differences')},
              'scope':'complete fixed-label connected amputated 1PI family; no graph dropped',
              'source_library':'path-b.json','numerical_accuracy_accepted':False,'unresolved_precision_cases':90,'strict_analytic_cases':24,'empirical_cases':48,
              'finite_linear_condition':'unspecified','reference_diagrams':'not supplied; external comparison pending',
              'taylor':'Two-point: total external-vector degree <=2, separate C2m/C2z. Three-point: degree zero.'}
    write_json(out/'report.json',manifest)
    for lang in ('cn','en'):
        (out/f'index-{lang}.html').write_text(page(manifest,lang),encoding='utf-8')
        latex_report(a,ordered,out,lang)
    return {'rows':len(rows),'columns':len(ordered),'cells':verification['checked_cells'],'all_equal':verification['all_equal']}


def latex_report(a,ordered,out,lang):
    cn=lang=='cn';T=lambda x,y:x if cn else y
    lines=[r'\documentclass[UTF8]{ctexart}' if cn else r'\documentclass{article}',
           r'\usepackage{amsmath,longtable,tikz,geometry}\geometry{margin=18mm}',r'\begin{document}',
           r'\section*{'+T('反项对应：','Counterterm correspondence: ')+f'L{a["loops"]}E{a["legs"]}'+'}',
           T('固定外腿标签；全部截肢连通 1PI 图。系数为精确有理数，积分保持符号形式。','Fixed external labels; complete amputated connected 1PI family. Exact rational coefficients; integrals remain symbolic.'),
           r'\[C(\gamma)=-t_\gamma\bar R(\gamma),\quad R(G)=(1-t_G)\bar R(G).\]',
           T('二点 Taylor 总向量次数不超过二；质量与动能分开。三点为零阶。未进行新增数值精度验收。','Two-point Taylor total vector degree is at most two, with separate mass and kinetic operators. Three-point degree is zero. No new numerical accuracy acceptance.')]
    for row in a['rows']:
        lines.extend([r'\subsection*{'+row['old_id']+'}',r'\input{'+row['old_id']+'-parent.tex}',
                      r'\begin{longtable}{lll}',r'Column & Coefficient & Forest \\\hline'])
        for i,t in enumerate(row['terms']):
            col='H'+str(ordered.index(t['column'])) if t['column'] not in ['bare']+a['overall_columns'] else t['column']
            nodes='; '.join(n['op']+'@'+','.join(map(str,n['edges'])) for n in t['nodes']) or 'empty'
            # A paragraph cell keeps long forest descriptions editable and inside margins.
            lines.append(col+' & $'+t['coefficient']+'$ & '+r'\parbox[t]{.58\textwidth}{\ttfamily\raggedright '+nodes+r'} \\')
        lines.extend([r'\end{longtable}',T('其余列为已核验的结构零；完整矩阵及两路径来源见 CSV 和 JSON。','All other columns are verified structural zeros; see CSV and JSON for the complete matrix and both paths.')])
    lines.append(r'\end{document}');(out/'latex'/f'report-{lang}.tex').write_text('\n'.join(lines),encoding='utf-8')


def page(data,lang):
    cn=lang=='cn';T=lambda x,y:x if cn else y
    labels={'bare':T('裸图','Bare graph'),'select':T('选择单元格','Select a cell'),'equal':T('精确相等','Exactly equal'),
      'zero':T('结构零：该原图与宿主之间两条路径均无贡献。','Structural zero: neither path contributes for this parent and host.'),
      'failed':T('不相等：保留失败记录','Mismatch: failure retained'),'trace':T('完整记录','Full trace'),'offline':T('请从启动程序打开以交互读取逐图记录；离线可直接查看 JSON、CSV 和 LaTeX。','Use the launcher for interactive row loading; JSON, CSV and LaTeX remain available offline.'),
      'a':T('A：裸图递归减除','A: recursive bare-graph subtraction'),'b':T('B：反项作用量生成','B: counterterm action generation'),
      'terms':T('逐项符号与系数','Symbolic terms and coefficients'),'nodes':T('投影算符与边集','Projectors and edge supports'),
      'henge':T('逐边 Henge / 骨架','Every-edge Henge / skeleton'),'details':T('选择一个符号项查看共享布局与来源','Select a symbolic term for shared geometry and provenance')}
    title=T('完整反项对应','Complete counterterm correspondence')+f' · L{data["loops"]}E{data["legs"]}'
    v=data['verification'];header=f'{len(data["rows"])} × {len(data["columns"])} = {v["checked_cells"]}'
    return '''<!doctype html><html lang="'''+('zh-CN' if cn else 'en')+'''"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'''+title+'''</title>
<style>body{max-width:1180px;margin:30px auto;padding:0 20px;background:#f2f6f3;color:#163b45;font:16px/1.6 system-ui}section,nav{background:white;border:1px solid #d4dfd9;border-radius:12px;padding:20px;margin:18px 0}select,button{font:inherit;padding:7px;max-width:95%}label{display:inline-block;margin:5px 16px 5px 0}a{color:#176856}.pair{display:flex;flex-wrap:wrap;gap:20px}.pair img{width:330px;max-width:100%}pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:440px;overflow:auto;background:#f5f7f6;padding:12px;font-size:12px}table{border-collapse:collapse;width:100%}th,td{padding:8px;border-bottom:1px solid #d4dfd9;text-align:left}.note{border-left:4px solid #af852f;padding-left:12px}.small{font-size:13px}summary{cursor:pointer;color:#176856}#cols{max-width:520px}</style>
<h1>'''+title+'''</h1><p>'''+header+' · '+T('所有格均已比较；下面可选择任意原图、任意列，包括零项。','Every cell has been compared; select any parent and any column, including zeros.')+'''</p>
<p><a href="index-'''+('en' if cn else 'cn')+'''.html">'''+('English' if cn else '中文')+'''</a> · <a href="../../workspace/index-'''+lang+'''.html">'''+T('统一工作台','Unified workbench')+'''</a> · <a href="table.csv">'''+T('完整矩阵 CSV','Full matrix CSV')+'''</a> · <a href="verification.json">'''+T('逐格验证','Cell verification')+'''</a> · <a href="report.json">'''+T('覆盖清单','Coverage manifest')+'''</a> · <a href="latex/report-'''+lang+'''.tex">LaTeX</a></p>
<nav><label>'''+T('原图','Parent')+''' <select id="rows"></select></label><label>'''+T('列','Column')+''' <select id="cols"></select></label><p id="state"></p></nav>
<section><h2 id="cell-title"></h2><div class="pair" id="images"></div><p id="result"></p><table><thead><tr><th>'''+labels['nodes']+'''</th><th>A</th><th>B</th></tr></thead><tbody id="terms"></tbody></table><label>'''+labels['details']+''' <select id="term-select"></select></label><div class="pair" id="term-images"></div><div id="traces"></div></section>
<section><h2>'''+labels['henge']+'''</h2><select id="henge-select"></select><div class="pair" id="henge-images"></div><details><summary>'''+labels['trace']+'''</summary><pre id="henge-trace"></pre></details><p id="links"></p></section>
<section class="note">'''+T('采用 +gφ³/3! 作用量、−g 裸顶点、固定外腿标签。二点 Taylor 是外动量向量总二阶，分别投影为质量 P⁰、动能 P²；三点为零阶。子图边界可能包含宿主环动量：只投影子图边界依赖，保留宿主传播子和环动量。嵌套算符按包含关系由内向外；重叠子图不能同时相乘。','Action +gφ³/3!, bare vertex −g, fixed external labels. Two-point Taylor has total vector degree two, split into mass P⁰ and kinetic P²; three-point degree zero. Subgraph boundary momenta may include host loop momenta: project source boundary dependence only, retaining host propagators and loop variables. Nested projectors act inside out; overlapping subgraphs cannot occur in one product.')+'''<p>'''+T('完整图集与符号对应不代表积分精度通过。162 组研究总账：严格解析 24、经验验收 48、待解决 90；一点反项有限条件未指定；参考图待外部比对。','Complete graph enumeration and symbolic correspondence are separate from integral accuracy. The research ledger has 24 strict analytic, 48 empirical and 90 unresolved cases out of 162; the finite linear condition is unspecified; reference diagrams await external comparison.')+'''</p></section>
<script>const D='''+json.dumps(data,ensure_ascii=False).replace('</',r'<\/')+''', L='''+json.dumps(labels,ensure_ascii=False)+''';const $=x=>document.getElementById(x);let R=null,active=[],ticket=0;
function option(value,text){const o=document.createElement('option');o.value=value;o.textContent=text;return o}
function picture(src){const im=document.createElement('img');im.src=src;return im}
function trace(title,value){const d=document.createElement('details'),s=document.createElement('summary'),p=document.createElement('pre');s.textContent=title;p.textContent=JSON.stringify(value,null,2);d.append(s,p);return d}
$('rows').replaceChildren(...D.rows.map(r=>option(r.old_id,r.old_id+' · 1/S='+r.weight)));$('cols').replaceChildren(...D.columns.map(c=>option(c.id,c.name)));
async function load(){const row=D.rows.find(r=>r.old_id===$('rows').value),t=++ticket;R=null;try{const response=await fetch(row.data);if(!response.ok)throw Error(response.status);const next=await response.json();if(t!==ticket)return;R=next;$('state').textContent=row.old_id+' · '+row.id+' · '+(R.all_equal?L.equal:L.failed);$('henge-select').replaceChildren(...R.A.henges.map(h=>option(h.minimum_edge,'e'+h.minimum_edge)));$('links').innerHTML='<a href="'+row.data+'">JSON</a> · <a href="'+row.layouts+'">Layout / maps</a>';render();henge()}catch(e){$('state').textContent=L.offline+' '+e}}
function render(){if(!R)return;const cid=$('cols').value,c=D.columns.find(c=>c.id===cid),cell=R.nonzero_cells.find(c=>c.column===cid);$('cell-title').textContent=R.A.old_id+' · '+c.name;$('images').replaceChildren(picture(R.figures.parent.svg));if(c.figure)$('images').append(picture(c.figure));$('result').textContent=cell?(cell.equal?L.equal:L.failed):L.zero;active=R.A.terms.map((t,i)=>({...t,index:i})).filter(t=>t.column===cid);const bt=(R.B?.terms||[]).filter(t=>t.column===cid);$('terms').replaceChildren(...active.map(t=>{const tr=document.createElement('tr');for(const value of [t.nodes.map(n=>n.op+' L'+n.loops+' {'+n.edges.join(',')+'}').join(' → ')||L.bare,t.coefficient,bt.find(x=>x.symbol===t.symbol)?.coefficient??'0']){const td=document.createElement('td');td.textContent=value;tr.append(td)}return tr}));$('term-select').replaceChildren(...active.map((t,i)=>option(i,(i+1)+': '+t.coefficient+' · '+t.symbol.slice(-12))));term()}
function term(){const t=active[Number($('term-select').value)];$('term-images').replaceChildren();$('traces').replaceChildren();if(!t)return;for(const k of ['t'+t.index+'-forest','t'+t.index+'-host'])if(R.figures[k])$('term-images').append(picture(R.figures[k].svg));const bt=(R.B?.terms||[]).filter(x=>x.column===t.column&&x.symbol===t.symbol);$('traces').append(trace(L.a,t),trace(L.b,bt))}
function henge(){if(!R)return;const i=Number($('henge-select').value);$('henge-images').replaceChildren(...['henge-e'+i,'skeleton-e'+i].filter(k=>R.figures[k]).map(k=>picture(R.figures[k].svg)));$('henge-trace').textContent=JSON.stringify(R.A.henges.find(h=>h.minimum_edge===i),null,2)}
$('rows').onchange=load;$('cols').onchange=render;$('term-select').onchange=term;$('henge-select').onchange=henge;const selected=new URLSearchParams(location.search).get('graph');if(D.rows.some(r=>r.old_id===selected))$('rows').value=selected;load();</script></html>'''
