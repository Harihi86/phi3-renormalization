"""One bilingual presentation layer over the exact closure dataset."""
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import html
import json
import math
from zipfile import ZipFile, ZIP_DEFLATED
from .graphs import Graph
from .reporting import svg as legacy_svg


from .ct_layout import layout, geometry, svg as shared_svg, tikz as shared_tikz, export_family


def draw(g,marked=(),selected=None,bridges=(),positions=None):
    return shared_svg(g,geometry(g,positions),marked,selected,bridges)


def tikz(g,marked=(),positions=None):
    return shared_tikz(g,geometry(g,positions),marked)


def coefficient_text(cell,names):
    if cell['zero']:return '0'
    parts=[]
    for symbol,value in cell['A'].items():
        if cell['column']!='C3L2':label=''
        elif symbol.startswith('T0[I:'):label=' T₀I'
        else:
            col=symbol.split(':')[2];label=' T₀'+names.get(col,'K')
        q=Fraction(value);parts.append(('− ' if q<0 else '+ ')+str(abs(q))+label)
    return ' '.join(parts).removeprefix('+ ')


def contraction_positions(row,term):
    g=Graph.from_dict(row['graph']);p=layout(g)
    sub=next(s for s in row['recursive_ast']['proper_parts'] if s['edge_ids']==term['embedded_edges'])
    vertices=sub['vertices'];raw={0:tuple(sum(p[v][k] for v in vertices)/len(vertices) for k in (0,1))}
    j=1
    for v in range(g.n):
        if v not in vertices:raw[j]=p[v];j+=1
    return {term['contraction_mapping']['vertices'][v]:xy for v,xy in raw.items()}


def write_reports(data,a,b,out):
    out=Path(out);fig=out/'figures';fig.mkdir(exist_ok=True)
    latex=out/'latex';latex.mkdir(exist_ok=True)
    columns=data['columns'];names={c['id']:c['name'] for c in columns};names.update(bare='I',C3L2='C3L2')
    cells={(c['row'],c['column']):c for c in data['verification']['cells']}
    for row in a['rows']:
        g=Graph.from_dict(row['graph']);(fig/(row['old_id']+'.svg')).write_text(draw(g),encoding='utf-8')
        (latex/(row['old_id']+'.tex')).write_text(tikz(g),encoding='utf-8')
        for h in row['henges']:
            marked=[e for block in h['blocks'] if block['loops'] for e in block['edge_ids']]
            (fig/f'{row["old_id"]}-henge-e{h["minimum_edge"]}.svg').write_text(draw(g,marked,h['minimum_edge'],h['bridges']),encoding='utf-8')
            vm=h['original_vertex_to_skeleton'];sg=Graph(len(h['blocks']),tuple(tuple(e) for e in h['skeleton_edges']),tuple((l,vm[v]) for l,v in g.external),tuple('H' if block['loops'] else 'phi3' for block in h['blocks']))
            pos={i:tuple(sum(layout(g)[v][k] for v in block['vertices'])/len(block['vertices']) for k in (0,1)) for i,block in enumerate(h['blocks'])}
            (fig/f'{row["old_id"]}-skeleton-e{h["minimum_edge"]}.svg').write_text(draw(sg,positions=pos),encoding='utf-8')
        for term in row['terms']:
            if term['kind']=='insertion':
                (fig/f'{row["old_id"]}-{names[term["column"]]}-subgraph.svg').write_text(draw(g,term['embedded_edges']),encoding='utf-8')
                (fig/f'{row["old_id"]}-{names[term["column"]]}-host.svg').write_text(draw(Graph.from_dict(term['host']),positions=contraction_positions(row,term)),encoding='utf-8')
    for c in columns:
        g=Graph.from_dict(c['host']);(fig/(c['name']+'.svg')).write_text(draw(g),encoding='utf-8');(latex/(c['name']+'.tex')).write_text(tikz(g),encoding='utf-8')
    local=Graph(1,(),((1,0),(2,0),(3,0)),('C3L2',))
    (fig/'C3L2.svg').write_text(draw(local),encoding='utf-8')
    (latex/'C3L2.tex').write_text(tikz(local),encoding='utf-8')
    export_family(a['rows'],columns,out)
    for lang in ('cn','en'):
        _html(data,a,b,out,lang,cells,names)
        _latex(data,out,lang,cells,names)
    for folder,name in [(latex,'latex-source.zip'),(fig,'figures-svg.zip')]:
        with ZipFile(out/name,'w',ZIP_DEFLATED) as archive:
            for path in sorted(folder.iterdir()):
                if path.suffix in ('.tex','.svg'):archive.write(path,path.name)


def _html(data,a,b,out,lang,cells,names):
    cn=lang=='cn';T=lambda zh,en:zh if cn else en;esc=html.escape
    cols=['bare']+[c['id'] for c in data['columns']]+['C3L2']
    b_rows={r['id']:r for r in b['rows']};b_cols={c['id']:c for c in b['columns']}
    column_info={c['id']:c for c in data['columns']}
    title=T('六维 φ³ · 两圈三点反项对应表','Six-dimensional φ³ · Two-loop three-point counterterms')
    nav=' · '.join([f'<a href="{file}" download>{name}</a>' for file,name in [
        ('correspondence.json',T('精确数据','Exact data')),('table.csv',T('对应表 CSV','Table CSV')),
        ('latex-source.zip',T('可编辑 LaTeX','Editable LaTeX')),('figures-svg.zip',T('全部 SVG','All SVGs')),('path-a.json',T('路径 A 记录','Path A trace')),
        ('path-b.json',T('路径 B 记录','Path B trace')),('verification.json',T('逐单元格核验','Cell verification'))]])
    headers=''.join(f'<th>{esc(names[c])}</th>' for c in cols)
    table=[];details=[];hengesections=[]
    for ri,row in enumerate(a['rows']):
        gid=row['id'];old=row['old_id'];g=Graph.from_dict(row['graph']);td=[]
        for ci,col in enumerate(cols):
            cell=cells[gid,col];key=f'cell-{ri}-{ci}';txt=coefficient_text(cell,names)
            thumb=''
            if not cell['zero'] and col in column_info:thumb=f'<img src="figures/{names[col]}.svg" alt="{names[col]}">'
            if col=='C3L2':thumb='<img src="figures/C3L2.svg" alt="C3L2">'
            td.append(f'<td><button class="cell {"zero" if cell["zero"] else "nonzero"}" data-open="{key}" aria-controls="{key}">{thumb}<span>{txt}</span></button></td>')
            ats=[t for t in row['terms'] if t['column']==col];bts=[t for t in b_rows[gid]['terms'] if t['column']==col]
            intro=T('两条路径均为结构零：此固定标签原图与该类型宿主之间没有可用嵌入或拼接。','Both paths give a structural zero: no embedding or gluing connects this labelled parent to the typed host.') if cell['zero'] else T('两条独立路径的图结构、符号与有理系数完全相等。','Both independent paths agree in graph structure, sign and rational coefficient.')
            content=f'<p class="ok">{intro}</p><p class="formula">{esc(txt)}</p><div class="pair"><figure>{draw(g,ats[0].get("embedded_edges",[]) if ats else [])}<figcaption>{old}</figcaption></figure>'
            if col in column_info:
                content+=f'<figure>{draw(Graph.from_dict(column_info[col]["host"]),positions=contraction_positions(row,ats[0]) if ats else None)}<figcaption>{names[col]}</figcaption></figure>'
            content+='</div>'
            content+=f'<p><a href="figures/{old}.svg" download>{T("下载原图 SVG","Download parent SVG")}</a>'
            if col in column_info:content+=f' · <a href="figures/{names[col]}.svg" download>{T("下载反项图 SVG","Download countergraph SVG")}</a>'
            content+='</p>'
            if ats and ats[0]['kind']=='insertion':
                ta=ats[0];bc=b_cols[col];w=bc['host_wick'];s=bc['source_wick'];factor=bts[0]['coefficient_factors']
                content+=f'<h3>{T("系数来源","Coefficient provenance")}</h3><p>{T("路径 A：每处子图贡献","Path A: each embedded subgraph contributes")} Sγ/S<sub>G</sub> = {ta["parent_embedding_factor"]}; {T("反项自身权重","countervertex source weight")} 1/Sγ = {bc["source_weight"]}.</p>'
                content+=f'<p>{T("路径 B 宿主 Wick 配对","Path B host Wick pairings")}: {w["wick_multiplicity"]} / ({w["species_factorial"]} × {w["vertex_factorial"]}) = {bc["host_weight"]}. '
                content+=f'{T("子图权重","Source weight")} {bc["source_weight"]}; {T("端口双射数","port bijections")} {bc["port_bijections"]}; {T("每次平均因子","factor per bijection")} {factor["port_average"]}. {T("带符号列和","Signed sum")}: {txt}.</p>'
                t=ta['taylor'];content+=f'<h3>{T("Taylor 变量与减除","Taylor variables and subtraction")}</h3><p>{T("总次数","Total degree")} {t["degree"]}; {T("独立边界动量","independent boundary momenta")}: {", ".join(t["independent_boundary_symbols"])}; {T("保留内部环动量","retained inner loop")}: {", ".join(t["inner_loop_basis"])}. '
                content+=T('外层环动量在子图中作为边界输入参与 Taylor 展开；宿主传播子保持不动。','An outer-loop momentum is a boundary input for the subgraph Taylor expansion; host propagators remain unchanged.')+'</p>'
                if ta['operator_type']=='C2z':content+='<p class="formula">K<sub>z</sub> = ∫ [−D(q)³ + (4/6)q²D(q)⁴]; &nbsp; C<sub>2z</sub>(p) = −g² K<sub>z</sub> p²/2.</p>'
                content+=f'<p>{T("展开后的裸图、端口、半边与环变量的精确映射见下方路径 B 记录。","Expanded parents, ports, half-edges and loop-variable maps are recorded in path B below.")}</p>'
            if col=='C3L2':content+=f'<p>{T("此格是多个局部符号项，不能把系数相加当作积分值。路径 B 独立用三点函数零外动量条件固定新两圈反项；保留两条环动量。","This cell contains distinct local symbolic terms, not a numerical integral. Path B fixes the new two-loop countervertex by the zero-external-momentum condition while retaining both loop momenta.")}</p>'
            if col=='bare':
                w=next((t.get('wick_record') for t in bts if 'wick_record' in t),None)
                if w:content+=f'<p>{T("生成泛函配对系数","Generating-functional Wick coefficient")}: −{w["wick_multiplicity"]}/({w["species_factorial"]} × {w["vertex_factorial"]}) = {row["bare_amplitude_coefficient"]}; S = {row["symmetry"]}.</p>'
            for label,trace in [('A',ats),('B',bts)]:content+=f'<details><summary>{T("路径","Path")} {label} · {T("完整机器记录","complete machine trace")}</summary><pre>{esc(json.dumps(trace,ensure_ascii=False,indent=2))}</pre></details>'
            details.append(f'<section class="detail" id="{key}" hidden><h2>{old} · {names[col]}</h2>{content}</section>')
        table.append(f'<tr><th><a href="#henge-{ri}">{old}</a><img class="parent" src="figures/{old}.svg" alt="{old}"><small>1/S = {row["positive_weight"]}</small></th>{"".join(td)}</tr>')
        hs=[]
        for h in row['henges']:
            eid=h['minimum_edge'];blocks='; '.join(f'{block["edge_ids"]} (ω={block["omega"]})' for block in h['blocks'] if block['loops']) or '∅'
            hs.append(f'<details><summary>e{eid} · {T("含圈块","loop blocks")}: {blocks}</summary><div class="pair"><img src="figures/{old}-henge-e{eid}.svg" alt="Henge"><img src="figures/{old}-skeleton-e{eid}.svg" alt="{T("骨架","Skeleton")}"></div><pre>{esc(json.dumps(h,ensure_ascii=False,indent=2))}</pre></details>')
        hengesections.append(f'<details class="card" id="henge-{ri}"><summary>{old} · Henge / BPH</summary><p><code>{gid}</code></p>{"".join(hs)}<details><summary>{T("实际递归减除程序","Actual recursive subtraction program")}</summary><pre>{esc(json.dumps(row["recursive_ast"],ensure_ascii=False,indent=2))}</pre></details></details>')
    countrows=''.join(f'<tr><td>{r["loops"]}</td><td>{r["legs"]}</td><td>{r["functional_classes"]}</td><td>{r["normalized_classes"]}</td><td>{r["connected_amputated_classes"]}</td><td>{r["one_pi_classes"]}</td><td>{r["unfixed_orbits"]}</td><td>{r["weight_sum"]}</td><td><a href="{r["audit_file"]}">{T("候选与理由","Candidates and reasons")}</a></td></tr>' for r in data['count_audit'])
    conventions=T('作用量含 +gφ³/3!，裸顶点为 −g，所有表项提出 g⁵。外腿标签固定；保留 tadpole。二点减除为向量外动量的总二阶（质量和动能分开），三点为零阶。局部积分系数保持符号形式。',
        'The action contains +gφ³/3!, the bare vertex is −g, and g⁵ is factored out of every cell. External labels are fixed; tadpoles are retained. Two-point Taylor subtraction has total vector-momentum degree two, split into mass and kinetic terms; three-point subtraction has degree zero. Local integrals remain symbolic.')
    boundary=T('未取得参考图：待外部逐图比对。历史“10”包含两圈一、二、三点的合计或另一项插入位置计数，不能当作两圈三点图数。162 组研究总账中严格解析 24 组、经验验收 48 组、待解决 90 组；本表生成不增加积分验收。',
        'Reference diagrams are unavailable: external graph-by-graph comparison is pending. The historical “10” counts either all two-loop one-, two- and three-point graphs or insertion positions in a different test. The research ledger records 24 strict analytic, 48 empirical and 90 unresolved cases out of 162; this exact build adds no integral acceptance.')
    title_lang='zh-CN' if cn else 'en'
    page=f'''<!doctype html><html lang="{title_lang}"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title>
<style>*{{box-sizing:border-box}}body{{margin:0;color:#163b45;background:#f2f5f3;font:15px/1.6 system-ui,sans-serif}}header,main{{max-width:1580px;margin:auto;padding:26px}}header{{padding-bottom:12px}}h1{{font-size:31px;margin:8px 0}}h2{{font-size:23px}}h3{{font-size:17px}}a{{color:#176856}}.eyebrow{{letter-spacing:.15em;font-size:12px;color:#63766e}}.stats{{display:flex;gap:16px;flex-wrap:wrap;margin:20px 0}}.stats div{{background:white;border:1px solid #d4dfd9;border-radius:10px;padding:16px 24px;min-width:160px}}.stats strong{{font-size:29px;display:block}}.stats span,small{{font-size:12px;color:#63766e}}.card,.detail{{background:white;border:1px solid #d4dfd9;border-radius:10px;padding:20px;margin:18px 0}}.table-wrap{{overflow:auto;background:white;border:1px solid #d4dfd9;border-radius:10px}}table{{width:100%;border-collapse:collapse;font-size:13px}}th,td{{border-bottom:1px solid #d4dfd9;border-right:1px solid #e6eeea;padding:6px;vertical-align:top;text-align:center}}thead th{{background:#e7efea;position:sticky;top:0;min-width:88px}}th:first-child{{min-width:130px}}.cell{{width:100%;min-height:92px;border:0;background:white;border-radius:6px;color:#163b45;cursor:pointer;font:14px/1.5 system-ui}}.cell img{{display:block;width:75px;height:59px;margin:auto}}.cell span{{display:block}}.cell:hover,.cell[aria-pressed=true]{{background:#e8f4ee;outline:2px solid #176856}}.zero{{color:#8c9b96}}.parent{{display:block;width:120px;height:95px}}.formula{{font:18px/1.8 Georgia,serif}}.pair{{display:flex;gap:22px;flex-wrap:wrap}}figure{{margin:0;max-width:350px}}figure svg,.pair>img{{width:330px;max-width:100%}}figcaption{{text-align:center}}pre{{max-height:360px;overflow:auto;text-align:left;background:#f6f8f7;border:1px solid #e1e9e5;padding:15px;font-size:12px}}summary{{cursor:pointer;color:#176856;font-weight:600;padding:8px 0}}.ok{{color:#176856}}.note{{border-left:4px solid #b58a32;padding-left:15px}}.legend{{font-size:13px;color:#63766e}}.detail:focus{{outline:2px solid #176856}}[hidden]{{display:none!important}}@media(max-width:600px){{header,main{{padding:14px}}h1{{font-size:25px}}}}@media print{{button{{color:black!important}}.table-wrap{{overflow:visible}}}}</style>
<header><p class="eyebrow">PHI³ LAB 4.3 · {T("精确组合闭环","EXACT COMBINATORIAL CLOSURE")}</p><h1>{title}</h1><p>{conventions}</p><nav>{nav} · <a href="index-{'en' if cn else 'cn'}.html">{'English' if cn else '中文'}</a></nav>
<div class="stats"><div><strong>7</strong><span>{T("固定标签裸图","Fixed-label bare graphs")}</span></div><div><strong>9</strong><span>{T("低阶反项插入列","Lower-order insertion columns")}</span></div><div><strong>77 / 77</strong><span>{T("逐格精确相等","Cells exactly equal")}</span></div><div><strong>32</strong><span>{T("非零符号项","Nonzero symbolic monomials")}</span></div></div></header>
<main><p>{T("点击任一系数或零项，查看图、阶乘、Wick 配对、端口拼接、Taylor 变量及两条路径记录。Mij 为边 ij 上的质量反项，Zij 为动能反项，Gi 为带外腿 i 的耦合反项。", "Select a coefficient or zero to inspect graphs, factorials, Wick pairings, port gluing, Taylor variables and both traces. Mij denotes a mass insertion on edge ij, Zij a kinetic insertion, and Gi a coupling insertion carrying external leg i.")}</p>
<div class="table-wrap"><table id="correspondence"><thead><tr><th>{T("原图","Parent")}</th>{headers}</tr></thead><tbody>{''.join(table)}</tbody></table></div>
<p class="legend">{T("图例：橙色为子图，红色虚线为所选删边，紫色为桥，金色方框为局部反项；无节点的交叉不是顶点。", "Legend: orange subgraph; red dashed selected/deleted edge; purple bridge; gold square local countervertex. Crossings without a node are not vertices.")}</p>
<div id="details" aria-live="polite">{''.join(details)}</div>
<section class="card"><h2>{T("从生成泛函到截肢 1PI","From the generating functional to amputated 1PI")}</h2><p class="formula">Z[J] → Z[J]/Z[0] → W[J] = log(Z[J]/Z[0]) → δ³W/δJ³|₀ → {T("截肢","amputation")} → 1PI</p>
<p>{T("实际枚举包含源与源直接配对、真空因子、非连通图和带内部桥的图。先分离真空分量，再连通化、截肢，最后筛去内部桥；令 J=0 并不删除 tadpole。", "The enumeration includes direct source-source pairings, vacuum factors, disconnected graphs and internal bridges. Vacuum normalization, connectedness, amputation and the internal-bridge filter are distinct operations; setting J=0 does not delete tadpoles.")}</p>
<div class="table-wrap"><table><thead><tr><th>L</th><th>E</th><th>Z</th><th>Z/Z₀</th><th>{T("连通截肢","Connected amputated")}</th><th>1PI</th><th>{T("无标签轨道","Unfixed orbits")}</th><th>Σ1/S</th><th>{T("追溯","Trace")}</th></tr></thead><tbody>{countrows}</tbody></table></div></section>
<section><h2>{T("逐边 Henge 与递归减除","Every-edge Henge and recursive subtraction")}</h2>{''.join(hengesections)}</section>
<section class="card"><h2>{T("验证边界与继续工作","Validation boundary and continuation")}</h2><p class="note">{boundary}</p><p>{T("一点反项有限条件未指定，目标三点 1PI 中的一点反项列因内部桥而排除。质量项、动能项和耦合项按符号积分基底精确比较，不能由数值接近替代。", "The finite linear-counterterm condition remains unspecified; one-point counterterm hosts are excluded from this three-point 1PI family by an internal bridge. Mass, kinetic and coupling terms are compared exactly in a symbolic integral basis.")}</p><p><a href="https://github.com/Harihi86/phi3-renormalization/blob/main/{'README_ZH.md' if cn else 'README.md'}">{T("使用与理论说明","Usage and theory notes")}</a> · <a href="https://arxiv.org/abs/hep-th/9612113">Kennedy 1996</a> · <a href="https://arxiv.org/abs/0712.1016">Kennedy–Binoth–Rippon 2007</a></p></section></main>
<script>document.addEventListener('click',ev=>{{const b=ev.target.closest('[data-open]');if(!b)return;document.querySelectorAll('.detail').forEach(x=>x.hidden=x.id!==b.dataset.open);document.querySelectorAll('[data-open]').forEach(x=>x.setAttribute('aria-pressed',x===b?'true':'false'));const d=document.getElementById(b.dataset.open);d.tabIndex=-1;d.scrollIntoView({{behavior:'smooth',block:'start'}});d.focus({{preventScroll:true}});}});</script></html>'''
    (out/f'index-{lang}.html').write_text(page,encoding='utf-8')


def _latex(data,out,lang,cells,names):
    cn=lang=='cn';T=lambda zh,en:zh if cn else en
    header=[r'\documentclass[UTF8]{ctexart}' if cn else r'\documentclass{article}',
            r'\usepackage{amsmath,amssymb,booktabs,longtable,tikz,graphicx}',
            r'\usepackage[margin=16mm]{geometry}',r'\begin{document}',
            r'\section*{'+T(r'六维 $\phi^3$ 两圈三点反项对应表',r'Six-dimensional $\phi^3$: two-loop three-point counterterms')+'}',
            T('固定外腿标签；截肢连通 1PI；所有表项提出 $g^5$。局部积分保留为符号。','Fixed external labels; amputated connected 1PI; $g^5$ factored out. Local integrals remain symbolic.'),
            r'\[S_{int}=g\phi^3/3!,\quad v_3=-g,\quad C_\gamma=-t_\gamma\bar R_\gamma,\quad R_G=(1-t_G)\bar R_G.\]',
            r'\[t_2: \deg(p)\le2,\qquad t_3: \deg(p_1,p_2)=0.\]',
            r'\[K_m=\int D(q)^2,\quad K_z=\int[-D(q)^3+\tfrac46q^2D(q)^4],\quad K_3=\int D(q)^3.\]',
            r'\[v_{2m}=-g^2K_m/2,\quad v_{2z}(p)=-g^2K_zp^2/2,\quad v_{3,1}=+g^3K_3.\]']
    cols=['bare']+[c['id'] for c in data['columns']]+['C3L2']
    header += [r'\begin{center}\resizebox{\textwidth}{!}{\begin{tabular}{l'+'c'*len(cols)+'}',r'\toprule',
               'Graph & '+' & '.join(names[c] for c in cols)+r' \\\midrule']
    for row in data['rows']:
        vals=[]
        for col in cols:
            cell=cells[row['id'],col]
            txt=coefficient_text(cell,names).replace('−','-').replace('T₀',r'T_0 ').replace(' ','')
            if col=='C3L2':txt=txt.replace('T_0',r'T_0\,')
            vals.append('$'+txt+'$')
        header.append(row['old_id']+' & '+' & '.join(vals)+r' \\')
    header += [r'\bottomrule\end{tabular}}\end{center}',T('每一格均由两条独立算法进行精确核验。零项为结构排除，不是近似零。','Every cell is checked by two independent algorithms. Zeros are structural, not approximate.'),
        T('参考图尚缺：待外部逐图比对。162 组研究总账：严格解析 24、经验验收 48、待解决 90。','External reference diagrams remain unavailable. The 162-case research ledger records 24 strict analytic, 48 empirical and 90 unresolved cases.')]
    for row in data['rows']:
        header += [r'\subsection*{'+row['old_id']+'}',r'\input{'+row['old_id']+'.tex}',
                   T('正权重','Positive weight')+': $1/S='+row['positive_weight']+'$.',
                   r'\begin{itemize}']
        for col in cols:
            cell=cells[row['id'],col]
            if cell['zero']:continue
            txt=coefficient_text(cell,names).replace('−','-').replace('T₀',r'T_0\,')
            header.append(r'\item '+names[col]+': $'+txt+'$')
        header.append(r'\end{itemize}')
    header.append(r'\section*{'+T('反项宿主图','Counterterm host graphs')+'}')
    for c in data['columns']:header += [r'\subsection*{'+c['name']+'}',r'\input{'+c['name']+'.tex}']
    header += [r'\end{document}']
    (out/'latex'/f'main-{lang}.tex').write_text('\n'.join(header),encoding='utf-8')
