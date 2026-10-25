"""Portable offline HTML/SVG and LaTeX exports; no CDN, server, or JS libraries."""
import html,json,math,os
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import quote
from .graphs import Graph
from .routing import momentum_text
from .provenance import external_evidence
from .two_loop_reporting import two_loop_panel,TWO_LOOP_SCRIPT


def svg(g):
    """Draw topology, not embedding physics: unmarked crossings are not vertices."""
    w,h=360,260;cx,cy=180,128
    if g.n==1:positions={0:(cx,cy)}
    else:
        positions={i:(cx+86*math.cos(-math.pi/2+2*math.pi*i/g.n),
                      cy+77*math.sin(-math.pi/2+2*math.pi*i/g.n)) for i in range(g.n)}
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" role="img" aria-label="Feynman multigraph">',
           '<style>.edge,.leg{stroke:#234b57;stroke-width:2;fill:none}.edge.highlight{stroke:#dc6a22;stroke-width:5}.vertex{fill:#234b57}.ct{fill:#f2ebd9;stroke:#926b1e}.ct-text{font-size:9px;text-anchor:middle;fill:#55401b}text{text-anchor:middle;font:13px system-ui;fill:#34525b}.edge-label{font-size:11px;paint-order:stroke;stroke:white;stroke-width:4;stroke-linejoin:round}.vertex-label{font-size:10px;fill:#80918c}</style>']
    grouped=defaultdict(list)
    for i,e in enumerate(g.edges):grouped[e].append(i)
    for (a,b),ids in grouped.items():
        x1,y1=positions[a];x2,y2=positions[b]
        for j,i in enumerate(ids):
            if a==b:
                radius=28+10*j
                path=f'M {x1-7} {y1} C {x1-radius*2} {y1-radius*2}, {x1+radius*2} {y1-radius*2}, {x1+7} {y1}'
                tx,ty=x1,y1-radius*1.5-5
            else:
                dx,dy=x2-x1,y2-y1;length=math.hypot(dx,dy)
                offset=(j-(len(ids)-1)/2)*34
                mx,my=(x1+x2)/2-offset*dy/length,(y1+y2)/2+offset*dx/length
                path=f'M {x1} {y1} Q {mx} {my} {x2} {y2}'
                tx,ty=(x1+x2)/2-offset*dy/(2*length),(y1+y2)/2+offset*dx/(2*length)
                tx+=8*dy/length;ty-=8*dx/length
            parts.append(f'<path class="edge" data-eid="{i}" d="{path}"/><text class="edge-label" x="{tx}" y="{ty}">e{i}</text>')
    attachments=defaultdict(list)
    for label,v in g.external:attachments[v].append(label)
    for v,labels in attachments.items():
        x,y=positions[v]
        angle=math.atan2(y-cy,x-cx) if g.n>1 else math.pi
        for j,label in enumerate(labels):
            a=angle+(j-(len(labels)-1)/2)*.4
            xx,yy=x+48*math.cos(a),y+48*math.sin(a)
            parts.append(f'<path class="leg" d="M{x} {y}L{xx} {yy}"/><text x="{xx+6*math.cos(a)}" y="{yy+6*math.sin(a)+4}">p{label}</text>')
    for v,(x,y) in positions.items():
        kind=g.vertex_kinds[v]
        if kind.startswith('C'):
            parts.append(f'<rect class="ct" x="{x-10}" y="{y-10}" width="20" height="20"/><text class="ct-text" x="{x}" y="{y+4}">{kind}</text>')
        else:parts.append(f'<circle class="vertex" cx="{x}" cy="{y}" r="4"/><text class="vertex-label" x="{x+9}" y="{y+15}">v{v}</text>')
    return ''.join(parts)+'</svg>'


def operator_html(term):
    result=[]
    for op in reversed(term['operators_inner_first']):
        label=op['subgraph'].replace('gamma','γ')
        result.append(f't<sub>{label}</sub><sup>{op["degree"]}</sup>')
    return ('− ' if term['sign']<0 else '+ ')+''.join(result)+'I<sub>G</sub>'


def monomial_html(powers):
    return ' '.join(f'x<sub>{i}</sub>'+ (f'<sup>{power}</sup>' if power>1 else '')
                    for i,power in enumerate(powers) if power) or '1'


def parametric_html(record):
    p=record['symanzik']
    u=' + '.join(monomial_html(t['powers']) for t in p['U'])
    f0=[]
    for t in p['F_kinematic']:
        q=html.escape(momentum_text(t['momentum'],p['external_basis']))
        f0.append(f'({q})<sup>2</sup> '+monomial_html(t['powers']))
    ident=record['id']
    return f'''<h4>参数积分：精确 Symanzik 多项式</h4>
<p class="formula">U = {u}</p><p class="formula">F = m² U (Σ<sub>e</sub> x<sub>e</sub>) + {" + ".join(f0) or '0'}</p>
<p class="meta">生成树 {p['spanning_trees_count']} · 生成双森林 {p['spanning_two_forests_count']} · U 次数 {record['loops']} · F 次数 {record['loops']+1}</p>
<p class="good">✓ 矩阵树有理点核对 · ✓ Euclidean 正值采样 · ✓ 齐次性</p>
<p><a href="inputs/{ident}.json" download>下载本图输入</a> · <a href="pysecdec/{ident}/integral.json" download>pySecDec 约定</a> · <a href="pysecdec/{ident}/generate.py" download>生成脚本</a> · <a href="pysecdec/{ident}/integrate.py" download>积分脚本</a></p>
<p class="muted">所有外动量流入，最后一条为其余的负和。F 中平方为 Euclidean 标量积。导出是裸 d=D−2ε 积分；脚本的生成、C++ 编译和积分需分别执行。<a href="pysecdec/README.txt">后端使用说明</a>。</p>'''


def overlap_panel(data):
    result=data.get('two_loop_overlap')
    if not result:return ''
    bound=result['analytic_bounds']
    rows=''.join(f'<tr><td>{v["order"]}</td><td>{v["value"]:.16g}</td></tr>' for v in result['quadrature'])
    payload=json.dumps({**result,'s':1,'m2':1,'order':16,'weighted_value':result['value']/2},ensure_ascii=False).replace('<','\\u003c')
    return f'''<section class="panel overlap" id="overlap-workbench"><p class="eyebrow">2.0 · 完整两圈重叠图基准</p><h2>BPH 递归：减除、有限积分与解析界</h2>
<p><a href="#L2E2G1">查看本图与嵌入子图</a> · <a href="../docs/overlap-analysis.html">完整解析推导</a> · <a href="../docs/overlap-analysis.html#momentum-proof">动量空间收敛界</a> · <a href="../analysis-project/main.pdf">完整 PDF</a></p>
<p class="formula">I = D(k)D(p−k)D(q)D(p−q)D(k−q)<br>R I = (1−t<sub>p</sub><sup>2</sup>)[I−D(k)³D(q)D(p−q)−D(q)³D(k)D(p−k)]</p>
<p>两个三角子图在中央边重叠，因此递归中只有两项子图反项；它们不相乘。整体二阶 Taylor 减除随后作用，合计六项。</p>
<div class="result-value">{result['value']:.16g}</div><p class="meta">已生成基准：p²=m²=1；测度 d⁶k/π³ d⁶q/π³；不含耦合和图权重。图权重 1/2 后为 {result['value']/2:.16g}。</p>
<details><summary>积分阶数与验证范围</summary><table><tr><th>每维阶数</th><th>完整有限积分</th></tr>{rows}</table><p>最后阶数差：{result['last_refinement_difference']:.2g}。另用不同坐标与五维有理 Taylor 余项核独立核对；原始证据见下方核查表。</p><p class="muted">阶数差和路线一致性是数值诊断，不是严格积分误差界。</p></details>
<h3>在本机重新计算</h3><form id="overlap-form"><div class="calculation-inputs"><label>外动量平方 p²<input id="overlap-s" name="s" type="number" min="0" step="any" value="1" required></label><label>共同质量平方 m²<input id="overlap-m2" name="m2" type="number" min="0" step="any" value="1" required></label><label>每维积分阶数<select id="overlap-order"><option>12</option><option selected>16</option><option>24</option><option>32</option></select></label></div><button type="submit" id="overlap-calculate">计算有限积分</button> <button type="button" id="overlap-download">下载当前数值 JSON</button></form>
<p id="overlap-live-status" class="muted" role="status" aria-live="polite">双击软件目录中的启动文件可开启本地计算；直接打开 HTML 仍可查看全部已生成结果。</p><pre id="overlap-live-result" hidden></pre><script type="application/json" id="overlap-initial">{payload}</script>
<h3>解析上界与紫外收敛</h3><p class="formula">∫ |合并后的参数核| ≤ min(C s²/m², s[7L/4+L²/4])<br>C = 11 log 2 / 24 − 29 / 384 ≈ {bound['quadratic_constant']:.10g}；L = log(1+s/(4m²))，s=p²。</p>
<p>小外动量界为 O(p⁴/m²)，大外动量界为 O(p² log²(p²/m²))。对 |p|≤P、m²≥m₀²，可直接代入 P² 和 m₀² 得到统一界。原始动量积分另按单硬圈与可比双圈区域证明绝对收敛，包含 k≈q。</p>
<p class="muted">上界控制此图的振幅与收敛，不能作为数值答案的误差条。结论限于本图、有质量 Euclidean 六维、零动量 Taylor 条件。</p></section>'''


def write_latex(data,path):
    lines=[r'\documentclass{article}',r'\usepackage{amsmath,amssymb,booktabs,longtable}',
           r'\usepackage[margin=20mm]{geometry}',r'\begin{document}',
           r'\section*{Massive $\phi^3$: BPH recursion benchmarks}',
           r'Fixed labelled external legs; amputated 1PI graphs. Formal operators, not evaluated integrals.',
           r'\begin{longtable}{lrrrrr}\toprule Graph & $L$ & $E$ & $I$ & $S$ & $\omega$\\\midrule']
    for r in data['graphs']:
        lines.append(f"{r['id']} & {r['loops']} & {r['external_legs']} & {len(r['edges'])} & {r['symmetry']} & {r['omega']} "+r'\\')
    lines.append(r'\bottomrule\end{longtable}')
    for r in data['graphs']:
        lines.append(r'\subsection*{'+r['id']+'}')
        expr=[]
        for term in r['forest_terms']:
            ops=''.join(r't_{'+op['subgraph'].replace('gamma',r'\gamma_')+'}^{'+str(op['degree'])+'}'
                        for op in reversed(term['operators_inner_first']))
            expr.append(('-' if term['sign']<0 else '+')+ops+r'I_G')
        lines.append(r'\[R_G I_G='+''.join(expr).lstrip('+')+r'.\]')
        lines.append('The rightmost Taylor operator acts first. All integrals are understood with a common UV regulator.\\')
    lines.extend([r'\end{document}'])
    Path(path).write_text('\n'.join(lines),encoding='utf-8')


def write_html(data,path):
    escape=html.escape
    cards=[]
    for r in data['graphs']:
        g=Graph.from_dict(r);subs=r['proper_divergent_subgraphs'];hs=r['henges']
        subrows=''.join(f'<tr><td><button data-highlight="{",".join(map(str,s["edge_ids"]))}">γ{j+1}</button></td><td>{s["edge_ids"]}</td><td>{s["external_legs"]}</td><td>{s["omega"]}</td><td>{s["symmetry"]}</td></tr>' for j,s in enumerate(subs))
        subtable='<p>无真 UV 发散子图。</p>' if not subs else '<table><tr><th>高亮</th><th>内部边</th><th>E</th><th>ω</th><th>Sγ</th></tr>'+subrows+'</table>'
        hrows=[]
        for h in hs:
            desc=[]
            for b in h['blocks']:
                if b['loops']:
                    desc.append(f"{b['edge_ids']} (L={b['loops']}, ω={b['omega']}, "+('需减除' if b['needs_uv_subtraction'] else '收敛块')+')')
            hrows.append(f'<tr><td>e{h["minimum_edge"]}</td><td>'+('; '.join(desc) or '仅裸顶点')+f'</td><td>{h["skeleton_loops"]}</td></tr>')
        routes=''.join(f'<li>e{i}: <code>{momentum_text(c,r["routing"]["basis"])}</code></li>' for i,c in r['routing']['edge_momenta'].items())
        terms=' '.join(operator_html(t) for t in r['forest_terms']).lstrip('+ ')
        overlapping=any(not any(i in f and j in f for f in r['proper_forests'])
                        for i in range(len(subs)) for j in range(i+1,len(subs)))
        uv='overlap' if overlapping else 'subdivergent' if subs else 'primitive'
        cards.append(f'''<article class="card" id="{r['id']}" data-uv="{uv}" data-loops="{r['loops']}" data-legs="{r['external_legs']}" data-id="{r['id']}">
<h3>{r['id']} <span class="badge">1/S = {r['weight']}</span></h3>{svg(g)}
<p class="meta">V={r['vertices']} · I={len(r['edges'])} · L={r['loops']} · E={r['external_legs']} · ω={r['omega']} · <a href="figures/{r['id']}.svg" download>SVG 图</a></p>
<details><summary>查看减除、参数多项式与后端输入</summary><p><a href="#{r['id']}">本图直达链接</a> · <button data-clear-highlight>清除高亮</button></p><h4>按边编号嵌入的真发散子图</h4>{subtable}
<h4>BPH 递归的显式展开（森林记法对照）</h4><p class="formula">R I<sub>G</sub> = {terms}</p><p class="muted">由 Bogoliubov 递归展开生成，另与森林枚举对照；只核查本图的组合项。右边算子先作用。此处是形式算子，没有一般积分求值或一般定理证明。</p>
<h4>Henge：包括收敛的含圈块</h4><table><tr><th>最小线</th><th>删线后的含圈 1PI 块</th><th>骨架 L</th></tr>{''.join(hrows)}</table>
<h4>自动动量路由</h4><ul>{routes}</ul><p class="muted">边朝向为小编号顶点到大编号顶点；所有外动量流入，最后一个外动量由守恒确定；一点函数的 p1=0。无点交叉不表示顶点。</p>
<p class="good">✓ 逐顶点动量守恒 · ✓ 本图递归/森林组合对照 · ✓ Henge 的 spinney 覆盖</p>{parametric_html(r)}</details></article>''')
    colrows=[]
    for j,c in enumerate(data['counterterm_columns'],1):
        entries='; '.join(f"{x['parent']} γ{x['edge_ids']}: {x['coefficient']}" for x in c['entries'])
        colrows.append(f'<tr><td>C{j}: {c["counterterm_type"]}</td><td>{svg(Graph.from_dict(c["cograph"]))}</td><td>{escape(entries)}</td><td>{c["total"]}</td><td>{c["expected_coefficient"]}</td></tr>')
    checkrows=''.join(f'<tr><td>L{x["loops"]} E{x["legs"]}</td><td>{x["pairings_visited"]}</td><td>{escape(str(x["weights"]))}</td><td>{"通过" if x["ok"] else "失败"}</td></tr>' for x in data['wick_validation'])
    summary=''.join(f'<tr><td>{x["loops"]}</td><td>{x["legs"]}</td><td>{x["count"]}</td><td>{x["weight_sum"]}</td></tr>' for x in data['summaries'])
    examples=''.join(f'<details><summary>{name}：{len(x["terms"])} 个实际减除项</summary><p class="formula">{escape(" ".join(x["operators"]))}</p><p>有理采样值：<code>{x["sum"]}</code></p><p class="muted">这只是声明路由下的被积函数采样，不是积分值或收敛证明。完整路由和实际分母见 examples.py。</p></details>' for name,x in data['representative_integrands'].items())
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>φ³ Lab · BPH 递归与参数积分</title>
<style>:root{--ink:#163b45;--accent:#176856;--bg:#f1f4f2;--line:#d4dfd9}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.7 -apple-system,BlinkMacSystemFont,'PingFang SC',sans-serif}header,main{max-width:1240px;margin:auto;padding:28px}h1{font-size:34px;margin:5px 0}h2{margin-top:40px}h3{margin-top:0}h4{margin-bottom:6px}p{margin:8px 0 14px}.muted,.meta{font-size:14px;color:#63766e}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:18px}.card,.panel{padding:23px;border:1px solid var(--line);background:white;border-radius:12px;margin-bottom:18px}.card svg{width:100%;max-height:260px}.card h3{display:flex;justify-content:space-between;gap:8px}.badge{font-size:13px;border:1px solid var(--line);padding:2px 8px;border-radius:30px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}.formula{font:17px/2 Georgia,serif;overflow-wrap:anywhere}.edge,.leg{stroke:#234b57;stroke-width:2;fill:none}.edge.highlight{stroke:#dc6a22;stroke-width:5}.vertex{fill:#234b57}.ct{fill:#f2ebd9;stroke:#926b1e}.ct-text{font-size:9px;text-anchor:middle;fill:#55401b}svg text{text-anchor:middle;font:13px system-ui;fill:#34525b}.edge-label{font-size:11px;paint-order:stroke;stroke:white;stroke-width:4;stroke-linejoin:round}.vertex-label{font-size:10px;fill:#80918c}button,select,input{font:inherit;border:1px solid var(--line);background:white;border-radius:6px;padding:6px 11px;color:var(--ink)}button,summary{cursor:pointer}summary{font-weight:600;color:var(--accent)}.good{color:var(--accent)}.filters{display:flex;gap:14px;flex-wrap:wrap;position:sticky;top:0;background:#f1f4f2f7;padding:14px 0;z-index:1}.columns svg{width:180px;height:130px}code{overflow-wrap:anywhere;font-size:13px}a{color:var(--accent)}.warning{border-left:4px solid #bf8f30;padding-left:15px}@media(max-width:600px){header,main{padding:16px}.grid{grid-template-columns:1fr}.columns{overflow:auto}h1{font-size:28px}}@media print{.filters{display:none}details{display:block}main{padding:0}.card{break-inside:avoid}}</style>
<header><p class="muted">PHI3LAB · 可复核两圈基准 · Python 3 标准库</p><h1>φ³ Lab：BPH 递归与积分互校</h1><p>固定标签 · 截肢 · 连通 1PI · 有质量 Euclidean φ³ · D=DIMENSION</p><p><a href="analysis.json">完整 JSON</a> · <a href="forests.tex">LaTeX 算子表</a> · <a href="../docs/modernization.html">现代文献与源码</a> · <a href="../README.md">使用说明</a></p><p class="warning">已完成六维两圈重叠二点图的递归减除、有限积分和本图解析界。其余图集保留精确组合数据；本软件不声称实现任意图、任意圈数的积分或一般证明。</p></header>
<main><section class="panel"><h2>图集总览</h2><table><tr><th>L</th><th>E</th><th>固定标签图数</th><th>Σ 1/S</th></tr>SUMMARY</table><p class="muted">tadpole、自环和平行边均保留；线性反项的有限重整化条件未在此选定。</p></section>
<h2>逐图检查</h2><div class="filters"><label>圈数 <select id="loops"><option value="all">全部</option><option>1</option><option>2</option></select></label><label>外腿 <select id="legs"><option value="all">全部</option><option>1</option><option>2</option><option>3</option></select></label><input id="search" placeholder="按图编号查找，如 L2E2"><span id="count"></span></div><div class="grid">CARDS</div>
<section class="panel"><h2>独立 Wick 配对检查</h2><p>枚举带标签半边的配对，独立比较按图自同构得到的 1/S。此检查不使用“图生成成功”等同于系数正确。</p><table><tr><th>图族</th><th>配对数</th><th>各图权重</th><th>结果</th></tr>CHECKS</table><p class="muted">若未运行 --wick，此表为空；请执行 validate 子命令取得独立配对验证。</p></section>
<section class="panel columns"><h2>低阶反项插入的列重组</h2><p>反项顶点已含自身 1/Sγ；一处裸图嵌入贡献 Sγ/S<sub>G</sub>。把同一带类型商图的项相加，检查是否等于该反项图自身的 1/S。C2 同时保留质量与动能两个局域结构。</p><table><tr><th>列</th><th>带类型的商图</th><th>父图与嵌入系数</th><th>列和</th><th>商图 1/S</th></tr>COLUMNS</table><p class="muted">此表核对两圈的一圈反项插入权重；整体两圈反项为各行 −t<sub>G</sub>R̄I 的加权和，未求其积分系数。不同 Henge 区域中同一 spinney 的出现次数不会被当成额外权重。</p></section>
<section class="panel"><h2>三个实际减除被积函数</h2>EXAMPLES<p>模块提供嵌套二点、嵌套三点和重叠二点例子，以及通用的传播子乘积 Taylor 系数递推。没有为全部森林自动生成任意路由的数值被积函数。</p></section>
<section class="panel"><h2>研究边界</h2><ul>LIMITS</ul><p>出处：<a href="https://arxiv.org/abs/hep-th/9612113">Kennedy 1996</a>；<a href="https://arxiv.org/abs/0712.1016">Kennedy–Binoth–Rippon 2007</a>；<a href="https://arxiv.org/abs/1402.2613">Borinsky 2014</a>；<a href="https://arxiv.org/abs/1506.03357">Gracey 2015</a>。实现与具体源码复用情况见调研文件。</p></section></main>
<script>function filter(){let n=0;const l=document.getElementById('loops').value,e=document.getElementById('legs').value,q=document.getElementById('search').value.toLowerCase();for(const c of document.querySelectorAll('.card')){const ok=(l==='all'||c.dataset.loops===l)&&(e==='all'||c.dataset.legs===e)&&c.dataset.id.toLowerCase().includes(q);c.hidden=!ok;if(ok)n++;}document.getElementById('count').textContent=n+' 张图';}for(const id of ['loops','legs','search'])document.getElementById(id).addEventListener('input',filter);document.addEventListener('click',ev=>{const b=ev.target.closest('[data-highlight]');if(!b)return;const ids=b.dataset.highlight.split(',');for(const edge of b.closest('.card').querySelectorAll('[data-eid]'))edge.classList.toggle('highlight',ids.includes(edge.dataset.eid));});filter();</script></html>'''
    evidence=external_evidence(Path(path).parent)
    rows=''.join(f'<tr><td><a href="{e["file"]}">{e["file"]}</a></td><td class="{e["status"]}">{e["description"]}</td></tr>' if e['status']!='not_run' else f'<tr><td>{e["file"]}</td><td>{e["description"]}</td></tr>' for e in evidence)
    poles=''.join(f'<tr><td>{b["id"]}</td><td>{escape(b["known_closed_form"])}</td><td>{b["simple_pole_coefficient"]}</td><td>通过</td></tr>' for b in data['one_loop_pole_benchmarks'])
    modern=f'''<section class="panel modern"><p class="eyebrow">PHI3LAB {data['software_version']} · 研究更新 2026-10-19</p><h2>以现代计算作为验证基准</h2>
<p>2025 年 φ³ 的 MS 重整化函数已推进到六圈。当前两圈软件的目标是有质量 BPH 的可复核计算链：从边嵌入子图和 BPH 递归减除，连接到参数积分与独立工具。</p>
<p><a href="../docs/modernization.html">阅读 2021—2026 文献与软件更新 →</a></p>
<div class="milestones"><div><strong>{len(data['graphs'])}</strong><span>精确 U / F 图</span></div><div><strong>pySecDec</strong><span>逐图输入与运行脚本</span></div><div><strong>3</strong><span>一圈解析 pole 基准</span></div></div>
<p class="muted">方法新颖性尚未确立；本图的两圈有限积分和解析界见上方专节，通用接口本身不构成方法创新。</p></section>
<section class="panel"><h2>本次源码的独立核查</h2><table><tr><th>证据</th><th>状态</th></tr>{rows}</table><p><a href="provenance.json">构建清单与文件校验值</a> · <a href="pysecdec/README.txt">运行后端的条件与归一化</a></p><p class="muted">pySecDec 核对范围为 U/F、Gamma 因子和指数；尚未生成 C++、编译或数值积分。核心计算仍只用 Python 标准库。</p></section>'''
    scope='<section class="panel"><h2>项目方法：BPH 递归重整化</h2><p>以 Henge 与 Bogoliubov 的递归 R 操作为主；本项目尚未核验的 Z 部分不作为证明前提。森林记法仅用于辅助展示和有限图集组合对照。软件测试不能替代一般收敛证明。</p><p><a href="../docs/bph-scope.html">导师的命名说明与本项目边界</a></p></section>'
    overlap_section=overlap_panel(data)
    if data.get('two_loop_finite') and overlap_section:
        start=overlap_section.index('<h3>在本机重新计算</h3>')
        end=overlap_section.index('<h3>解析上界与紫外收敛</h3>')
        overlap_section=overlap_section[:start]+overlap_section[end:]
        overlap_section='<details class="panel"><summary>重叠图的递归表达式与专项推导</summary>'+overlap_section+'</details>'
    page=page.replace('<main>', '<main>'+two_loop_panel(data)+overlap_section+scope+modern,1)
    if data['finite_one_loop_benchmarks']:
        finite_rows=''.join(f'<tr><td>{b["name"]}</td><td>{b["value"]:.15g}</td><td>{b["last_refinement_difference"]:.2g}</td></tr>' for b in data['finite_one_loop_benchmarks'])
        finite_section=f'''<section class="panel" id="finite-benchmarks"><h2>实际完成的一圈有限 BPH 积分</h2><p>在 m²=s<sub>i</sub>=1 的 Euclidean 点，计算 bubble 的 (1−t<sub>p</sub><sup>2</sup>)J 和 triangle 的 (1−t<sup>0</sup>)J。</p><table><tr><th>图</th><th>有限数值</th><th>32→64 阶差</th></tr>{finite_rows}</table><p>测度 d⁶k/π³，不含 1/S 与耦合。使用本地 Gauss–Legendre 参数积分；另与精确 Beta/Dirichlet 矩的独立收敛级数核对。</p><p class="muted">阶数差是稳定性诊断，不是严格误差上界。这些是本地一圈结果；两圈重叠图结果见上方专节。pySecDec 后端数值积分尚未执行。<a href="../docs/one-loop-benchmarks.html">推导、约定与复查方法</a>。</p></section>'''
        page=page.replace('<section class="panel"><h2>三个实际减除被积函数</h2>',finite_section+'<section class="panel"><h2>三个实际减除被积函数</h2>')
    page=page.replace('<section class="panel"><h2>三个实际减除被积函数</h2>',f'<section class="panel"><h2>一圈裸积分的解析极点基准</h2><p>测度 d<sup>d</sup>k / π<sup>d/2</sup>，d={data["dimension"]}−2ε；不含图权重和耦合。由 F 多项式的 simplex 矩积分与闭式公式独立核对。</p><table><tr><th>图</th><th>1/ε 系数</th><th>记录采样值</th><th>核查</th></tr>{poles}</table><p class="muted">此表计算裸 DR 极点；未求裸 Laurent 展开的有限系数。另一表单列 Taylor 减除后的有限 BPH 积分。记录样点 m²=1，bubble 的 p₁,E²=1/2。</p></section><section class="panel"><h2>三个实际减除被积函数</h2>')
    page=page.replace('<input id="search"', '<label>UV 结构 <select id="uv"><option value="all">全部</option><option value="primitive">无真 UV 子图</option><option value="subdivergent">含子图</option><option value="overlap">重叠发散</option></select></label><input aria-label="按图编号查找" id="search"')
    page=page.replace('<span id="count"></span>','<button id="reset">重置</button><span id="count" role="status" aria-live="polite"></span>')
    script='''<script>
const ids=['loops','legs','uv','search'];
function filter(){let n=0;const l=document.getElementById('loops').value,e=document.getElementById('legs').value,u=document.getElementById('uv').value,q=document.getElementById('search').value.trim().toLowerCase();for(const c of document.querySelectorAll('.card')){const ok=(l==='all'||c.dataset.loops===l)&&(e==='all'||c.dataset.legs===e)&&(u==='all'||(u==='subdivergent'?c.dataset.uv!=='primitive':c.dataset.uv===u))&&c.dataset.id.toLowerCase().includes(q);c.hidden=!ok;if(ok)n++;}document.getElementById('count').textContent=n+' / TOTAL 张图';}
for(const id of ids)document.getElementById(id).addEventListener('input',filter);
document.getElementById('reset').addEventListener('click',()=>{for(const id of ids)document.getElementById(id).value=id==='search'?'':'all';filter();});
document.addEventListener('click',ev=>{const b=ev.target.closest('[data-highlight],[data-clear-highlight]');if(!b)return;const card=b.closest('.card');const selected=b.hasAttribute('data-highlight')&&b.getAttribute('aria-pressed')!=='true';const edges=selected?b.dataset.highlight.split(','):[];for(const edge of card.querySelectorAll('[data-eid]'))edge.classList.toggle('highlight',edges.includes(edge.dataset.eid));for(const button of card.querySelectorAll('[data-highlight]'))button.setAttribute('aria-pressed',button===b&&selected?'true':'false');});
function locate(){const id=decodeURIComponent(location.hash.slice(1));const card=document.getElementById(id);if(!card||!card.classList.contains('card'))return;document.getElementById('reset').click();card.querySelector('details').open=true;requestAnimationFrame(()=>card.scrollIntoView({block:'start'}));}
window.addEventListener('hashchange',locate);filter();locate();
const form=document.getElementById('overlap-form');
if(form){let current=JSON.parse(document.getElementById('overlap-initial').textContent);
form.addEventListener('submit',async ev=>{ev.preventDefault();const status=document.getElementById('overlap-live-status'),out=document.getElementById('overlap-live-result'),button=document.getElementById('overlap-calculate');
if(location.protocol==='file:'){status.textContent='请先双击软件目录中的启动文件，开启本地计算服务。';return;}
button.disabled=true;status.textContent='正在本机计算并比较两个积分阶数…';
try{const query=new URLSearchParams({s:document.getElementById('overlap-s').value,m2:document.getElementById('overlap-m2').value,order:document.getElementById('overlap-order').value});const response=await fetch('/api/overlap?'+query);if(!response.headers.get('content-type')?.includes('application/json'))throw new Error('当前为静态预览，请双击启动文件开启本地计算服务。');const r=await response.json();if(!response.ok)throw new Error(r.error||'计算未完成');current=r;out.hidden=false;out.textContent='p² = '+r.s+'，m² = '+r.m2+'，阶数 '+r.order+'\\n有限积分 = '+r.value+'\\n乘图权重 1/2 = '+r.weighted_value+'\\n阶数差 = '+r.refinement_difference+'\\n振幅的解析上界 = '+(r.analytic_bounds.absolute_integral_bound??r.analytic_bounds.decimal_upper_bounds.absolute_integral_bound);status.textContent=(r.refinement_stable?'阶数比较通过。':'阶数比较尚未通过，请提高积分阶数。')+' 用时 '+r.seconds.toFixed(2)+' 秒；阶数差不代表严格误差界。';}
catch(error){status.textContent=error.message;}finally{button.disabled=false;}});
document.getElementById('overlap-download').addEventListener('click',()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(current,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='phi3lab-overlap-result.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});}
</script>'''.replace('TOTAL',str(len(data['graphs'])))
    script=script.replace('</script>',TWO_LOOP_SCRIPT+'</script>')
    page=page[:page.index('<script>')]+script+'</html>'
    for a,b in [('DIMENSION',str(data['dimension'])),('SUMMARY',summary),('CARDS',''.join(cards)),
                ('CHECKS',checkrows),('COLUMNS',''.join(colrows)),('EXAMPLES',examples),
                ('LIMITS',''.join('<li>'+escape(x)+'</li>' for x in data['limits']))]:page=page.replace(a,b)
    insertion=''.join(f'<p>L{x["loops"]}：{x["insertion_positions"]} 个逐线插入位置 → {x["three_point_graphs"]} 张带标签三点图，权重匹配：{x["ok"]}。</p>' for x in data.get('zero_leg_insertion_validation',[]))
    page=page.replace('<h2>研究边界</h2>', '<h2>裸图零动量外腿插入互校</h2>'+insertion+'<p>这里只验证图与权重；两点和三点的 Taylor 阶数不同，不能据此断言重整化后的质量导数恒等式成立。</p><h2>研究边界</h2>')
    page=page.replace('</style>', '.card{min-width:0;scroll-margin-top:100px}details{overflow-x:auto}button[aria-pressed=true]{background:#e2efe8;border-color:#176856}.modern,.overlap{border-top:4px solid #176856}.modern h2,.overlap h2{margin-top:6px}.eyebrow{font-size:12px;letter-spacing:.06em}.milestones,.calculation-inputs{display:flex;gap:20px;flex-wrap:wrap;margin:20px 0}.calculation-inputs label{display:grid;gap:6px;flex:1;min-width:160px}.calculation-inputs input,.calculation-inputs select{width:100%}.result-value{font:32px/1.6 ui-monospace,monospace;overflow-wrap:anywhere}#overlap-live-result{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf4ee;padding:16px;font-size:14px}#overlap-calculate{background:#176856;color:white}#overlap-calculate:disabled{opacity:.6}.milestones strong{display:block;font-size:25px}.milestones span{font-size:13px;color:#63766e}.passed{color:#176856}.failed,.stale,.invalid{color:#9c591c}@media(max-width:600px){.grid{grid-template-columns:minmax(0,1fr)}.filters{position:static}.card,.panel{padding:16px}.milestones{gap:18px}table{font-size:12px}td,th{padding:7px 4px}input{min-width:0;width:100%}.card h3{flex-wrap:wrap}.card{scroll-margin-top:16px}.result-value{font-size:24px}}'+ '</style>',1)
    page=page.replace('</style>','[hidden]{display:none!important}.table-scroll{overflow-x:auto}#two-loop-live-result{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf4ee;padding:16px;font-size:14px}#two-loop-calculate{background:#176856;color:white}#two-loop-graph{display:block;width:100%;margin-top:6px}.all-two-loop{border-top:4px solid #176856}.all-two-loop h2{margin-top:6px}</style>',1)
    page=page.replace('已完成六维两圈重叠二点图的递归减除、有限积分和本图解析界。其余图集保留精确组合数据；本软件不声称实现任意图、任意圈数的积分或一般证明。','已按统一零动量 BPH 递归方案扩展六维两圈图集，逐图列出有限积分及解析界。此处不外推为任意圈数的一般定理。')
    page=page.replace('tadpole、自环和平行边均保留；线性反项的有限重整化条件未在此选定。','tadpole、自环和平行边均保留；本计算采用完整零动量 Taylor 条件，一点图余项为零，裸 tadpole 并不为零。')
    page=page.replace('方法新颖性尚未确立；本图的两圈有限积分和解析界见上方专节，通用接口本身不构成方法创新。','方法新颖性尚未确立；支持图集的两圈有限积分和拓扑解析界见上方专节，通用接口本身不构成方法创新。')
    page=page.replace('pySecDec 核对范围为 U/F、Gamma 因子和指数；尚未生成 C++、编译或数值积分。核心计算仍只用 Python 标准库。','U/F 与约定核查、上游实际数值互校分别记录。原始三点图的 native 对照见 native_primitive_comparison.json；是否对应当前源码，以上表状态为准。核心可用 Python 标准库运行，NumPy 为可选加速。')
    page=page.replace('这些是本地一圈结果；两圈重叠图结果见上方专节。pySecDec 后端数值积分尚未执行。','这些是本地一圈结果；全部两圈结果见上方工作台。独立上游数值证据另列，不把本地 BPH 计算冒充上游输出。')
    project=Path(__file__).resolve().parents[1]
    report_dir=Path(path).resolve().parent
    automatic_page=report_dir/'automatic/index.html'
    if not automatic_page.is_file():automatic_page=project/'results/automatic/index.html'
    if automatic_page.is_file():
        target=quote(os.path.relpath(automatic_page,report_dir).replace(os.sep,'/'))
        page=page.replace('</header>','<p><a href="'+target+'">新增：从图构造 BPH 有限积分 · 一至四圈 698 张图 →</a></p></header>',1)
    for label,target in [('../docs/research.html',project/'docs/research.html'),('../README.md',project/'README.md'),('../docs/modernization.html',project/'docs/modernization.html'),('../docs/one-loop-benchmarks.html',project/'docs/one-loop-benchmarks.html'),('../docs/bph-scope.html',project/'docs/bph-scope.html'),('../docs/overlap-analysis.html',project/'docs/overlap-analysis.html'),('../analysis-project/main.pdf',project/'analysis-project/main.pdf'),('../full-analysis-project/main.pdf',project/'full-analysis-project/main.pdf')]:
        page=page.replace('href="'+label,'href="'+quote(os.path.relpath(target,report_dir).replace(os.sep,'/')))
    for name in ('two-loop-scheme.html','iteration-ledger.json','primitive-analysis.html','nested-bubble-analysis.html','nested-triangle-analysis.html','extreme-scales.html'):
        target=project/'docs'/name
        page=page.replace('href="../docs/'+name,'href="'+quote(os.path.relpath(target,report_dir).replace(os.sep,'/')))
    Path(path).write_text(page,encoding='utf-8')
