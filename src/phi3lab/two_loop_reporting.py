"""Offline all-graph results and the corresponding local calculation form."""
import html
import json


FAMILIES = {
    'one_point_local': ('一点图：局域完整减除', 'two-loop-scheme.html'),
    'overlap': ('重叠三角子图', 'overlap-analysis.html'),
    'nested_bubble_two_point': ('二点图：bubble 插入', 'nested-bubble-analysis.html'),
    'nested_bubble_three_point': ('三点图：bubble 插入', 'nested-bubble-analysis.html'),
    'primitive_three_point': ('三点图：无真发散子图', 'primitive-analysis.html'),
    'nested_triangle': ('三点图：triangle 插入', 'nested-triangle-analysis.html'),
}


def two_loop_panel(data):
    records = data.get('two_loop_finite', [])
    if not records:
        return ''
    rows, options = [], []
    selected = 'L2E2G1' if any(r['graph'] == 'L2E2G1' for r in records) else records[0]['graph']
    for r in records:
        title, proof = FAMILIES[r['family']]
        bound = r['analytic_bounds'].get('absolute_integral_bound')
        bound_text = format(bound, '.8g') if isinstance(bound, (int, float)) else '见解析式'
        rows.append(f'''<tr><td><a href="#{r['graph']}">{r['graph']}</a><br><small>{html.escape(title)}</small></td>
<td>{r['value']:.15g}</td><td>{r['symmetry_weight']}</td><td>{r['weighted_value']:.15g}</td><td>{bound_text}<br><a href="../docs/{proof}">推导</a></td></tr>''')
        options.append(f'<option value="{r["graph"]}"'+(' selected' if r['graph'] == selected else '')+f'>{r["graph"]} · {html.escape(title)}</option>')
    payload = json.dumps(records, ensure_ascii=False).replace('<', '\\u003c')
    return f'''<section class="panel all-two-loop" id="two-loop-workbench"><p class="eyebrow">全部两圈图 · BPH 递归 · 同一减除条件</p>
<h2>逐图有限积分与解析界</h2><p>此表为六维、有共同正质量的计算结果。二点图取 p²=1，三点图取 p₁²=p₂²=p₃²=1，质量平方 m²=1。测度为 d⁶k/π³ d⁶q/π³；顶点因子 (−g)<sup>E+2</sup> 与 (4π)<sup>−6</sup> 单列，未暗中乘入数值。</p>
<div class="table-scroll"><table><tr><th>图与拓扑</th><th>无权重 J<sub>R</sub></th><th>1/S</th><th>乘图权重后</th><th>|J<sub>R</sub>| 上界</th></tr>{''.join(rows)}</table></div>
<p class="muted">一点图的零来自所声明的完整 Taylor 减除条件；裸 tadpole 并不为零。上界控制振幅及收敛，阶数差仅诊断求积稳定性。<a href="../docs/two-loop-scheme.html">统一方案与覆盖范围</a> · <a href="../full-analysis-project/main.pdf">全部拓扑的完整证明 PDF</a> · <a href="../docs/extreme-scales.html">极端尺度改进与独立验证</a> · <a href="../docs/iteration-ledger.json">迭代记录</a>。</p>
<h3>选择一张图，在本机重新计算</h3><form id="two-loop-form"><label>两圈图<select id="two-loop-graph">{''.join(options)}</select></label>
<div class="calculation-inputs"><label id="two-loop-label-s1">外动量平方 p₁²<input id="two-loop-s1" type="number" min="0" step="any" value="1" required></label>
<label id="two-loop-label-s2">外动量平方 p₂²<input id="two-loop-s2" type="number" min="0" step="any" value="1" required></label>
<label id="two-loop-label-s3">外动量平方 p₃²<input id="two-loop-s3" type="number" min="0" step="any" value="1" required></label>
<label>共同质量平方 m²<input id="two-loop-m2" type="number" min="0" step="any" value="1" required></label>
<label>每维积分阶数<input id="two-loop-order" type="number" min="8" max="32" step="1" value="16" required></label></div>
<p id="two-loop-kinematic-note" class="muted"></p><button type="submit" id="two-loop-calculate">计算这张图</button> <button type="button" id="two-loop-download">下载已显示结果 JSON</button></form>
<p id="two-loop-live-status" role="status" aria-live="polite" class="muted">下方显示预先生成的实际计算；启动本地服务后可更改参数重算。</p>
<pre id="two-loop-live-result"></pre><script type="application/json" id="two-loop-initial">{payload}</script></section>'''


TWO_LOOP_SCRIPT = r'''
const allForm=document.getElementById('two-loop-form');
if(allForm){
 const benchmarks=JSON.parse(document.getElementById('two-loop-initial').textContent);
 const byGraph=new Map(benchmarks.map(r=>[r.graph,r]));
 let current;
 const input=id=>document.getElementById('two-loop-'+id);
 const status=input('live-status'),out=input('live-result');
 function show(r){current=r;const b=r.analytic_bounds;
 const bound=b.absolute_integral_bound??b.decimal_upper_bounds?.absolute_integral_bound??b.formula??'见本图推导';
 out.textContent=r.graph+'\n外动量不变量 = '+JSON.stringify(r.external_invariants)+'，m² = '+r.m2+
 '\n无权重有限积分 = '+r.value+'\n乘图权重 '+r.symmetry_weight+' = '+r.weighted_value+
 '\n振幅的解析上界 = '+bound+'\n积分阶数 = '+r.order+'，阶数差 = '+r.refinement_difference+
 '\n测度 d⁶k/π³ d⁶q/π³；耦合、顶点符号和 (4π)⁻⁶ 单列。'+
 (r.kinematics_status?.boundary_projection?'\n输入位于浮点 Gram 边界，已使用上列有效物理不变量；原请求为 '+JSON.stringify(r.requested_external_invariants):'');
 if(r.extreme_scale_method){out.textContent+='\n极端尺度计算：'+r.extreme_scale_method;}
 if(r.precision_scope){out.textContent+='\n'+r.precision_scope.message;}
 if(r.precision_diagnostics){out.textContent+='\n当前求积路线 = '+r.precision_diagnostics.production_route+'（'+r.precision_diagnostics.integration_dimension+' 维）';out.textContent+='\n两种坐标求积的差 = '+r.precision_diagnostics.coordinate_route_difference+
 '\n两种细化外推的差 = '+r.precision_diagnostics.extrapolated_route_difference+'（诊断，非严格误差界）';}}
 function selectGraph(){const r=byGraph.get(input('graph').value);const n=r.external_invariants.length;
 for(let j=1;j<=3;j++){input('label-s'+j).hidden=j>n;input('s'+j).disabled=j>n;input('s'+j).value='1';}
 input('m2').value='1';input('order').value=r.order;input('order').max=r.order_limits.max;
 input('order').disabled=n===0;
 input('kinematic-note').textContent=n===3?'三条外动量流入并满足守恒；三个平方须构成合法 Euclidean Gram 矩阵。':n===1?'p₁ = −p₂，只需一个动量平方。':'一点图的外动量由守恒固定为零，完整 Taylor 减除后的余项精确为零。';
 show(r);status.textContent='已显示本图预先生成的实际计算。修改参数后点击计算。';}
 input('graph').addEventListener('change',selectGraph);
 for(const name of ['s1','s2','s3','m2','order'])input(name).addEventListener('input',()=>{status.textContent='输入已修改；下方仍是上一次结果，点击计算后更新。';});
 allForm.addEventListener('submit',async ev=>{ev.preventDefault();
 if(location.protocol==='file:'){status.textContent='请先双击软件目录中的启动文件，开启本地计算服务。';return;}
 const graph=input('graph').value,r0=byGraph.get(graph),n=r0.external_invariants.length;
 const q=new URLSearchParams({graph,m2:input('m2').value});
 if(n)q.set('order',input('order').value);
 for(let j=1;j<=n;j++)q.set('s'+j,input('s'+j).value);
 const controls=Array.from(allForm.querySelectorAll('input,select,button'));
 const disabledBefore=controls.map(control=>control.disabled);
 controls.forEach(control=>control.disabled=true);status.textContent='正在本机计算两个阶数，请保持页面开启…';
 try{const response=await fetch('/api/two-loop?'+q);
 if(!response.headers.get('content-type')?.includes('application/json'))throw new Error('当前是静态预览，请双击启动文件开启本地服务。');
 const r=await response.json();if(!response.ok)throw new Error(r.error||'计算未完成');show(r);
 status.textContent=r.precision_scope?r.precision_scope.message:r.family==='one_point_local'?'完整 Taylor 减除的精确恒等式；无需数值积分。':
 (r.refinement_stable?'两个阶数的数值比较通过。':(r.order>=r.order_limits.max?'已达本图阶数上限，细化尚未通过；当前精度未确认，请查看独立复算记录。':'两个阶数尚未达到设定的一致性阈值，请提高积分阶数。'))+' 用时 '+r.seconds.toFixed(2)+' 秒。阶数差不是严格误差界。';
 }catch(error){status.textContent=error.message;}finally{controls.forEach((control,index)=>control.disabled=disabledBefore[index]);}});
 input('download').addEventListener('click',()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(current,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='phi3lab-'+current.graph+'-BPH.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
 selectGraph();
}
'''
