"""Exact, fail-closed comparison of the independent counterterm paths."""
from collections import defaultdict
from fractions import Fraction
from .graphs import Graph, canonical, symmetry
from .ct_graph import identity, marked_identity, momentum_transport, isomorphism


def validate_path(data):
    if len({r['id'] for r in data['rows']})!=len(data['rows']):raise ValueError('Duplicate graph row')
    columns={c['id']:c for c in data['columns']}
    for cid,c in columns.items():
        host=Graph.from_dict(c['host'])
        if identity(host)!=cid:raise ValueError('Countergraph identity/label/type mismatch')
        if c['operator_type'] not in host.vertex_kinds:raise ValueError('Local operator mismatch')
        if 'host_wick' in c and Fraction(c['host_weight'])!=Fraction(1,symmetry(host)):
            raise ValueError('Wick host weight differs from automorphisms')
    for r in data['rows']:
        g=Graph.from_dict(r['graph'])
        if identity(g)!=r['id']:raise ValueError('Bare graph identity/label/type mismatch')
        for t in r['terms']:
            Fraction(t['coefficient'])
            if t['column'] not in ('bare','C3L2') and t['column'] not in columns:
                raise ValueError('Missing host column')
            if t['kind']=='bare' and t['operator']!=f'I:{r["id"]}':raise ValueError('Bare symbol mismatch')
            if t['kind']=='insertion':
                if t['marked_id']!=marked_identity(g,t['embedded_edges']):raise ValueError('Embedded edge orbit mismatch')
                c=columns[t['column']]
                expected=f'{t["operator_type"]}:{c["source"]}:{c["id"]}:{t["marked_id"]}'
                if t['operator']!=expected:raise ValueError('Insertion symbol mismatch')
                if 'gluing' in t:
                    gl=t['gluing']; raw=Graph.from_dict(gl['raw_expanded_graph'])
                    if identity(raw)!=r['id']:raise ValueError('Glued parent mismatch')
                    m=isomorphism(raw,g)
                    if m!=gl['expanded_to_canonical']:raise ValueError('Half-edge map mismatch')
                    if momentum_transport(raw,g,m)!=gl['loop_momentum_change']:raise ValueError('Momentum map mismatch')
            if t['kind']=='overall' and t['operator']!=f'T0[{t["projection_of"]}]':
                raise ValueError('Overall projection mismatch')


def compare_paths(a,b):
    validate_path(a);validate_path(b)
    ar={r['id']:r for r in a['rows']}; br={r['id']:r for r in b['rows']}
    ac={c['id']:c for c in a['columns']}; bc={c['id']:c for c in b['columns']}
    differences=[];cells=[]
    if set(ar)!=set(br):differences.append({'reason':'parent_graph_set','A':sorted(ar),'B':sorted(br)})
    if set(ac)!=set(bc):differences.append({'reason':'countergraph_set','A':sorted(ac),'B':sorted(bc)})
    cols=['bare']+sorted(ac.keys()|bc.keys())+['C3L2']
    for gid in sorted(ar.keys()|br.keys()):
        for col in cols:
            sums=[]; trace=[]
            for data in (ar,br):
                terms=[t for t in data.get(gid,{}).get('terms',[]) if t['column']==col]
                accum=defaultdict(Fraction)
                for t in terms:accum[t['operator']]+=Fraction(t['coefficient'])
                sums.append({k:str(v) for k,v in sorted(accum.items()) if v})
                trace.append(terms)
            ok=sums[0]==sums[1]
            cell={'row':gid,'column':col,'A':sums[0],'B':sums[1],'equal':ok,
                  'zero':not sums[0] and not sums[1],
                  'zero_reason':'no embedding and no action gluing yields this labelled parent/typed host' if not sums[0] and not sums[1] else None,
                  'A_term_count':len(trace[0]),'B_term_count':len(trace[1])}
            cells.append(cell)
            if not ok:differences.append(cell)
    return {'exact':True,'all_equal':not differences,'rows':len(ar),'columns':len(cols),
            'checked_cells':len(cells),'nonzero_cells':sum(not c['zero'] for c in cells),
            'symbolic_monomials':sum(len(c['A']) for c in cells),'cells':cells,'differences':differences,
            'numeric_tolerance_used':False}


def bridge_oracle(vertices,edges):
    """Independent endpoint reachability, retaining parallel-edge identities."""
    bridges=[]
    for i,(start,end) in enumerate(edges):
        if start==end:continue
        seen={start};todo=[start]
        while todo:
            v=todo.pop()
            for j,(a,b) in enumerate(edges):
                if j==i:continue
                if a==v and b not in seen:seen.add(b);todo.append(b)
                if b==v and a not in seen:seen.add(a);todo.append(a)
        if end not in seen:bridges.append(i)
    return bridges


def check_henge(g,records):
    checks=[]
    for h in records:
        ids=h['remaining_original_edges'];edges=[g.edges[i] for i in ids]
        bridges=[ids[i] for i in bridge_oracle(range(g.n),edges)]
        checks.append(sorted(bridges)==sorted(h['bridges']))
        checks.append(set(h['original_edge_map'])==set(map(str,range(len(g.edges)))))
        for i in range(len(g.edges)):
            mapping=h['original_edge_map'][str(i)]
            if 'skeleton_edge' in mapping:
                a,b=g.edges[i]; vm=h['original_vertex_to_skeleton']
                get=lambda v:vm[v] if v in vm else vm[str(v)]
                checks.append(list(h['skeleton_edges'][mapping['skeleton_edge']])==[get(a),get(b)])
            else:checks.append(h['blocks'][mapping['block']]['edge_ids'][mapping['local_edge']]==i)
    return all(checks)
