"""Independent validation and exact cell comparison for recursive CT series."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from .graphs import Graph,subgraph,is_1pi
from .ct_identity import graph_id
from .ct_series_model import signature,operators,vertex_type,glue_raw,canonicalize


def validate_source_library(b,families):
    """Check references against B's saved lower-loop family projections."""
    for ref in b.get('source_library',{}).values():
        family=families[ref['loops'],ref['legs']]
        terms=[t for row in family['rows'] if row['id']==ref['source_id'] for t in row['terms']
               if t['symbol']==ref['symbol'] and t['column']==vertex_type(ref['operator'],ref['loops'])]
        if len(terms)!=1 or any(terms[0][k]!=ref[k] for k in ('coefficient','nodes','projection_of')):
            raise ValueError('Saved lower-loop source projection differs from reference')
    return True


def validate(data):
    rows=data['rows'];columns={c['id']:c for c in data['columns']}
    if len(columns)!=len(data['columns']) or len({r['id'] for r in rows})!=len(rows):raise ValueError('Duplicate identities')
    local=set(data['overall_columns'])
    if local!={vertex_type(op,data['loops']) for op in operators(data['legs'])}:raise ValueError('Local operator columns mismatch')
    for cid,c in columns.items():
        if graph_id(Graph.from_dict(c['host']))!=cid:raise ValueError('Host identity/type/label mismatch')
    for row in rows:
        g=Graph.from_dict(row['graph'])
        if graph_id(g)!=row['id'] or g.loops!=data['loops'] or len(g.external)!=data['legs']:raise ValueError('Parent identity or family mismatch')
        if data['path']=='A':
            from .ct_bph import taylor_trace
            def check_expansion(parent,exp):
                if len(exp['top'])!=len(exp['children']):raise ValueError('Recursive child count mismatch')
                for sub,child in zip(exp['top'],exp['children']):
                    source,ports=subgraph(parent,sub['edge_ids'])
                    if child['source_graph']!=source.to_dict() or child['source_id']!=graph_id(source):raise ValueError('Recursive source mismatch')
                    import json
                    if child['taylor']!=json.loads(json.dumps(taylor_trace(parent,sub))):raise ValueError('Recursive Taylor trace mismatch')
                    check_expansion(source,child['inner_terms'])
            for expansion in row['recursive_ast']['terms']:check_expansion(g,expansion)
        for t in row['terms']:
            if t['column'] not in local|{'bare'}|set(columns):raise ValueError('Unknown column')
            if signature(g,t['nodes'])!=t['symbol']:raise ValueError('Typed projector structure mismatch')
            parts=[];seen=set()
            for n in t['nodes']:
                ids=n['edges'];key=tuple(sorted(ids))
                if not ids or len(ids)!=len(set(ids)) or key in seen or any(type(i)is not int or i<0 or i>=len(g.edges) for i in ids):raise ValueError('Invalid or duplicate projector support')
                seen.add(key);sg,_=subgraph(g,ids)
                if sg.loops!=n['loops'] or n['op'] not in operators(len(sg.external)) or not is_1pi(range(sg.n),sg.edges):raise ValueError('Invalid local source or loop grade')
                parts.append((set(ids),{v for i in ids for v in g.edges[i]}))
            for (a,va),(b,vb) in combinations(parts,2):
                if not (a<b or b<a or va.isdisjoint(vb)):raise ValueError('Forbidden overlapping product')
            if t['kind']=='overall':
                roots=[n for n in t['nodes'] if len(n['edges'])==len(g.edges)]
                if len(roots)!=1 or vertex_type(roots[0]['op'],roots[0]['loops'])!=t['column']:raise ValueError('Overall projection column mismatch')
                if data['path']=='A':
                    expected={'operator':roots[0]['op'],'degree':2 if data['legs']==2 else 0,
                              'operator_degree':2 if roots[0]['op']=='C2z' else 0,
                              'external_vectors':[f'p{i}' for i in range(1,data['legs'])],
                              'retained_loops':[f'q{i+1}' for i in range(data['loops'])]}
                    if t['taylor']!=expected:raise ValueError('Overall Taylor variables or degree mismatch')
            contributions=t.get('contributions',[])
            if contributions and sum((Fraction(c['coefficient']) for c in contributions),Fraction())!=Fraction(t['coefficient']):raise ValueError('Contribution sum mismatch')
            for c in contributions:
                if 'gluing_events' not in c:continue
                host=Graph.from_dict(columns[t['column']]['host']);current=host;nodes=[]
                coefficient=Fraction(c['host_weight'])*c['host_bare_sign']
                for event in c['gluing_events']:
                    if current.to_dict()!=event['host']:raise ValueError('Broken gluing chain')
                    source=Graph.from_dict(event['source'])
                    if graph_id(source)!=event['source_id']:raise ValueError('Source identity mismatch')
                    ref=data['source_library'].get(event['source_key'])
                    if not ref or any(ref[k]!=event[v] for k,v in [('source_id','source_id'),('symbol','source_symbol'),('coefficient','source_coefficient'),('nodes','source_nodes')]):
                        raise ValueError('Counterterm source library reference mismatch')
                    if signature(source,event['source_nodes'])!=ref['symbol']:raise ValueError('Source forest identity mismatch')
                    current,new,expected=glue_raw(current,event['vertex'],source,event['source_nodes'],event['ports'])
                    for key,value in expected.items():
                        # JSON roundtrip equates tuple/list and integer-key dictionaries.
                        import json
                        if json.dumps(value,sort_keys=True)!=json.dumps(event[key],sort_keys=True):
                            if json.loads(json.dumps(value))!=json.loads(json.dumps(event[key])):raise ValueError('Gluing map mismatch: '+key)
                    nodes+=new;coefficient*=Fraction(event['source_coefficient'])*Fraction(event['port_average'])
                parent,ns,canonical_map=canonicalize(current,nodes)
                if graph_id(parent)!=row['id'] or signature(parent,ns)!=t['symbol']:raise ValueError('Expanded decorated parent mismatch')
                import json
                if json.loads(json.dumps(canonical_map))!=c['canonical_map']:raise ValueError('Canonical or momentum mapping mismatch')
                if coefficient!=Fraction(c['coefficient']):raise ValueError('Gluing factors mismatch')
    return True


def compare(a,b):
    # Round-trip first: the same verifier is used on memory and saved artifacts.
    import json
    a=json.loads(json.dumps(a));b=json.loads(json.dumps(b));validate(a);validate(b)
    ar={r['id']:r for r in a['rows']};br={r['id']:r for r in b['rows']}
    ac={c['id'] for c in a['columns']};bc={c['id'] for c in b['columns']}
    differences=[];cells=[]
    if set(ar)!=set(br):differences.append({'kind':'graph_set','A':sorted(ar),'B':sorted(br)})
    if ac!=bc:differences.append({'kind':'host_set','A':sorted(ac),'B':sorted(bc)})
    cols=['bare']+sorted(ac|bc)+sorted(set(a['overall_columns'])|set(b['overall_columns']))
    for gid in sorted(ar.keys()|br.keys()):
        for col in cols:
            values=[]
            for rs in (ar,br):
                accum=defaultdict(Fraction)
                for t in rs.get(gid,{}).get('terms',[]):
                    if t['column']==col:accum[t['symbol']]+=Fraction(t['coefficient'])
                values.append({k:str(v) for k,v in sorted(accum.items()) if v})
            ok=values[0]==values[1];zero=not values[0] and not values[1]
            cell={'row':gid,'column':col,'A':values[0],'B':values[1],'equal':ok,'zero':zero,
                  'zero_reason':'No A recursive embedding and no B action-source expansion' if zero else None}
            cells.append(cell)
            if not ok:differences.append(cell)
    return {'all_equal':not differences,'exact':True,'numeric_tolerance_used':False,'rows':len(ar),'columns':len(cols),
            'checked_cells':len(cells),'nonzero_cells':sum(not c['zero'] for c in cells),'cells':cells,'differences':differences}
