"""Path B: expand the counterterm action, Wick pair, then glue local sources.

NO imports of graphs.generate/divergent_subgraphs/contract, ct_bph,
renormalization, insertion_columns, or target tables. New two-loop local
coefficients are fixed by the zero-momentum condition on this path's own
bare and lower-counterterm action diagrams.
"""
from fractions import Fraction
from itertools import permutations
from math import factorial
from .graphs import Graph, canonical
from .ct_graph import identity, marked_identity, isomorphism, momentum_transport
from .ct_wick import enumerate_wick, action_sectors, SPECS


def glue(host,vertex,source,port_permutation):
    boundary=[]
    for label,v in host.external:
        if v==vertex:boundary.append(('external',label,0))
    for i,(a,b) in enumerate(host.edges):
        if a==vertex:boundary.append(('edge',i,0))
        if b==vertex:boundary.append(('edge',i,1))
    if len(boundary)!=len(source.external):raise ValueError('Arity mismatch')
    if sorted(port_permutation)!=list(range(len(boundary))):raise ValueError('Invalid port bijection')
    old={v:i for i,v in enumerate(v for v in range(host.n) if v!=vertex)}
    offset=host.n-1
    source_attach={label:v+offset for label,v in source.external}
    targets={port:source_attach[source.external[port_permutation[i]][0]] for i,port in enumerate(boundary)}
    edges=[]
    for i,(a,b) in enumerate(host.edges):
        edges.append((targets[('edge',i,0)] if a==vertex else old[a],
                      targets[('edge',i,1)] if b==vertex else old[b]))
    embedded=list(range(len(edges),len(edges)+len(source.edges)))
    edges.extend((a+offset,b+offset) for a,b in source.edges)
    ext=tuple((label,targets[('external',label,0)] if v==vertex else old[v]) for label,v in host.external)
    parent=Graph(len(old)+source.n,tuple(edges),ext)
    _,cg=canonical(parent); mapping=isomorphism(parent,cg)
    return cg,{'boundary_slots':boundary,'source_ports_in_slot_order':[source.external[j][0] for j in port_permutation],
               'raw_expanded_graph':parent.to_dict(),'expanded_to_canonical':mapping,
               'loop_momentum_change':momentum_transport(parent,cg,mapping),
               'embedded_edges':[mapping['edges'][str(i)] for i in embedded],
               'host_edges_to_parent':{str(i):mapping['edges'][str(i)] for i in range(len(host.edges))},
               'source_edges_to_parent':{str(i):mapping['edges'][str(e)] for i,e in enumerate(embedded)}}


def build_path_b():
    sectors=action_sectors(); rows={}; columns={}; sector_records=[]; audits={}; sources={}
    def row(g):
        gid=identity(g)
        if gid not in rows:rows[gid]={'id':gid,'graph':g.to_dict(),'terms':[]}
        return rows[gid]
    for legs in (1,2,3):
        d=enumerate_wick(('phi3',)*legs,tuple(range(1,legs+1)))
        sources[legs]=[r for r in d['records'] if r['accepted']]
        audits[f'source_L1E{legs}']=d
    local_host=None
    for sector in sectors:
        kinds=sector['kinds']; d=enumerate_wick(kinds)
        key='+'.join(kinds); audits[key]=d
        accepted=[r for r in d['records'] if r['accepted']]
        sector_records.append({**sector,'action_monomial':key,'candidate_classes':len(d['records']),
                               'accepted_graphs':[r['id'] for r in accepted],
                               'reason':None if accepted else 'no_connected_1PI_Wick_graph; tree_edges_or_C1_leaf_are_bridges'})
        for host_record in accepted:
            host=Graph.from_dict(host_record['graph']); host_weight=Fraction(host_record['positive_weight'])
            ct_vertices=[v for v,k in enumerate(host.vertex_kinds) if k!='phi3']
            if not ct_vertices:
                row(host)['terms'].append({'column':'bare','operator':f'I:{identity(host)}','coefficient':str((-1)**len(kinds)*host_weight),
                    'kind':'bare','wick_record':host_record,'sign_factors':{'bare_vertices':(-1)**len(kinds)}})
                continue
            if host.vertex_kinds==('C3L2',):
                local_host=host_record; continue
            if len(ct_vertices)!=1:raise ValueError('Unexpected accepted multiple CT host')
            v=ct_vertices[0];op=host.vertex_kinds[v];arity=SPECS[op][0]
            if SPECS[op][1]!=1:raise ValueError('Unexpected local countervertex at this order')
            cid=identity(host)
            for source_record in sources[arity]:
                source=Graph.from_dict(source_record['graph']); sw=Fraction(source_record['positive_weight'])
                source_sign=(-1)**source.n * -1
                host_bare_sign=(-1)**(host.n-1)
                permutations_count=factorial(arity)
                columns[cid]={'id':cid,'operator_type':op,'host':host.to_dict(),'source':identity(source),
                              'host_wick':host_record,'source_wick':source_record,'source_weight':str(sw),
                              'host_weight':str(host_weight),'source_vertex_sign':source_sign,
                              'source_action_sign':-source_sign,'port_bijections':permutations_count}
                for permutation in permutations(range(arity)):
                    g,mapping=glue(host,v,source,permutation)
                    mark=marked_identity(g,mapping['embedded_edges'])
                    coeff=host_bare_sign*source_sign*host_weight*sw/permutations_count
                    row(g)['terms'].append({'column':cid,'operator':f'{op}:{identity(source)}:{cid}:{mark}',
                        'coefficient':str(coeff),'kind':'insertion','operator_type':op,'marked_id':mark,
                        'embedded_edges':mapping['embedded_edges'],'gluing':mapping,
                        'coefficient_factors':{'host_weight':str(host_weight),'source_weight':str(sw),
                          'port_average':str(Fraction(1,permutations_count)),
                          'host_bare_sign':host_bare_sign,'countervertex_sign':source_sign},
                        'host_wick_id':host_record['id'],'source_wick_id':source_record['id']})
    if local_host is None:raise ArithmeticError('Action did not generate the two-loop local vertex')
    for r in rows.values():
        # Independent renormalization condition Gamma_R^(3)(0)=(-g): the
        # order-g^5 coefficient at zero vanishes. This fixes -delta_g^(2).
        for t in list(r['terms']):
            r['terms'].append({'column':'C3L2','operator':f'T0[{t["operator"]}]',
                'coefficient':str(-Fraction(t['coefficient'])*Fraction(local_host['positive_weight'])),
                'kind':'overall','projection_of':t['operator'],'local_host_wick':local_host['id'],
                'derivation':'minus zero-momentum projection of this path B action contribution'})
    return {'path':'B','generation':'typed counterterm action -> independent Wick -> all port bijections',
            'local_host':local_host,'rows':list(rows.values()),'columns':list(columns.values()),
            'sectors':sector_records,'audits':audits}
