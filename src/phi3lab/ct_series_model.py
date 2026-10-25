"""Shared representation of typed Taylor forests, not either generation path."""
from collections import defaultdict
from functools import lru_cache
import hashlib,json
from .graphs import Graph
from .higher_graphs import canonical
from .ct_identity import graph_id,mapping
from .ct_graph import momentum_transport

def operators(legs):return {1:('C1',),2:('C2m','C2z'),3:('C3',)}.get(legs,())
def vertex_type(op,loops):return op if loops==1 else f'{op}L{loops}'

@lru_cache(maxsize=65536)
def _signature(g,nodes):
    # Incidence encoding treats parallel propagators as distinct edge vertices.
    # A projector node is connected to every edge in its support; containment
    # and overlap are thus preserved without enumerating vertex relabellings.
    edges=[];kinds=list(g.vertex_kinds)
    for i,(a,b) in enumerate(g.edges):
        v=len(kinds);kinds.append('propagator');edges.extend([(a,v),(b,v)])
    for op,loop,ids in nodes:
        v=len(kinds);kinds.append(f'projector:{op}:L{loop}')
        edges.extend((v,g.n+i) for i in ids)
    aux=Graph(len(kinds),tuple(edges),g.external,tuple(kinds))
    key,_=canonical(aux)
    return 'forest_'+hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()

def signature(g,nodes):
    return _signature(g,tuple(sorted((n['op'],n['loops'],tuple(sorted(n['edges']))) for n in nodes)))

def canonicalize(g,nodes):
    _,cg=canonical(g);m=mapping(g,cg)
    mapped=[{**n,'edges':sorted(m['edges'][str(i)] for i in n['edges'])} for n in nodes]
    return cg,mapped,{'raw_graph':g.to_dict(),'mapping':m,'loop_momentum_change':momentum_transport(g,cg,m)}

def glue_raw(host,vertex,source,source_nodes,ports):
    from .routing import momentum_routing
    boundary=[('external',l,0) for l,v in host.external if v==vertex]
    for i,(a,b) in enumerate(host.edges):
        if a==vertex:boundary.append(('edge',i,0))
        if b==vertex:boundary.append(('edge',i,1))
    if len(boundary)!=len(source.external) or sorted(ports)!=list(range(len(boundary))):raise ValueError('Port bijection mismatch')
    keep=[v for v in range(host.n) if v!=vertex];vm={v:i for i,v in enumerate(keep)};offset=len(keep)
    targets={slot:source.external[ports[j]][1]+offset for j,slot in enumerate(boundary)}
    edges=[(targets[('edge',i,0)] if a==vertex else vm[a],targets[('edge',i,1)] if b==vertex else vm[b]) for i,(a,b) in enumerate(host.edges)]
    edge_offset=len(edges);edges.extend((a+offset,b+offset) for a,b in source.edges)
    external=tuple((l,targets[('external',l,0)] if v==vertex else vm[v]) for l,v in host.external)
    raw=Graph(len(keep)+source.n,tuple(edges),external,tuple(host.vertex_kinds[v] for v in keep)+source.vertex_kinds)
    nodes=[{**n,'edges':[edge_offset+i for i in n['edges']]} for n in source_nodes]
    route=momentum_routing(host)
    boundary_vectors=[route['external_momenta'][str(i)] if kind=='external' else
                      [(-1 if end==0 else 1)*x for x in route['edge_momenta'][str(i)]] for kind,i,end in boundary]
    if any(sum(r[i] for r in boundary_vectors) for i in range(len(route['basis']))):raise ArithmeticError('Host boundary momentum is not conserved')
    trace={'host':host.to_dict(),'vertex':vertex,'source':source.to_dict(),'source_id':graph_id(source),
           'source_nodes':source_nodes,'ports':list(ports),'boundary_slots':boundary,
           'host_vertices_to_expanded':vm,'source_vertices_to_expanded':{v:v+offset for v in range(source.n)},
           'source_edges_to_expanded':{i:i+edge_offset for i in range(len(source.edges))},'result':raw.to_dict(),
           'host_routing':route,'source_routing':momentum_routing(source),
           'source_boundary_in_host_basis':{str(source.external[ports[j]][0]):vec for j,vec in enumerate(boundary_vectors)},
           'taylor_rule':'Project only source boundary dependence; retain all host loop momenta and propagators',
           'operator_basis':{'C2m':'P^0','C2z':'P^2; coefficient of total vector degree two, not degree two in P^2','C3':'P^0'}}
    return raw,nodes,trace

def aggregate(g,terms):
    from fractions import Fraction
    grouped={}
    for t in terms:
        key=(t['column'],signature(g,t['nodes']))
        if key not in grouped:grouped[key]={**t,'coefficient':Fraction(),'contributions':[]}
        r=grouped[key];r['coefficient']+=Fraction(t['coefficient'])
        r['contributions'].append({k:v for k,v in t.items() if k not in ('column','nodes','kind')})
    return [{**r,'coefficient':str(r['coefficient']),'symbol':symbol} for (_,symbol),r in grouped.items() if r['coefficient']]
