"""Path A: actual bare-graph Bogoliubov recursion with Henge and Taylor traces."""
from fractions import Fraction
from itertools import combinations
from functools import lru_cache
from .graphs import (Graph, generate, canonical, symmetry, divergent_subgraphs,
                     subgraph, contract, henges, connected)
from .automatic_bph import compile_bph, _compatible_routing
from .renormalization import spinneys, forests, structural_checks
from .routing import momentum_routing
from .ct_graph import marked_identity, isomorphism, external_orbit
from .ct_identity import graph_id as identity


def local_operators(legs):
    return {1:['C1'],2:['C2m','C2z'],3:['C3']}.get(legs,[])


@lru_cache(maxsize=8192)
def recursion(g):
    """C(G)=-t_G [I(G)+sum_S prod C(gamma) I(G/S)], calculated recursively.

    The AST is formed from the input graph, never from forest words. Forests
    are only an independently enumerated check on the resulting supports.
    """
    subs=divergent_subgraphs(g); terms=[{'spinney':[],'children':[], 'cograph':g.to_dict()}]
    for s in spinneys(subs):
        if not s: continue
        selected=[subs[i] for i in s]
        children=[]
        for i in s:
            sg,ports=subgraph(g,subs[i]['edge_ids'])
            children.append({'embedding':subs[i], 'boundary_ports':ports, 'node':recursion(sg)})
        terms.append({'spinney':list(s),'children':children,'cograph':contract(g,selected).to_dict()})
    omega=6*g.loops-2*len(g.edges)
    return {'graph':g.to_dict(),'id':identity(g),'omega':omega,'proper_parts':subs,
            'rbar':terms,'counterterm':{'operation':'minus_Taylor','degree':omega,'argument':'rbar'} if omega>=0 else None,
            'renormalized':{'operation':'one_minus_Taylor','degree':omega,'argument':'rbar'} if omega>=0 else {'operation':'identity','argument':'rbar'}}


def recursive_supports(node):
    result={frozenset()}
    for term in node['rbar'][1:]:
        choices=[frozenset()]
        for child in term['children']:
            parent_edges=child['embedding']['edge_ids']
            nested=recursive_supports(child['node'])
            mapped=[frozenset([tuple(parent_edges)]) | frozenset(tuple(parent_edges[e] for e in ids) for ids in word)
                    for word in nested]
            choices=[a|b for a in choices for b in mapped]
        result.update(choices)
    return result


def taylor_trace(g, sub):
    sg,ports=subgraph(g,sub['edge_ids'])
    rows,chords,tree=_compatible_routing(g,[sub])
    names=[f'q{i+1}' for i in range(g.loops)]+[f'p{k}' for k,_ in g.external[:-1]]
    # Original physical source momenta, expressed in the compatible basis.
    unit=lambda j:[int(i==j) for i in range(len(names))]
    physical={k:unit(g.loops+i) for i,(k,_) in enumerate(g.external[:-1])}
    physical[g.external[-1][0]]=[-sum(row[i] for row in physical.values()) for i in range(len(names))]
    boundary=[]
    for port in ports:
        typ,index,end,vertex=port
        coefficients=physical[index] if typ=='external' else [(-1 if end==0 else 1)*x for x in rows[index]]
        boundary.append({'port':list(port),'coefficients':coefficients})
    if any(sum(b['coefficients'][i] for b in boundary) for i in range(len(names))):
        raise ArithmeticError('Boundary momentum not conserved')
    inner=[j for j,e in enumerate(chords) if e in sub['edge_ids']]
    outer=[j for j in range(g.loops) if j not in inner]
    return {'subgraph':sg.to_dict(),'boundary_ports':boundary,'parent_basis':names,
        'parent_edge_routing':rows,'tree_edges':tree,'chord_edges':chords,
        'inner_loop_basis':[names[i] for i in inner],'outer_loop_basis':[names[i] for i in outer],
        'independent_boundary_symbols':[f'P{i+1}' for i in range(len(ports)-1)],
        'dependent_boundary':'P_last = -sum(P_independent)',
        'subgraph_own_routing':momentum_routing(sg),
        'degree':sub['omega'],'subtraction_point':'all independent boundary P=0; m^2>0 fixed',
        'operation':'sum_{k=0}^omega [t^k] I_gamma(k_inner,t*P); co-graph is not scaled',
        'outside_propagators_retained':[i for i in range(len(g.edges)) if i not in sub['edge_ids']],
        'operators':local_operators(len(ports)),
        'two_point_projection':{'mass':'I_gamma(0)',
          'linear':'Odd under simultaneous inversion of internal loop momenta; zero after regulated integration',
          'kinetic':'p^2/(2D) sum_mu partial_mu partial_mu I_gamma(0), D=6; degree 2 in p, degree 1 in p^2'} if len(ports)==2 else None,
        'outside_loop_as_boundary':any(any(b['coefficients'][i] for i in outer) for b in boundary),
        'regulator':'common UV regulator until subtractions combined; no numerical counterterm integral asserted'}


def henge_trace(g):
    out=[]
    for h in henges(g):
        vertex_map={v:i for i,b in enumerate(h['blocks']) for v in b['vertices']}
        skeleton_original=h['bridges']+[h['minimum_edge']]
        records=[]
        for b in h['blocks']:
            item=dict(b)
            if b['loops']:
                sg,ports=subgraph(g,b['edge_ids'])
                item.update(graph=sg.to_dict(),boundary_ports=ports,
                            recursive_node=recursion(sg),
                            zero_reason=None if b['needs_uv_subtraction'] else 'overall_power_counting_convergent_block_retained')
            records.append(item)
        out.append({**h,'blocks':records,'deleted_edge':list(g.edges[h['minimum_edge']]),
                    'remaining_original_edges':[i for i in range(len(g.edges)) if i!=h['minimum_edge']],
                    'skeleton_original_edges':skeleton_original,
                    'original_vertex_to_skeleton':vertex_map,
                    'original_edge_map':{str(i):({'skeleton_edge':skeleton_original.index(i)} if i in skeleton_original else
                         {'block':next(j for j,b in enumerate(h['blocks']) if i in b['edge_ids']),
                          'local_edge':next(b['edge_ids'].index(i) for b in h['blocks'] if i in b['edge_ids'])}) for i in range(len(g.edges))},
                    'sector_rule':'e_min is the least momentum-norm line; tie by original edge id; Henge blocks recursively renormalized',
                    'multiplicity_rule':'Sectors partition the integral; never multiply a subtraction by its Henge occurrence count'})
    return out


def build_path_a():
    rows=[]; columns={}
    for old,g in enumerate(generate(2,3),1):
        gid=identity(g); s=symmetry(g); node=recursion(g); terms=[]
        def add(column,operator,coefficient,**trace):
            terms.append({'column':column,'operator':operator,'coefficient':str(coefficient),**trace})
        add('bare',f'I:{gid}',-Fraction(1,s),kind='bare',coupling_power=5)
        for term in node['rbar'][1:]:
            # The enumerated two-loop E3 scope has primitive one-loop children.
            # Reject unsupported products instead of silently flattening them.
            if len(term['children'])!=1: raise ValueError('Two-loop E3 spinney unexpectedly non-singleton')
            child=term['children'][0]; sub=child['embedding']; sg=Graph.from_dict(child['node']['graph'])
            if child['node']['proper_parts']:raise ValueError('Unexpected nested child in the two-loop E3 scope')
            for op in local_operators(sub['external_legs']):
                raw_host=contract(g,[sub]); kinds=tuple(op if k.startswith('C') else k for k in raw_host.vertex_kinds)
                typed=Graph(raw_host.n,raw_host.edges,raw_host.external,kinds); _,host=canonical(typed)
                cid=identity(host); source=identity(sg); mark=marked_identity(g,sub['edge_ids'])
                columns[cid]={'id':cid,'operator_type':op,'host':host.to_dict(),'source':source,'source_weight':str(Fraction(1,symmetry(sg)))}
                add(cid,f'{op}:{source}:{cid}:{mark}',Fraction(1,s),kind='insertion',
                    source=source,source_symmetry=symmetry(sg),parent_symmetry=s,
                    parent_embedding_factor=str(Fraction(symmetry(sg),s)),
                    embedded_edges=sub['edge_ids'],marked_id=mark,operator_type=op,
                    host=host.to_dict(),contraction_mapping=isomorphism(typed,host),
                    taylor=taylor_trace(g,sub),sign_factors={'bare_vertices':-1,'subtraction':-1})
        lower=list(terms)
        for t in lower:
            add('C3L2',f'T0[{t["operator"]}]',-Fraction(t['coefficient']),kind='overall',
                projection_of=t['operator'],sign_factors={'rbar_coefficient':t['coefficient'],'overall_subtraction':-1},
                taylor={'degree':0,'variables':['p1','p2'],'p3':'-p1-p2','retained_loop_variables':['q1','q2']})
        hs=henge_trace(g)
        forest_supports={frozenset(tuple(node['proper_parts'][j]['edge_ids']) for j in f) for f in forests(node['proper_parts'])}
        checks=structural_checks(g,node['proper_parts'],hs)
        checks['new_recursive_ast_equals_forest_oracle']=recursive_supports(node)==forest_supports
        rows.append({'id':gid,'old_id':f'L2E3G{old}','graph':g.to_dict(),'symmetry':s,'positive_weight':str(Fraction(1,s)),
                     'bare_amplitude_coefficient':str(-Fraction(1,s)),'coupling_power':5,
                     'external_permutation_orbit':external_orbit(g),'unfixed_id':identity(g,False),
                     'recursive_ast':node,'henges':hs,'executable_taylor_program':compile_bph(g).to_dict(),
                     'checks':checks,'terms':terms})
    return {'path':'A','generation':'degree-constrained bare adjacency; actual recursive C=-T Rbar',
            'rows':rows,'columns':list(columns.values())}
