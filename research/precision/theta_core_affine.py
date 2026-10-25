"""Exact original-edge loop changes for decorated K2,3 kernels."""
from itertools import product,combinations
import sympy as s

def verify(raw,tree):
    gd=raw['gaussian_full_descriptor'];L=gd['loops'];edges=gd['graph']['edges'];gram=s.Matrix([[s.Rational(v) for v in row] for row in raw['external_gram']]);parts=gd['proper_subgraphs'];z=s.symbols('z:'+str(len(parts)));nodes=[];parents={};parent_arcs={}
    def arcs(node):return node['arcs'] if node['kind']=='cycle' else [a for path in node['paths'] for a in path['arcs']]
    def visit(node,parent):
        for child in node['children']:
            rep=min(v for e in child['edge_ids'] for v in edges[e]);parent_arcs[id(child)]=next(a for a in arcs(node) if rep in a['child_representatives']);visit(child,node)
        if parent is not None:nodes.append(node);parents[id(node)]=parent
    visit(tree,None);assert len(nodes)==L-2;choices=[[a['edges'][0] for a in n['arcs']] for n in nodes];root_edges=[path['arcs'][0]['edges'][0] for path in tree['paths'][:2]];records=[]
    def same(a,b):return a==b or a==[-v for v in b]
    for term in gd['terms']:
        variables=term['jet_subgraphs'];word=term['subgraphs_inner_first'];assert term['jet_orders']==[2]*len(variables) and term['sign']==(-1)**len(word)
        route=[[s.Rational(v['coefficient'])*s.prod(z[i]**power for i,power in zip(variables,v['powers'])) for v in row] for row in term['scaled_routing']];found=None
        for chosen_children in product(*choices):
            for signs in product((-1,1),repeat=2):
                chosen=[*chosen_children,*root_edges];rows=[route[e] for e in chosen_children]+[[sg*v for v in route[e]] for sg,e in zip(signs,root_edges)];change=s.Matrix([r[:L] for r in rows]);det=s.factor(change.det())
                if det not in (-1,1):continue
                inv=change.inv();offset=s.Matrix([r[L:] for r in rows]);converted=[]
                for row in route:
                    loop=s.Matrix([row[:L]])*inv;ext=s.Matrix([row[L:]])-loop*offset;converted.append([s.expand(v) for v in list(loop)+list(ext)])
                ok=True
                for i,node in enumerate(nodes):
                    own=[0]*(L+2);own[i]=1;selected=next(a for a in node['arcs'] if chosen[i] in a['edges'])
                    for arc in node['arcs']:
                        flow=converted[arc['edges'][0]];ok &= all(same(converted[e],flow) for e in arc['edges'])
                        if arc is selected:ok &= same(flow,own)
                        else:
                            part=node['source_part_index'];scale=z[part] if part in word else 1;parent=converted[parent_arcs[id(node)]['edges'][0]]
                            ok &= any(same(flow,[own[j]+sg*scale*parent[j] for j in range(L+2)]) for sg in (-1,1))
                shifts=[]
                for j,path in enumerate(tree['paths']):
                    own=[0]*(L+2)
                    if j<2:own[L-2+j]=1
                    else:own[L-2]=own[L-1]=-1
                    first,second=path['arcs'];ok &= all(same(converted[e],own) for e in first['edges']);flow=converted[second['edges'][0]];pivot=next(i for i,v in enumerate(own[:L]) if v);sg=s.factor(flow[pivot]/own[pivot])
                    if sg not in (-1,1):ok=False;break
                    ok &= all(flow[i]==sg*own[i] for i in range(L)) and all(same(converted[e],flow) for e in second['edges']);shifts.append(s.Matrix([[v/sg for v in flow[L:]]]))
                if ok:
                    ok=len(shifts)==3 and all(not(v.free_symbols&set(z)) for row in shifts for v in row)
                    ok &= sum(shifts,s.zeros(1,2))==s.zeros(1,2)
                    ok &= all(s.expand((v*gram*v.T)[0])==1 for v in shifts)
                    ok &= all(s.expand((a*gram*b.T)[0])==s.Rational(-1,2) for a,b in combinations(shifts,2))
                if ok:
                    found={'word':word,'sign':term['sign'],'chosen_original_edges':chosen,'root_row_signs':signs,'jacobian_determinant':int(det),'selected_rows':[[str(v) for v in r] for r in rows],'converted_routing':[[str(v) for v in r] for r in converted],'root_external_shifts':[[str(v) for v in r] for r in shifts],'all_parent_boundary_momenta_retained':True,'sum_root_branch_momenta_zero':True,'external_Gram_diagonal':1,'external_Gram_offdiagonal':'-1/2','extra_routing_multiplicity':1};break
            if found:break
        assert found is not None,('No original-edge affine witness',word);records.append(found)
    return records
