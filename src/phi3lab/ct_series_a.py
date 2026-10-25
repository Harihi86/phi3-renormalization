"""Path A: actual recursive spinneys, with typed local Taylor projectors."""
from fractions import Fraction
from itertools import product
from functools import lru_cache
from .graphs import Graph,divergent_subgraphs,subgraph,contract,generate as low_generate
from .higher_graphs import generate as high_generate,symmetry,canonical
from .renormalization import spinneys,forests
from .ct_identity import graph_id,mapping
from .ct_bph import taylor_trace,henge_trace
from .ct_series_model import operators,vertex_type,aggregate


@lru_cache(maxsize=8192)
def recursive_expansion(g):
    subs=divergent_subgraphs(g);terms=[{'nodes':[],'top':[],'children':[]}]
    for spinney in spinneys(subs):
        if not spinney:continue
        alternatives=[]
        for index in spinney:
            sub=subs[index];sg,ports=subgraph(g,sub['edge_ids']);child=recursive_expansion(sg)
            choices=[]
            for c in child['terms']:
                for op in operators(len(sg.external)):
                    nested=[{**n,'edges':[sub['edge_ids'][i] for i in n['edges']]} for n in c['nodes']]
                    root={'edges':sub['edge_ids'],'op':op,'loops':sg.loops}
                    choices.append({'nodes':nested+[root],'top':{**sub,'op':op},
                                    'child':{'source_id':graph_id(sg),'source_graph':sg.to_dict(),'boundary_ports':ports,
                                             'taylor':taylor_trace(g,sub),'inner_terms':c}})
            alternatives.append(choices)
        for selected in product(*alternatives):
            terms.append({'nodes':[n for x in selected for n in x['nodes']],
                          'top':[x['top'] for x in selected],'children':[x['child'] for x in selected]})
    return {'graph':g.to_dict(),'id':graph_id(g),'proper_parts':subs,'terms':terms}


def build_family(loops,legs,include_henges=True):
    rows=[];columns={}
    for index,original in enumerate((low_generate if loops<=2 else high_generate)(loops,legs),1):
        _,g=canonical(original);node=recursive_expansion(g);weight=Fraction(1,symmetry(g));terms=[]
        for expansion in node['terms']:
            tops=expansion['top'];nodes=expansion['nodes'];host_id='bare';trace={}
            if tops:
                raw=contract(g,tops)
                kinds=tuple(vertex_type(s['op'],s['loops']) for s in tops)+raw.vertex_kinds[len(tops):]
                raw=Graph(raw.n,raw.edges,raw.external,kinds);_,host=canonical(raw);host_id=graph_id(host)
                columns[host_id]={'id':host_id,'host':host.to_dict(),'local_vertices':[{'kind':k,'vertex':v} for v,k in enumerate(host.vertex_kinds) if k!='phi3']}
                trace={'contraction_mapping':mapping(raw,host),'top_subgraphs':tops,
                       'taylor':[taylor_trace(g,s) for s in tops],'recursive_children':expansion['children']}
            coefficient=((-1)**(g.n+len(nodes)))*weight
            terms.append({'column':host_id,'nodes':nodes,'coefficient':str(coefficient),'kind':'bare' if not nodes else 'insertion',
                          'factors':{'bare_sign':(-1)**g.n,'subtraction_sign':(-1)**len(nodes),'parent_weight':str(weight)},**trace})
        lower=aggregate(g,terms);overall=[]
        for op in operators(legs):
            col=vertex_type(op,loops)
            for t in lower:
                overall.append({'column':col,'nodes':t['nodes']+[{'edges':list(range(len(g.edges))),'op':op,'loops':loops}],
                                'coefficient':str(-Fraction(t['coefficient'])),'kind':'overall','projection_of':t['symbol'],
                                'taylor':{'operator':op,'degree':2 if legs==2 else 0,'operator_degree':2 if op=='C2z' else 0,
                                          'external_vectors':[f'p{i}' for i in range(1,legs)],
                                          'retained_loops':[f'q{i+1}' for i in range(loops)]}})
        supports={frozenset(tuple(n['edges']) for n in t['nodes']) for t in terms}
        oracle={frozenset(tuple(node['proper_parts'][j]['edge_ids']) for j in f) for f in forests(node['proper_parts'])}
        if supports!=oracle:raise ArithmeticError('Actual recursion differs from forest oracle')
        hs=henge_trace(g) if include_henges else None
        if hs is not None:
            from .ct_verify import check_henge
            if not check_henge(g,hs):raise ArithmeticError('Henge mapping mismatch')
        rows.append({'id':graph_id(g),'old_id':f'L{loops}E{legs}G{index}','graph':g.to_dict(),'symmetry':symmetry(g),
                     'positive_weight':str(weight),'recursive_ast':node,'terms':lower+aggregate(g,overall),
                     'henges':hs,'henge_mapping_verified':hs is not None,'forest_oracle_equal':True})
    return {'path':'A','loops':loops,'legs':legs,'rows':rows,'columns':list(columns.values()),
            'overall_columns':[vertex_type(op,loops) for op in operators(legs)],
            'method':'actual recursive C=-T Rbar over vertex-disjoint spinneys; forest oracle only afterwards'}
