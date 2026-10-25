"""Path B: typed action orbit sums and independent recursive CT source library.

No bare A catalogue, divergent subgraph finder, contractions or A coefficients.
Only this path's own lower-loop subtraction conditions define its CT sources.
"""
from fractions import Fraction
from itertools import permutations
from math import factorial
from .graphs import Graph
from .ct_typed import enumerate_action,sectors,specifications
from .ct_identity import graph_id
from .ct_series_model import operators,vertex_type,glue_raw,canonicalize,aggregate


class ActionLibrary:
    def __init__(self):self.families={};self.sources={}

    def family(self,loops,legs):
        key=(loops,legs)
        if key in self.families:return self.families[key]
        specs=specifications();rows={};columns={};audit=[];local={};used_sources={}
        def row(g):
            gid=graph_id(g)
            if gid not in rows:rows[gid]={'id':gid,'graph':g.to_dict(),'terms':[]}
            return rows[gid]
        for sector in sectors(loops,legs):
            d=enumerate_action(sector['kinds'],tuple(range(1,legs+1)))
            accepted=[r for r in d['records'] if r['accepted']]
            audit.append({**sector,'candidate_classes':len(d['records']),'accepted':len(accepted),
                          'exclusions':{reason:sum(r['reason']==reason for r in d['records']) for reason in {r['reason'] for r in d['records']}},
                          'algorithm':d['algorithm'],'denominator':d['species_factorial']*d['vertex_factorial'],
                          'records':d['records'],'raw_terms':d['raw_terms']})
            for rec in accepted:
                host=Graph.from_dict(rec['graph']);cts=[v for v,k in enumerate(host.vertex_kinds) if k!='phi3']
                if not cts:
                    row(host)['terms'].append({'column':'bare','nodes':[],'coefficient':str((-1)**host.n*Fraction(rec['positive_weight'])),
                                               'kind':'bare','wick_record':rec});continue
                if host.n==1 and specs[host.vertex_kinds[0]]['loops']==loops:
                    local[host.vertex_kinds[0]]=rec;continue
                if any(specs[host.vertex_kinds[v]]['arity']==1 for v in cts):
                    raise ValueError('A C1 1PI host requires an explicit unresolved finite-condition treatment')
                cid=graph_id(host);columns[cid]={'id':cid,'host':host.to_dict(),'host_wick':rec,
                    'local_vertices':[{'kind':host.vertex_kinds[v],'vertex':v} for v in cts]}
                def expand(current,nodes,coefficient,events):
                    vertices=[v for v,k in enumerate(current.vertex_kinds) if k!='phi3']
                    if not vertices:
                        g,ns,canonical_map=canonicalize(current,nodes)
                        row(g)['terms'].append({'column':cid,'nodes':ns,'coefficient':str(coefficient),'kind':'insertion',
                            'host_wick_id':cid,'host_weight':rec['positive_weight'],'host_bare_sign':(-1)**sum(k=='phi3' for k in host.vertex_kinds),
                            'gluing_events':events,'canonical_map':canonical_map});return
                    v=vertices[0];kind=current.vertex_kinds[v];spec=specs[kind]
                    if spec['loops']>=loops:raise ValueError('Recursive CT source did not decrease loop order')
                    self.family(spec['loops'],spec['arity'])
                    for source in self.sources[(spec['loops'],spec['arity'],spec['operator'])]:
                        sg=Graph.from_dict(source['graph']);average=Fraction(1,factorial(spec['arity']))
                        source_key=graph_id(sg)+':'+source['symbol']
                        used_sources[source_key]={'source_id':graph_id(sg),'symbol':source['symbol'],'coefficient':source['coefficient'],
                            'nodes':source['nodes'],'loops':spec['loops'],'legs':spec['arity'],'operator':spec['operator'],
                            'projection_of':source['projection_of'],'family_record':f'library/L{spec["loops"]}E{spec["arity"]}.json'}
                        for ports in permutations(range(spec['arity'])):
                            raw,new_nodes,trace=glue_raw(current,v,sg,source['nodes'],ports)
                            trace.update(source_symbol=source['symbol'],source_key=source_key,source_coefficient=source['coefficient'],port_average=str(average))
                            expand(raw,nodes+new_nodes,coefficient*Fraction(source['coefficient'])*average,events+[trace])
                expand(host,[],Fraction(rec['positive_weight'])*(-1)**sum(k=='phi3' for k in host.vertex_kinds),[])
        if set(local)!={vertex_type(op,loops) for op in operators(legs)}:raise ArithmeticError('Missing local action vertex')
        for r in rows.values():
            g=Graph.from_dict(r['graph']);r['terms']=aggregate(g,r['terms']);lower=list(r['terms']);overall=[]
            for op in operators(legs):
                col=vertex_type(op,loops);w=Fraction(local[col]['positive_weight'])
                if w!=1:raise ArithmeticError('Unexpected local host weight')
                pieces=[]
                for t in lower:
                    piece={'column':col,'nodes':t['nodes']+[{'edges':list(range(len(g.edges))),'op':op,'loops':loops}],
                           'coefficient':str(-w*Fraction(t['coefficient'])),'kind':'overall','projection_of':t['symbol'],
                           'derivation':'minus local projection of this path B action sum','local_wick':local[col]}
                    pieces.append(piece)
                for t in aggregate(g,pieces):
                    self.sources.setdefault((loops,legs,op),[]).append({'graph':g.to_dict(),**t})
                    overall.append(t)
            r['terms']+=overall
        result={'path':'B','loops':loops,'legs':legs,'rows':list(rows.values()),'columns':list(columns.values()),
                'overall_columns':list(local),'sectors':audit,'source_library':used_sources,'method':'independent typed action orbit counting, source expansion and recursive local conditions',
                'finite_linear_condition':'unspecified; no accepted C1 hosts in target family'}
        self.families[key]=result;return result
