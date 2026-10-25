"""Exact typed action expansion by adjacency orbits, without Wick-pair scans.

For each adjacency matrix, Gaussian pair multiplicity is the product of
remaining stub factorials divided by parallel-edge factorials and loop
2^l l!. Future identical vertices are ordered with an explicit multinomial
orbit multiplier, not silently discarded. No A catalogue or subgraphs enter.
"""
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import product
from math import factorial,prod
from functools import lru_cache
from .graphs import Graph,components,connected,is_1pi
from .higher_graphs import canonical


def specifications(max_loop=4):
    result={'phi3':{'arity':3,'loops':0,'coupling_power':1,'operator':'bare'}}
    for loop in range(1,max_loop+1):
        for op,arity in [('C1',1),('C2m',2),('C2z',2),('C3',3)]:
            name=op if loop==1 else f'{op}L{loop}'
            result[name]={'arity':arity,'loops':loop,'coupling_power':2*loop+arity-2,'operator':op}
    return result


def sectors(loops,legs):
    if type(loops)is not int or not 1<=loops<=4 or legs not in (1,2,3):raise ValueError('Use loops 1..4, legs 1..3')
    specs=specifications(loops);names=tuple(specs);power=2*loops+legs-2;out=[]
    def visit(i,left,grade,kinds):
        if i==len(names):
            if left or not kinds:return
            stubs=sum(specs[k]['arity'] for k in kinds)-legs
            if stubs<0 or stubs%2:return
            topology=stubs//2-len(kinds)+1
            if topology<0 or topology+grade!=loops:return
            out.append({'kinds':tuple(kinds),'coupling_power':power,'topological_loops':topology,
                        'counterterm_loops':grade,'total_loops':loops})
            return
        name=names[i];s=specs[name]
        for count in range(left//s['coupling_power']+1):
            if grade+count*s['loops']>loops:break
            visit(i+1,left-count*s['coupling_power'],grade+count*s['loops'],kinds+[name]*count)
    visit(0,power,0,[])
    return out


def attachments(kinds,labels,specs):
    pools={k:tuple(v for v,t in enumerate(kinds) if t==k) for k in set(kinds)}
    arities=[specs[k]['arity'] for k in kinds]
    for attach in product(range(len(kinds)),repeat=len(labels)):
        used=defaultdict(list);valid=True
        for v in attach:
            kind=kinds[v]
            if v not in used[kind]:
                if v!=pools[kind][len(used[kind])]:valid=False;break
                used[kind].append(v)
        counts=Counter(attach)
        if not valid or any(counts[v]>r for v,r in enumerate(arities)):continue
        orbit=prod(factorial(len(pools[k]))//factorial(len(pools[k])-len(vs)) for k,vs in used.items())
        stubs=prod(factorial(r)//factorial(r-counts[v]) for v,r in enumerate(arities))
        yield attach,tuple(r-counts[v] for v,r in enumerate(arities)),orbit,stubs


def adjacency_orbits(degree,colors):
    """Yield a degree solution and number of labelled matrices it represents."""
    n=len(degree);matrix=[[0]*n for _ in range(n)];remaining=list(degree)
    def row(i,multiplier):
        if i==n:
            if not any(remaining):yield tuple(tuple(r) for r in matrix),multiplier
            return
        original=remaining[i]
        groups=defaultdict(list)
        for v in range(i+1,n):groups[(colors[v],tuple(matrix[a][v] for a in range(i)))].append(v)
        previous={b:a for cell in groups.values() for a,b in zip(cell,cell[1:])}
        for loops in range(original//2+1):
            matrix[i][i]=loops;need=original-2*loops
            def fill(j,left):
                if j==n:
                    if left:return
                    factor=prod(factorial(len(cell))//prod(factorial(c) for c in Counter(matrix[i][v] for v in cell).values()) for cell in groups.values())
                    remaining[i]=0
                    yield from row(i+1,multiplier*factor)
                    remaining[i]=original
                    return
                maximum=min(left,remaining[j])
                if j in previous:maximum=min(maximum,matrix[i][previous[j]])
                future=sum(remaining[j+1:])
                for count in range(max(0,left-future),maximum+1):
                    matrix[i][j]=matrix[j][i]=count;remaining[j]-=count
                    yield from fill(j+1,left-count)
                    remaining[j]+=count
                matrix[i][j]=matrix[j][i]=0
            yield from fill(i+1,need)
        matrix[i][i]=0
    yield from row(0,1)


@lru_cache(maxsize=128)
def enumerate_action(kinds,labels=(1,2,3)):
    kinds=tuple(kinds);labels=tuple(labels);specs=specifications();n=len(kinds)
    if not kinds or any(k not in specs for k in kinds):raise ValueError('Unknown or empty action monomial')
    species=prod(factorial(c) for c in Counter(kinds).values())
    vertices=prod(factorial(specs[k]['arity']) for k in kinds)
    records={};raw=[];total=0
    for attach,degree,attachment_orbit,stub_choices in attachments(kinds,labels,specs):
        if sum(degree)%2:continue
        colors=tuple((kinds[v],tuple(l for l,w in zip(labels,attach) if w==v)) for v in range(n))
        for matrix,matrix_orbit in adjacency_orbits(degree,colors):
            denominator=prod(factorial(matrix[a][b])*(2**matrix[a][b] if a==b else 1)
                             for a in range(n) for b in range(a,n))
            numerator=prod(factorial(r) for r in degree)
            if numerator%denominator:raise ArithmeticError('Nonintegral Gaussian pairing multiplicity')
            pairings=numerator//denominator
            weight=pairings*matrix_orbit*attachment_orbit*stub_choices;total+=weight
            edges=tuple((a,b) for a in range(n) for b in range(a,n) for _ in range(matrix[a][b]))
            g=Graph(n,edges,tuple(zip(labels,attach)),kinds);key,cg=canonical(g)
            cpts=components(range(n),edges);vacuum=sum(not set(attach).intersection(c) for c in cpts)
            conn=len(cpts)==1;onepi=conn and is_1pi(range(n),edges)
            if key not in records:
                records[key]={'graph':cg.to_dict(),'wick_multiplicity':0,'positive_weight':'0','raw_ids':[],
                    'species_factorial':species,'vertex_factorial':vertices,'connected':conn,'one_pi':onepi,
                    'vacuum_components':vacuum,'self_loops':sum(a==b for a,b in edges),'accepted':onepi,
                    'reason':'accepted_amputated_1PI' if onepi else 'vacuum_factor' if vacuum else 'disconnected' if not conn else 'internal_bridge_1PR'}
            r=records[key];r['raw_ids'].append(len(raw));r['wick_multiplicity']+=weight
            raw.append({'graph':g.to_dict(),'attachment_orbit':attachment_orbit,'stub_attachment_choices':stub_choices,
                        'adjacency_orbit':matrix_orbit,'remaining_degrees':degree,'pairing_numerator':numerator,
                        'parallel_loop_denominator':denominator,'pairings_per_matrix':pairings,'wick_multiplicity':weight})
    for r in records.values():r['positive_weight']=str(Fraction(r['wick_multiplicity'],species*vertices))
    return {'kinds':kinds,'labels':labels,'records':[records[k] for k in sorted(records)],'raw_terms':raw,
            'full_wick_multiplicity':total,'species_factorial':species,'vertex_factorial':vertices,
            'algorithm':'degree-matrix orbits with explicit multinomial and Gaussian stub factorials',
            'source_source_pairs_included':False,
            'scope':'amputated source attachments; retains vacuum, disconnected, bridge and self-loop candidates',
            'pairings_explicitly_visited':0}


def enumerate_functional(kinds,labels):
    """Full Gaussian derivative audit, including direct source-source pairs.

    There are at most four partial source matchings for E<=3. The expensive
    part uses the same exact matrix-orbit count, never individual Wick pairs.
    These disconnected terms are audited before connected/amputated filtering.
    """
    from .ct_wick import source_matchings
    cases=[];total=0;counts=Counter()
    for pairs,attached in source_matchings(tuple(labels)):
        data=enumerate_action(tuple(kinds),tuple(attached));records=[]
        for r in data['records']:
            conn=not pairs and r['connected'];accepted=conn and r['one_pi']
            record={**r,'source_pairs':pairs,'connected':conn,'one_pi':accepted,'accepted':accepted,
                    'reason':'accepted_amputated_1PI' if accepted else 'vacuum_factor' if r['vacuum_components'] else 'disconnected' if not conn else 'internal_bridge_1PR'}
            records.append(record);counts['functional']+=1
            if not r['vacuum_components']:counts['normalized']+=1
            if conn:counts['connected_amputated']+=1
            if accepted:counts['one_pi']+=1
        total+=data['full_wick_multiplicity']
        cases.append({'source_pairs':pairs,'attached_labels':attached,'full_wick_multiplicity':data['full_wick_multiplicity'],
                      'records':records,'raw_terms':data['raw_terms']})
    stubs=sum(specifications()[k]['arity'] for k in kinds)+len(labels)
    expected=prod(range(1,stubs,2)) if stubs%2==0 else 0
    if total!=expected:raise ArithmeticError('Full Gaussian pairing sum does not equal (N-1)!!')
    return {'kinds':kinds,'labels':labels,'cases':cases,'counts':dict(counts),'full_wick_multiplicity':total,
            'gaussian_double_factorial':expected,'source_source_pairs_included':True,
            'pairings_explicitly_visited':0,'source_factorial':'fixed derivatives; no additional external E! division'}
