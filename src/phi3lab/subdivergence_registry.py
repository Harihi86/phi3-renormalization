"""Structural coverage and exact unit-point external-permutation certificates.

The reduction concerns unweighted amplitudes with equal positive masses and
unit symmetric invariants only. It does not identify asymmetric amplitudes.
Every member retains an explicit vertex, edge and external-label bijection.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations
from .automatic_bph import compile_bph
from .higher_graphs import canonical, generate


def classify_parts(parts):
    edges=[set(p['edge_ids']) for p in parts]
    vertices=[set(p['vertices']) for p in parts]
    nested=[];overlap=[];disjoint=[]
    for a,b in combinations(range(len(parts)),2):
        A,B=edges[a],edges[b]
        if A<B:nested.append([a,b])
        elif B<A:nested.append([b,a])
        elif vertices[a].isdisjoint(vertices[b]):disjoint.append([a,b])
        else:overlap.append([a,b])
    def depth(i):return 1+max((depth(j) for j,A in enumerate(edges) if A<edges[i]),default=0)
    if nested and overlap:category='mixed'
    elif overlap:category='overlap'
    elif nested:category='nested'
    elif disjoint:category='disjoint'
    elif parts:category='single'
    else:category='primitive'
    return {'category':category,'nested_pairs':nested,'overlap_pairs':overlap,
            'disjoint_pairs':disjoint,'has_nested':bool(nested),'has_overlap':bool(overlap),
            'maximum_nesting_depth':max((depth(i) for i in range(len(parts))),default=0),
            'proper_uv_subgraphs':len(parts),'quadratic_parts':sum(p['omega']==2 for p in parts),
            'logarithmic_parts':sum(p['omega']==0 for p in parts)}


def _adjacency(g):
    a=[[0]*g.n for _ in range(g.n)]
    for v,w in g.edges:
        a[v][w]+=1
        if v!=w:a[w][v]+=1
    return a


def external_unlabelled_isomorphism(source,target):
    """Independent exact backtracking witness, checked against all multiedges."""
    if source.n!=target.n or len(source.edges)!=len(target.edges) or len(source.external)!=len(target.external):
        raise ValueError('Graphs have different discrete sizes')
    A,B=_adjacency(source),_adjacency(target)
    pa,pb=Counter(v for _,v in source.external),Counter(v for _,v in target.external)
    colors=lambda g,m,p:[(g.vertex_kinds[v],p[v],m[v][v],tuple(sorted(m[v]))) for v in range(g.n)]
    ca,cb=colors(source,A,pa),colors(target,B,pb)
    candidates={v:[w for w in range(target.n) if ca[v]==cb[w]] for v in range(source.n)}
    mapping={};used=set()
    def search():
        if len(mapping)==source.n:return True
        options=[]
        for v in range(source.n):
            if v in mapping:continue
            viable=[w for w in candidates[v] if w not in used and all(A[v][q]==B[w][r] for q,r in mapping.items())]
            options.append((len(viable),-sum(A[v][q] for q in mapping),v,viable))
        _,_,v,viable=min(options)
        for w in viable:
            mapping[v]=w;used.add(w)
            if search():return True
            used.remove(w);del mapping[v]
        return False
    if not search():raise ValueError('No external-unlabelled graph isomorphism')
    vertex_map=[mapping[v] for v in range(source.n)]
    queues=defaultdict(deque)
    for e,pair in enumerate(target.edges):queues[pair].append(e)
    edge_map=[]
    for v,w in source.edges:
        pair=tuple(sorted((mapping[v],mapping[w])))
        if not queues[pair]:raise ArithmeticError('Multiedge bijection failed')
        edge_map.append(queues[pair].popleft())
    if any(queues.values()) or len(set(edge_map))!=len(source.edges):raise ArithmeticError('Incomplete edge bijection')
    target_ports={v:label for label,v in target.external}
    if len(target_ports)!=len(target.external):raise ValueError('This certificate requires distinct external attachment vertices')
    external_map={label:target_ports[mapping[v]] for label,v in source.external}
    E=len(source.external)
    if E not in (2,3):raise ValueError('Numerical orbit certificates use two or three external legs')
    gram=[[Fraction(1) if i==j else Fraction(-1,E-1) for j in range(E)] for i in range(E)]
    labels=[label for label,_ in source.external];target_labels=[label for label,_ in target.external]
    perm=[target_labels.index(external_map[label]) for label in labels]
    assert all(gram[i][j]==gram[perm[i]][perm[j]] for i in range(E) for j in range(E))
    return {'vertex_map':vertex_map,'edge_map':edge_map,'external_label_map':{str(k):v for k,v in external_map.items()},
            'unit_full_external_gram':[[str(x) for x in row] for row in gram],
            'edge_parameter_jacobian':1,'graph_bijection_verified':True,'external_gram_invariant_verified':True,
            'scope':'Equal common mass; unit symmetric Euclidean invariants; unweighted finite full zero-momentum recursive BPH amplitude. Graph/parameter permutation preserves the common original-edge regulator and conjugates every total Taylor subtraction. Not valid for generic unequal invariants.'}


def build_registry():
    rows=[];orbits={};counts=Counter();representatives=[]
    for loops in (3,4):
      for legs in (1,2,3):
        for index,g in enumerate(generate(loops,legs),1):
            ident=f'L{loops}E{legs}G{index}'
            c=compile_bph(g);structure=classify_parts(c.parts)
            row={'id':ident,'loops':loops,'legs':legs,'index':index,'graph':g.to_dict(),
                 **structure,'recursive_terms':len(c.terms),'exact_zero':legs==1}
            counts[(loops,legs,structure['category'])]+=1
            if legs==1:
                row.update(representative=ident,representation='exact_one_point_Taylor_identity')
            else:
                key=(loops,legs,canonical(g,fixed=False)[0])
                if key not in orbits:
                    orbits[key]=(ident,g,structure['category']);representatives.append(ident)
                    row.update(representative=ident,representation='must_execute')
                else:
                    rep,rg,category=orbits[key]
                    assert category==structure['category']
                    row.update(representative=rep,representation='exact_symmetric_point_isomorphism',
                               isomorphism=external_unlabelled_isomorphism(g,rg))
            rows.append(row)
    return {'schema_version':1,'loops':[3,4],'mass_squared':1,'symmetric_invariants':{'E2':[1],'E3':[1,1,1]},
            'total_graphs':len(rows),'nonzero_graphs':sum(not r['exact_zero'] for r in rows),
            'exact_zero_graphs':sum(r['exact_zero'] for r in rows),
            'independent_representatives':len(representatives),'representatives':representatives,
            'counts':[{'loops':k[0],'legs':k[1],'category':k[2],'graphs':v} for k,v in sorted(counts.items())],
            'graphs':rows,'numerical_integrations_executed_by_registry':0,
            'scope':'Structural census and exact discrete symmetry witnesses only; it does not claim any numerical integral has been executed.'}
