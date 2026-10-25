"""Pair mixed E2 forests with per-forest unimodular bases and fixed edge identities."""
from itertools import combinations,product
import sympy as S
import primitive_triangle_factor_v3 as P
from nested_bubble_disjoint_triangle_factor_v1 import convert

def recognize(gd,z,bubble):
    # Identify the two propagators once in the bare graph. A parent's T0 can
    # degenerate one choice of host basis, so choose a unimodular basis per pair.
    import weighted_bubble_e2_factor_v2 as original
    bare=next(t for t in gd['terms'] if not t['subgraphs_inner_first'])
    single=dict(gd);single['terms']=[bare]
    be,host,duplicates,bare_choices,_=original.recognize(single,z,bubble)
    terms={frozenset(t['subgraphs_inner_first']):t for t in gd['terms']}
    by_word={}
    for word in sorted((w for w in terms if bubble not in w),key=lambda w:(len(w),sorted(w))):
        mate=word|{bubble};assert mate in terms
        source=[P.scaled_routes(terms[w],z) for w in (word,mate)]
        accepted=None
        for be0,kedge in product(be,duplicates):
            for others in combinations([e for e in host if e!=kedge],2):
                choices=(be0,kedge,*others);pair=[]
                for w,rows in zip((word,mate),source):
                    witness=convert(rows,choices)
                    if witness is None:break
                    rr=witness['routes'];scale=z[bubble] if bubble in w else S.Integer(1)
                    if any(rr[e][0] or any(v.has(*z) for v in rr[e]) for e in host):break
                    if any(rr[e][0] not in (-1,1) for e in be):break
                    delta=rr[be[0]][:,1:]/rr[be[0]][0]-rr[be[1]][:,1:]/rr[be[1]][0]
                    if delta[0] not in (-scale,scale) or any(delta[j] for j in range(1,4)):break
                    if any(abs(rr[e][1])!=1 or any(rr[e][j] for j in (0,2,3,4)) for e in duplicates):break
                    witness['basis_edges']=list(choices);pair.append(witness)
                if len(pair)==2 and all(pair[0]['routes'][e]==pair[1]['routes'][e] for e in host):
                    accepted=pair;break
            if accepted is not None:break
        assert accepted is not None,('No paired weighted-edge basis',gd['graph'],sorted(word))
        by_word[word],by_word[mate]=accepted
    assert len(by_word)==len(gd['terms'])
    return be,host,duplicates,bare_choices,[by_word[frozenset(t['subgraphs_inner_first'])] for t in gd['terms']]


def minors(routes,columns):
    matrix=S.Matrix(routes);terms=[]
    for ids in combinations(range(len(routes)),len(columns)):
        determinant=matrix.extract(ids,columns).det()
        assert determinant.is_Integer
        if determinant:terms.append([sum(1<<i for i in ids),int(determinant**2)])
    assert terms;return terms

def polynomial(terms,alpha):return S.expand(sum(co*S.prod(alpha[e] for e in range(len(alpha)) if mask&(1<<e)) for mask,co in terms))

def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,4,2,2,11)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert raw['external_gram']==[[1.0]] and raw['positive_invariants']==[1.0]
    parts=gd['proper_subgraphs'];bubbles=[i for i,p in enumerate(parts) if (p['loops'],p['external_legs'],p['omega'],len(p['edge_ids']))==(1,2,2,2)]
    assert len(bubbles)==1;bubble=bubbles[0];assert all(p['external_legs']==3 and p['omega']==0 for i,p in enumerate(parts) if i!=bubble)
    sets=[set(p['edge_ids']) for p in parts];words=set()
    for mask in range(1<<len(parts)):
        indices=[i for i in range(len(parts)) if mask&(1<<i)]
        if all(not sets[i]&sets[j] or sets[i]<sets[j] or sets[j]<sets[i] for i,j in combinations(indices,2)):words.add(frozenset(indices))
    assert len(words)==len(gd['terms'])==len(pd['states'])
    assert words=={frozenset(t['subgraphs_inner_first']) for t in gd['terms']}=={frozenset(t['word']) for t in pd['states']}
    z=S.symbols('z:'+str(len(parts)));alpha=S.symbols('a:11');be,host,duplicates,choices,witnesses=recognize(gd,z,bubble)
    lam=sum(alpha[e] for e in be);records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=term['subgraphs_inner_first'];ps=next(s for s in pd['states'] if frozenset(s['word'])==frozenset(word))
        assert term['sign']==ps['sign']==(-1)**len(word)
        active=[bubble] if bubble in word else [];assert term['jet_subgraphs']==ps['variables']==active and term['jet_orders']==[2]*len(active)
        for i,j in combinations(word,2):assert not sets[j]<sets[i]
        rows=w['routes'];matrix=sum((a*r.T*r for a,r in zip(alpha,rows)),S.zeros(5))
        H=sum((alpha[e]*rows[e][:,1:].T*rows[e][:,1:] for e in host),S.zeros(4));V=S.zeros(4)
        V[0,0]=alpha[be[0]]*alpha[be[1]]*(z[bubble]**2 if bubble in word else 1)
        reduced=lam*H+V
        assert matrix[0,0]==lam and (lam*matrix[1:,1:]-matrix[1:,:1]*matrix[:1,1:]-reduced).applyfunc(S.expand)==S.zeros(4)
        U=S.Poly(reduced[:3,:3].det(method='domain-ge'),*alpha,*z).exquo(S.Poly(lam**2,*alpha,*z)).as_expr()
        F=S.Poly(reduced.det(method='domain-ge'),*alpha,*z).exquo(S.Poly(lam**3,*alpha,*z)).as_expr()
        old=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,z))),S.zeros(5));counts=[]
        for name,direct,new in [('U',old[:4,:4].det(method='domain-ge'),U),('F',old.det(method='domain-ge'),F)]:
            positive=P.positive_polynomial(ps[name],ps,alpha,z)
            assert S.expand(direct-positive)==0 and S.expand(new-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha,*z).terms()))
        states[frozenset(word)]=(H,rows)
        records.append({'word':word,'sign':term['sign'],'jet_orders':term['jet_orders'],'basis_edges':w['basis_edges'],'jacobian_determinant':int(w['det']),
            'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(4)],'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(4)],
            'routing':[[str(v) for v in r] for r in rows],'U_F_monomials':counts,'three_exact_U_F_constructions_equal':True})
    baseedges=[e for e in host if e!=duplicates[1]];assert len(baseedges)==8;weighted=baseedges.index(duplicates[0]);beta=S.symbols('b:8');paired=[]
    for word in sorted((w for w in words if bubble not in w),key=lambda w:(len(w),sorted(w))):
        assert word|{bubble} in words
        H,rows=states[word];H2,rows2=states[word|{bubble}];assert H==H2 and all(rows[e]==rows2[e] for e in host)
        rr=[[int(v) for v in rows[e][:,1:]] for e in baseedges]
        loop_matrix=sum((a*S.Matrix([r]).T*S.Matrix([r]) for a,r in zip(beta,rr)),S.zeros(4))
        U,F=minors(rr,[0,1,2]),minors(rr,[0,1,2,3])
        assert polynomial(U,beta)==S.expand(loop_matrix[:3,:3].det(method='domain-ge'))
        assert polynomial(F,beta)==S.expand(loop_matrix.det(method='domain-ge'))
        paired.append({'word':sorted(word),'paired_bubble_word':sorted(word|{bubble}),'sign':(-1)**len(word),'routes':rr,'U':U,'F':F,
            'weighted_edge':weighted,'positive_minor_polynomials_equal_gaussian':True})
    assert len(paired)*2==len(records)
    return {'id':raw['id'],'states':records,'paired_states':paired,'bubble_index':bubble,'bubble_edges':be,'original_host_edges':host,
        'duplicated_host_edges':duplicates,'duplicate_identity_rule':'the two edges are identified in the bare graph and retained in every forest; extra equal routes after vertex projection remain separate propagators','effective_edges':baseedges,'weighted_edge':weighted,'effective_loops':3,'effective_edge_count':8,
        'exact_bubble_boundary':'first host loop k in each paired unimodular basis; original edge identities retained across every parent forest','per_forest_basis':True,
        'weighted_edge_factor':'B(k^2)/(1+k^2)^2','all_overlap_compatibility_and_local_orders_checked':True,
        'complete_kernel':'sum_f sign(f) J(F_f/U_f;M,d)/U_f^3',
        'M':'sum of eight effective Schwinger parameters','d':'parameter of the weighted edge','projective_simplex_volume':'1/7!',
        'all_four_loops_integrated_analytically':True,'candidate_dimension':7,'whole_graph_acceptance':False}
