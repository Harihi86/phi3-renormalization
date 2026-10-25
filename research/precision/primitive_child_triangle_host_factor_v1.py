"""Primitive three-loop vertex in a one-loop triangle: exact forest reduction."""
from fractions import Fraction as Q
from itertools import combinations,permutations
import sympy as S
import primitive_triangle_factor_v3 as P


def basis(gd):
    selected=gd['proper_subgraphs'][0]['edge_ids'];outside=[e for e in range(12) if e not in selected]
    routes=[P.scaled_routes(term,[]) for term in gd['terms']]
    for own in combinations(selected,3):
        for other in outside:
            choices=(*own,other);witnesses=[]
            for rows in routes:
                change=S.Matrix.vstack(*(rows[e][:,:4] for e in choices));det=change.det()
                if det not in (-1,1):break
                offset=S.Matrix.vstack(*(rows[e][:,4:] for e in choices));inverse=change.inv();converted=[]
                for row in rows:
                    loop=row[:,:4]*inverse;external=row[:,4:]-loop*offset
                    converted.append(S.Matrix([[S.expand(v) for v in list(loop)+list(external)]]))
                if any(v for e in outside for v in converted[e][:,:3]):break
                if any(converted[e][3] not in (-1,1) for e in outside):break
                witnesses.append({'change':change,'offset':offset,'det':det,'routing':converted})
            if len(witnesses)==len(routes):return choices,witnesses
    raise AssertionError('No common three-child-loop unimodular basis')


def polynomial(rows,alpha):return S.expand(sum(co*S.prod(alpha[e] for e in range(len(alpha)) if mask&(1<<e)) for mask,co in rows))


def child_polynomials(routes):
    matrix=S.Matrix(routes);U=[];N=[[[] for _ in range(3)] for _ in range(3)]
    for subset in combinations(range(9),3):
        d=matrix.extract(subset,range(3)).det()
        if d:assert abs(d)==1;U.append([sum(1<<e for e in subset),1])
    for subset in combinations(range(9),4):
        ds=[matrix.extract(subset,[0,1,2,3+j]).det() for j in range(3)];mask=sum(1<<e for e in subset)
        for i in range(3):
            for j in range(3):
                if ds[i]*ds[j]:N[i][j].append([mask,int(ds[i]*ds[j])])
    assert U and N[0][0] and all(co>0 for _,co in N[0][0])
    return U,N


def grouped(routes,U):
    buckets={}
    for e,row in enumerate(routes):
        first=next(v for v in row[:3] if v);assert first in (-1,1)
        buckets.setdefault(tuple(first*v for v in row[:3]),[]).append(e)
    groups=list(buckets.values());assert sorted(map(len,groups))==[1,1,1,2,2,2]
    a=S.symbols('a:6');x=S.symbols('x:3');alpha=[None]*9;fractions=[]
    for g,edges in enumerate(groups):
        if len(edges)==1:alpha[edges[0]]=a[g]
        else:
            j=len(fractions);fractions.append({'group':g,'edges':edges,'index':j,'weighted':False})
            alpha[edges[0]]=a[g]*x[j];alpha[edges[1]]=a[g]*(1-x[j])
    ug=polynomial(U,alpha);assert not any(ug.has(v) for v in x)
    records=[]
    for powers,co in S.Poly(ug,*a).terms():
        assert co==1 and all(p in (0,1) for p in powers);records.append([sum(1<<j for j,p in enumerate(powers) if p),1])
    assert len(records)==16;degrees=list(map(len,groups));masks=[m for m,_ in records]
    nu={mask:min((mask&m).bit_count() for m in masks) for mask in range(64)}
    gamma={mask:sum(degrees[j] for j in range(6) if mask&(1<<j))-3*nu[mask] for mask in range(1,63)}
    assert min(gamma.values())>=1;H={0:Q(1)}
    for size in range(1,6):
        for mask in range(1,63):
            if mask.bit_count()==size:H[mask]=sum((H[mask^(1<<j)] for j in range(6) if mask&(1<<j)),Q())/gamma[mask]
    normalization=sum((H[63^(1<<j)] for j in range(6)),Q());literal=Q()
    for order in permutations(range(6)):
        mask=63;mass=Q(1);suffix=[]
        for e in order[:-1]:mask^=1<<e;mass/=gamma[mask];suffix.append(mask)
        assert any(all((m&s).bit_count()==nu[s] for s in suffix) for m in masks)
        literal+=mass
    assert literal==normalization
    return {'groups':groups,'fractions':fractions,'degrees':degrees,'U':records,'gamma':gamma,'nu':nu,'H':H,
            'normalization':normalization,'sector_orders':720,'all_nested_leading_monomials':True,
            'integrated_fraction_indices':[0,1,2],'exact_child_U_grouping':True}


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,4,3,0,12)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert len(gd['proper_subgraphs'])==1
    part=gd['proper_subgraphs'][0];assert (part['loops'],part['external_legs'],part['omega'],len(part['edge_ids']))==(3,3,0,9)
    assert len(gd['terms'])==len(pd['states'])==2 and {tuple(t['subgraphs_inner_first']) for t in gd['terms']}=={(),(0,)}
    assert {tuple(s['word']) for s in pd['states']}=={(),(0,)}
    assert raw['external_gram']==[[1.,-.5],[-.5,1.]] and raw['positive_invariants']==[1.,1.,1.]
    choices,witnesses=basis(gd);alpha=S.symbols('a:12');selected=part['edge_ids'];outside=[e for e in range(12) if e not in selected]
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]]);records=[];states={}
    for term,w in zip(gd['terms'],witnesses):
        word=tuple(term['subgraphs_inner_first']);ps=next(s for s in pd['states'] if tuple(s['word'])==word)
        assert term['sign']==ps['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==term['jet_orders']==ps['variables']==[]
        routes=w['routing'];child=sum((alpha[e]*routes[e].T*routes[e] for e in selected),S.zeros(6))
        host=sum((alpha[e]*routes[e][:,3:].T*routes[e][:,3:] for e in outside),S.zeros(3))
        A=child[:3,:3];B=child[:3,3:];Uc=S.expand(A.det());N=(Uc*child[3:,3:]-B.T*A.adjugate()*B).applyfunc(S.expand)
        effective=Uc*host+N;reduced_U=S.expand(effective[0,0])
        fn=S.expand(sum(gram[i,j]*effective.extract([0,1+i],[0,1+j]).det() for i in range(2) for j in range(2)))
        reduced_F=S.Poly(fn,*alpha).exquo(S.Poly(Uc,*alpha)).as_expr()
        original=sum((a*r.T*r for a,r in zip(alpha,P.scaled_routes(term,[]))),S.zeros(6))
        ug=S.expand(original[:4,:4].det(method='domain-ge'))
        fg=S.expand(sum(gram[i,j]*original.extract([0,1,2,3,4+i],[0,1,2,3,4+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        counts=[]
        for name,direct,reduced in (('U',ug,reduced_U),('F',fg,reduced_F)):
            positive=P.positive_polynomial(ps[name],ps,alpha,[])
            assert S.expand(direct-positive)==0 and S.expand(reduced-positive)==0,(raw['id'],word,name)
            counts.append(len(S.Poly(positive,*alpha).terms()))
        if word:assert N==S.zeros(3),'Child T0 left an external or host-loop component'
        states[word]=(host,Uc,N,routes)
        records.append({'word':list(word),'sign':term['sign'],'basis_edges':list(choices),'jacobian_determinant':int(w['det']),
            'loop_change':[[str(v) for v in w['change'].row(j)] for j in range(4)],
            'external_offset':[[str(v) for v in w['offset'].row(j)] for j in range(4)],
            'routing':[[str(v) for v in row] for row in routes],'U_F_monomials':counts,'three_exact_U_F_constructions_equal':True})
    host,Uc,N,routes=states[()];h2,u2,n2,projected=states[(0,)]
    assert host==h2 and Uc==u2 and N!=S.zeros(3)
    assert [routes[e][:,:3] for e in selected]==[projected[e][:,:3] for e in selected]
    assert all(v==0 for e in selected for v in projected[e][:,3:])
    child_routes=[[int(v) for v in routes[e]] for e in selected]
    host_routes=[[int(v/routes[e][3]) for v in routes[e][:,3:]] for e in outside]
    for i,j in combinations(range(3),2):
        delta=S.Matrix([host_routes[i][1:]])-S.Matrix([host_routes[j][1:]])
        assert (delta*gram*delta.T)[0]==1
    U,polys=child_polynomials(child_routes);betas=[alpha[e] for e in selected]
    assert S.expand(polynomial(U,betas)-Uc)==0
    assert all(S.expand(polynomial(polys[i][j],betas)-N[i,j])==0 for i in range(3) for j in range(3))
    return {'id':raw['id'],'states':records,'child_edges':selected,'host_edges':outside,'child_routes':child_routes,
            'host_routes':host_routes,'child_U':U,'child_variance_numerator':polys,'groups':grouped(child_routes,U),
            'child_loop_matrix_unchanged':True,'child_projected_variance_zero':True,'host_unchanged_under_child_T0':True,
            'all_three_host_pair_invariants_one':True,'complete_integral_dimension':11,
            'radial_density':'[log1p(Q0/(1+r))-(1+r*ell)^(-3)*log1p(Qr/(1+r))]/(r*Uc^3)',
            'radial_scale':'r=lambda*sum(child_a); W=N/(Uc*sum(child_a)); ell=W00',
            'whole_integral_evaluated':False}


def scalar_identities():
    r,m,q,d,U,V=S.symbols('r m q d U V',positive=True)
    direct=S.log(1+q/(m+r))/U**3-S.log(1+(q+d)/(m+r))/V**3
    stable=(1/U**3-1/V**3)*S.log(1+q/(m+r))-S.log(1+d/(m+r+q))/V**3
    assert S.simplify(S.diff(direct-stable,d))==0 and S.simplify((direct-stable).subs(d,0))==0
    # Nine child parameters and three host parameters supply rho^11 lambda^8;
    # the six-dimensional determinants supply rho^-12 lambda^-9.
    assert (12-1)-3*4==-1 and (9-1)-3*3==-1
    return {'complete_common_scale_log_identity':True,'common_scale_power':-1,'relative_scale_power':-1,
            'child_shape_measure_retained':True,'host_simplex_volume':'1/2','no_extra_factorial_or_symmetry_weight':True}
