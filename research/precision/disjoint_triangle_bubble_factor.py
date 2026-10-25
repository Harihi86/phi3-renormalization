"""Bind a disjoint primitive vertex and bubble to one remaining loop.

All four original forest states are retained. The two local loop integrations
are proved by exact Gaussian/positive-polynomial identities. Overall T0 acts
on physical external momenta at fixed remaining loop momentum.
"""
from itertools import product,combinations
import sympy as s
import primitive_triangle_factor_v3 as P


def convert(rows,choices):
    change=s.Matrix.vstack(*(rows[e][:,:3] for e in choices))
    det=s.factor(change.det())
    if det not in (-1,1):return None
    offset=s.Matrix.vstack(*(rows[e][:,3:] for e in choices))
    inv=change.inv();result=[]
    for row in rows:
        loop=row[:,:3]*inv;external=row[:,3:]-loop*offset
        result.append(s.Matrix([[s.expand(v) for v in list(loop)+list(external)]]))
    return {'det':det,'change':change,'offset':offset,'routes':result}


def recognize(gd,z):
    parts=gd['proper_subgraphs'];assert len(parts)==2
    tri=P.find_triangle(gd);bubble=1-tri
    tp,bp=parts[tri],parts[bubble];te=tp['edge_ids'];be=bp['edge_ids']
    assert (bp['loops'],bp['external_legs'],bp['omega'],len(be))==(1,2,2,2)
    assert tp['omega']==0 and not set(te)&set(be)
    edges=gd['graph']['edges'];assert set(edges[be[0]])==set(edges[be[1]]) and len(set(edges[be[0]]))==2
    host=sorted(set(range(gd['edges']))-set(te)-set(be));assert len(host)==4
    source=[P.scaled_routes(term,z) for term in gd['terms']]
    for choices in product(te,be,host):
        witnesses=[]
        for term,rows in zip(gd['terms'],source):
            w=convert(rows,choices)
            if w is None:break
            rr=w['routes'];word=term['subgraphs_inner_first']
            if any(rr[e][0] or rr[e][1] or rr[e][2] not in (-1,1) or any(v.has(*z) for v in rr[e]) for e in host):break
            if any(rr[e][0] not in (-1,1) or rr[e][1] or any(v.has(z[bubble]) for v in rr[e]) for e in te):break
            if any(rr[e][1] not in (-1,1) or rr[e][0] or rr[e][3] or rr[e][4] for e in be):break
            shifts=[rr[e][:,2:]/rr[e][1] for e in be]
            factor=z[bubble] if bubble in word else s.Integer(1)
            delta=shifts[1]-shifts[0]
            if delta[0] not in (-factor,factor) or delta[1] or delta[2]:break
            # Exactly two unshifted host propagators; preserve both.
            norms=[rr[e][:,2:]/rr[e][2] for e in host]
            if sum(row[1]==row[2]==0 for row in norms)!=2:break
            witnesses.append(w)
        if len(witnesses)==len(source):return tri,bubble,te,be,host,choices,witnesses
    raise AssertionError('No common original-edge basis isolates both child loops and their host')


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['overall_degree'],gd['edges'])==(6,3,3,0,9)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    parts=gd['proper_subgraphs'];z=s.symbols('z:'+str(len(parts)))
    tri,bub,te,be,he,choices,witnesses=recognize(gd,z)
    words={tuple(t['subgraphs_inner_first']) for t in gd['terms']}
    assert words=={(),(tri,),(bub,),tuple(sorted((tri,bub)))} and len(gd['terms'])==4
    assert {tuple(st['word']) for st in pd['states']}==words and len(pd['states'])==4
    alpha=s.symbols('a:9');lt=sum(alpha[e] for e in te);lb=sum(alpha[e] for e in be)
    gram=s.Matrix([[s.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram.shape==(2,2) and gram.det()>0
    assert all(v==1 for v in raw['positive_invariants'])
    states=[];matrices={};bare_shifts=None;host_shifts=None
    for term,w in zip(gd['terms'],witnesses):
        word=tuple(term['subgraphs_inner_first']);state=next(st for st in pd['states'] if tuple(st['word'])==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        active={i for i in word if parts[i]['omega']>0}
        assert set(term['jet_subgraphs'])==active==set(state['variables'])
        assert term['jet_orders']==[parts[i]['omega'] for i in term['jet_subgraphs']]
        rows=w['routes'];matrix=s.zeros(5);host=s.zeros(3)
        for e,row in enumerate(rows):
            matrix+=alpha[e]*row.T*row
            if e in he:host+=alpha[e]*row[:,2:].T*row[:,2:]
        assert matrix[:2,:2]==s.diag(lt,lb),'Local loops are coupled'
        vt=s.zeros(3);vb=s.zeros(3)
        ts=[rows[e][:,2:]/rows[e][0] for e in te]
        bs=[rows[e][:,2:]/rows[e][1] for e in be]
        for i,j in combinations(range(3),2):
            delta=ts[i]-ts[j];vt+=alpha[te[i]]*alpha[te[j]]*delta.T*delta
        delta=bs[0]-bs[1];vb=alpha[be[0]]*alpha[be[1]]*delta.T*delta
        vt=vt.applyfunc(s.expand);vb=vb.applyfunc(s.expand)
        schur=matrix[2:,2:]-matrix[2:,:2]*s.diag(1/lt,1/lb)*matrix[:2,2:]
        reduced=host+vt/lt+vb/lb
        assert (schur-reduced).applyfunc(s.cancel)==s.zeros(3)
        scale=z[bub]**2 if bub in word else 1
        expected=s.zeros(3);expected[0,0]=alpha[be[0]]*alpha[be[1]]*scale
        assert vb==expected,'Bubble boundary must be the retained loop q'
        if tri in word:assert vt==s.zeros(3)
        else:assert vt!=s.zeros(3) and not any(v.has(z[bub]) for v in vt)
        original=s.zeros(5)
        for a,row in zip(alpha,P.scaled_routes(term,z)):original+=a*row.T*row
        ug=s.expand(original[:3,:3].det(method='domain-ge'))
        fg=s.expand(sum(gram[i,j]*original.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        ur=s.cancel(lt*lb*reduced[0,0])
        fr=s.cancel(lt*lb*sum(gram[i,j]*reduced.extract([0,1+i],[0,1+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        counts=[]
        for name,original_poly,reduced_poly in (('U',ug,ur),('F',fg,fr)):
            pos=P.positive_polynomial(state[name],state,alpha,z)
            assert s.expand(original_poly-pos)==0 and s.expand(reduced_poly-pos)==0,(name,word)
            counts.append(len(s.Poly(pos,*alpha,*z).terms()))
        normalized=[rows[e][:,2:]/rows[e][2] for e in he]
        if host_shifts is None:host_shifts=normalized
        assert normalized==host_shifts
        if not word:bare_shifts=ts
        matrices[word]=(host,vt,vb)
        states.append({'word':list(word),'sign':term['sign'],'jet_orders':term['jet_orders'],
                       'basis_edges':list(choices),'jacobian_determinant':int(w['det']),
                       'loop_change':[[str(v) for v in w['change'].row(i)] for i in range(3)],
                       'external_offset':[[str(v) for v in w['offset'].row(i)] for i in range(3)],
                       'edge_routing':[[str(v) for v in row] for row in rows],
                       'U_F_monomials':counts,'three_exact_polynomial_constructions_equal':True})
    h0,t0,b0=matrices[()]
    for word,(h,t,b) in matrices.items():
        assert h==h0
        assert t==(s.zeros(3) if tri in word else t0)
        assert (b-(z[bub]**2 if bub in word else 1)*b0).applyfunc(s.expand)==s.zeros(3)
    assert bare_shifts is not None
    triangle_pairs=[]
    for i,j in combinations(range(3),2):
        triangle_pairs.append({'simplex_pair':[i,j],'momentum':[str(v) for v in bare_shifts[i]-bare_shifts[j]]})
    host_external=[[str(v) for v in row[1:]] for row in host_shifts if row[1] or row[2]]
    assert len(host_external)==2
    return {'triangle_index':tri,'bubble_index':bub,'triangle_edges':te,'bubble_edges':be,'host_edges':he,
            'states':states,'triangle_pairs':triangle_pairs,'host_external_shifts':host_external,
            'two_child_loops_eliminated':True,'outer_loop_preserved':True,'all_four_forests_retained':True,
            'whole_integral_not_computed':True,'normalization':gd['normalization']}


def finite_scale_identities():
    """Derivatives plus zero values bind the convergent Schwinger differences."""
    a,V=s.symbols('a V',positive=True)
    local=(1+V)*s.log(1+V)-V
    assert s.diff(local,V,2)==1/(1+V)
    assert local.subs(V,0)==0 and s.diff(local,V).subs(V,0)==0
    assert s.integrate(s.exp(-a*(1+V)),(a,0,s.oo))==1/(1+V)
    vertex=-s.log(1+V)
    assert s.diff(vertex,V)==-1/(1+V) and vertex.subs(V,0)==0
    assert s.integrate(-s.exp(-a*(1+V)),(a,0,s.oo))==-1/(1+V)
    x,t,z=s.symbols('x t z')
    assert s.expand(s.series(s.exp(-a*x*(1-x)*t*z*z),z,0,3).removeO()-(1-a*x*(1-x)*t*z*z))==0
    assert s.integrate((x*(1-x))**2/2,(x,0,1))==s.Rational(1,60)
    A,d,e=s.symbols('A d e')
    difference=-(A*(d+e)+d*e)/(A**4*(A+d)*(A+e))
    assert s.cancel(1/(A*A*(A+d)*(A+e))-1/A**4-difference)==0
    u,v=s.symbols('u v')
    assert s.cancel((1+u)*(1+(v-u)/(1+u))-(1+v))==0
    return {'bubble_subtraction':'integral_0^1 ((1+x*(1-x)*t)*log1p(x*(1-x)*t)-x*(1-x)*t) dx',
            'bubble_leading_coefficient':'1/60','triangle_subtraction':'-integral_simplex log1p(V)',
            'bubble_taylor_total_degree':2,'triangle_taylor_total_degree':0,'overall_taylor_total_degree':0,
            'stable_host_difference_exact':True,'log1p_product_identity_exact':True,
            'three_point_zero_sets_only_physical_external_momenta':True}


def kernel(ctx,binding,t,qp1,qp2,simplex,gram, *, bubble_value):
    """Stable overall-subtracted local kernel before the remaining measures.

    simplex is (x0,x1,x2) with sum one. q.p1 and q.p2 are not set to zero
    by the subgraph bubble subtraction. bubble_value must be B_R(q^2).
    The caller still owes the simplex and full six-dimensional q integrals.
    """
    def number(value):
        if isinstance(value,str) and '/' in value:
            num,den=value.split('/');return ctx.mpf(num)/ctx.mpf(den)
        return ctx.mpf(str(value)) if isinstance(value,s.Basic) else ctx.mpf(value)
    t,u,v=map(number,(t,qp1,qp2));x=list(map(number,simplex))
    assert t>=0 and len(x)==3 and all(y>=0 for y in x)
    assert abs(sum(x)-1)<=8*ctx.eps
    g11,g12,g22=map(number,(gram[0][0],gram[0][1],gram[1][1]))
    def norm_ext(a,b):return a*a*g11+2*a*b*g12+b*b*g22
    deltas=[]
    for a,b in binding['host_external_shifts']:
        a,b=number(a),number(b);deltas.append(2*(a*u+b*v)+norm_ext(a,b))
    A=1+t;d,e=deltas
    assert A+d>0 and A+e>0
    H=1/(A*A*(A+d)*(A+e));H0=1/A**4
    difference=-(A*(d+e)+d*e)/(A**4*(A+d)*(A+e))
    V0=ctx.mpf(0);DV=ctx.mpf(0)
    for row in binding['triangle_pairs']:
        i,j=row['simplex_pair'];c,a,b=map(number,row['momentum']);weight=x[i]*x[j]
        V0+=weight*c*c*t
        DV+=weight*(2*c*(a*u+b*v)+norm_ext(a,b))
    assert V0>=0 and 1+V0+DV>0
    logarithm_difference=ctx.log1p(DV/(1+V0))
    stable=-bubble_value*(difference*ctx.log1p(V0)+H*logarithm_difference)
    direct=-bubble_value*(H*ctx.log1p(V0+DV)-H0*ctx.log1p(V0))
    return {'stable':stable,'direct':direct,'H':H,'H0':H0,'H_difference':difference,
            'V0':V0,'V_difference':DV,'log_difference':logarithm_difference}
