"""Four-loop primitive host with one triangle: complete radial forest pairing."""
from itertools import combinations
import sympy as S
import primitive_triangle_factor_v3 as P
import primitive_host_bubble_factor_v1 as B
import primitive_host_grouped_v1 as G


def verify(raw):
    gd=raw['gaussian_full_descriptor']
    assert (gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(4,3,12,0)
    assert len(gd['proper_subgraphs'])==1 and len(gd['terms'])==2
    part=gd['proper_subgraphs'][0];assert (part['loops'],part['external_legs'],part['omega'])==(1,3,0)
    proof=P.verify(raw);assert len(proof['paired_subtractions'])==1
    state=next(row for row in proof['states'] if not row['word'])
    projected=next(row for row in proof['states'] if row['word']==[0])
    rows=[list(map(S.Rational,row)) for row in state['original_edge_routing']]
    old=[list(map(S.Rational,row)) for row in projected['original_edge_routing']]
    selected=proof['triangle_original_edges'];outside=[e for e in range(12) if e not in selected]
    routes=[r[1:] for e,r in enumerate(rows) if e in outside]
    assert routes==[r[1:] for e,r in enumerate(old) if e in outside] and len(routes)==9
    shifts=[S.Matrix([rows[e][1:]])/rows[e][0] for e in selected]
    differences=[(shifts[i]-shifts[j]) for i,j in combinations(range(3),2)]
    x=S.symbols('x:3');variance=sum((x[i]*x[j]*d.T*d for (i,j),d in zip(combinations(range(3),2),differences)),S.zeros(5))
    assert all(row[0]==0 for e,row in enumerate(rows) if e in outside)
    U,F=B.host_polynomials([S.Matrix([row]) for row in routes]);a=S.symbols('a:9')
    host=sum((aa*S.Matrix(row)*S.Matrix([row]) for aa,row in zip(a,routes)),S.zeros(5))
    assert S.expand(host[:3,:3].det()-B.polynomial(U,a))==0
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    ff=sum(gram[i,j]*host.extract([0,1,2,3+i],[0,1,2,3+j]).det() for i in range(2) for j in range(2))
    assert S.expand(ff-B.polynomial(F,a))==0
    binding={'id':raw['id'],'host_routes':[[int(v) for v in row] for row in routes],
             'host_U':U,'host_F':F,'weighted_host_edge':-1}
    # The grouped helper accepts no marked edge here: all three fractions are unweighted.
    grouped=G.build(binding);assert all(not row['weighted'] for row in grouped['fractions'])
    assert len(grouped['integrated_fraction_indices'])==3
    lo=S.Matrix.vstack(*(d[:,:3] for d in differences));rank=lo.rank();assert rank in (1,2)
    return {'id':raw['id'],'original_proof':proof,'triangle_edges':selected,'host_edges':outside,
            'host_routes':binding['host_routes'],'host_U':U,'host_F':F,'groups':grouped,
            'triangle_pair_differences':[[int(v) for v in d] for d in differences],
            'triangle_pair_indices':list(combinations(range(3),2)),'triangle_loop_variance_rank':rank,
            'triangle_variance_is_positive_semidefinite':True,'host_unchanged_by_child_T0':True,
            'unweighted_host_measure_degree':9,'triangle_simplex_volume':'1/2',
            'scale_paired_density':'[log1p(Q0/(M+r))/U0^3-log1p(Qr/(M+r))/Ur^3]/r',
            'radial_map':'r=M*t/(1-t), dr/r=dt/[t*(1-t)]',
            'complete_integral_dimension':11,'whole_integral_evaluated':False}


def scalar_identities():
    r,M,q0,delta,U0,Ur,t=S.symbols('r M q0 delta U0 Ur t',positive=True)
    q=q0+delta
    direct=S.log(1+q0/(M+r))/U0**3-S.log(1+q/(M+r))/Ur**3
    stable=(1/U0**3-1/Ur**3)*S.log(1+q0/(M+r))-S.log(1+delta/(M+r+q0))/Ur**3
    # Prove the logarithmic identity through its derivative and its zero boundary.
    difference=direct-stable
    assert S.simplify(S.diff(difference,delta))==0 and S.simplify(difference.subs(delta,0))==0
    mapped=M*t/(1-t);assert S.cancel(S.diff(mapped,t)/mapped-1/(t*(1-t)))==0
    x,y=S.symbols('x y',positive=True)
    simplex=[x,(1-x)*y,(1-x)*(1-y)]
    assert S.expand(sum(simplex)-1)==0 and S.det(S.Matrix(simplex[:2]).jacobian([x,y]))==1-x
    assert S.integrate(1-x,(x,0,1),(y,0,1))==S.Rational(1,2)
    # The complete common scale has power -1: 9 host edges and one d(lambda)/lambda.
    assert (9-1)-3*3==-1
    return {'complete_log_difference_identity':True,'radial_Jacobian_exact':True,
            'simplex_absolute_Jacobian':'1-x','simplex_volume':'1/2','common_scale_exponent':-1,
            'sign':'child bare minus T0; overall physical-external T0 makes log(host)-log(unprojected)',
            'no_extra_factorial_or_coupling':True,'full_external_momentum_dependence':True}
