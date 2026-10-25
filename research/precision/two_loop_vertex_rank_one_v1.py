"""L3E3G35: a primitive two-loop vertex inside a one-loop host.

Reuses only the proved scalar finite log pair. Original graph routings and
both U/F polynomials are bound afresh for this different topology.
"""
from fractions import Fraction as Q
import sympy as S
import triangle_vertex_rank_one_v1 as SECTORS
import primitive_triangle_factor_v3 as P
PERMUTATIONS=SECTORS.PERMUTATIONS
norm=SECTORS.norm


def verify(raw):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert raw['id']=='L3E3G35' and gd['loops']==3 and gd['edges']==9
    assert gd['dimension']==6 and gd['external_legs']==3 and gd['overall_degree']==0
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert tuple(map(Q,raw['positive_invariants']))==(Q(1),)*3
    gram=S.Matrix([[S.Rational(v) for v in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    parts=gd['proper_subgraphs'];assert len(parts)==1
    part=parts[0]
    assert (part['loops'],part['external_legs'],part['omega'])==(2,3,0)
    assert part['edge_ids']==list(range(6))
    assert len(gd['terms'])==len(pd['states'])==2
    alpha=S.symbols('a0:9',positive=True);w=S.symbols('w0:3',positive=True);f=S.symbols('f0:3',positive=True)
    sub={alpha[i]:w[i]*f[i] for i in range(3)}|{alpha[i+3]:w[i]*(1-f[i]) for i in range(3)}
    psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2]
    h=[w[i]*f[i]*(1-f[i]) for i in range(3)];kappa=S.prod(w)/psi
    x,y=f[0]-f[1],f[0]-f[2]
    a=h[0]+h[1]+kappa*x*x;b=h[0]+h[2]+kappa*y*y;c=h[0]+kappa*x*y
    determinant=h[0]*h[1]+h[0]*h[2]+h[1]*h[2]+kappa*(h[0]*(y-x)**2+h[1]*y*y+h[2]*x*x)
    assert S.cancel(a*b-c*c-determinant)==0
    expected=S.Matrix([[a,c,c],[c,b,b],[c,b,b]])
    records=[];hosts=[];words=[]
    for term in gd['terms']:
        word=term['subgraphs_inner_first'];words.append(word)
        state=next(r for r in pd['states'] if r['word']==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==state['variables']==term['jet_orders']==[]
        assert all(v['powers']==[] for row in term['scaled_routing'] for v in row)
        route=[S.Matrix([[S.Rational(v['coefficient']) for v in row]]) for row in term['scaled_routing']]
        original=sum((t*r.T*r for t,r in zip(alpha,route)),S.zeros(5))
        u=S.expand(original[:3,:3].det(method='domain-ge'))
        ff=S.expand(sum(gram[i,j]*original.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        assert S.expand(u-P.positive_polynomial(state['U'],state,alpha,[]))==0
        assert S.expand(ff-P.positive_polynomial(state['F'],state,alpha,[]))==0
        # q_original = q_host + p1 has unit Jacobian, for both forest states.
        converted=[S.Matrix([[r[0],r[1],r[2],r[3]+r[2],r[4]]]) for r in route]
        assert all(r[0]==r[1]==0 for r in converted[6:])
        child=sum((t*r.T*r for t,r in zip(alpha[:6],converted[:6])),S.zeros(5))
        host=sum((t*r[:,2:].T*r[:,2:] for t,r in zip(alpha[6:],converted[6:])),S.zeros(3))
        hosts.append(host)
        matrix=child[:2,:2];uc=S.factor(matrix.det())
        assert S.expand(uc.subs(sub)-psi)==0
        variance=child[2:,2:]-child[2:,:2]*matrix.inv()*child[:2,2:]
        desired=expected if word==[] else S.zeros(3)
        assert word in ([],[0])
        assert (variance.subs(sub)-desired).applyfunc(S.cancel)==S.zeros(3)
        records.append({'word':word,'sign':term['sign'],'original_U_F_exact':True,
            'child_loops_integrated':2,'child_Schur_form_exact':True,
            'unit_Jacobian_q_shift':[1,0],'converted_original_routing':[[str(v) for v in r] for r in converted]})
    assert sorted(words)==[[],[0]] and hosts[0]==hosts[1]
    h0,h1=S.symbols('h0 h1',positive=True);hh=[h0,h1,1-h0-h1]
    host=hosts[0].subs(dict(zip(alpha[6:],hh)));assert host[0,0]==1
    mu=-host[0,1:]
    assert mu==S.Matrix([[-hh[0]-hh[2],-hh[0]]])
    q0=S.expand(sum(gram[i,j]*(host[1:,1:]-host[1:,:1]*host[:1,1:])[i,j] for i in range(2) for j in range(2)))
    assert S.expand(q0-(hh[0]*hh[1]+hh[0]*hh[2]+hh[1]*hh[2]))==0
    # Abstract scalar Schur complement proves the same rank-one formula,
    # independently of the rational theta coefficient expressions above.
    aa,bb,cc,sigma,mm0,mm1,qq=S.symbols('a b c sigma m0 m1 Q',positive=True)
    pp=S.Matrix([[1,1]]);mu=S.Matrix([[mm0,mm1]])
    relative=1+sigma*aa;z=mu-sigma*cc*pp
    qnew=qq+norm(mm0,mm1)+sigma*bb-norm(z[0],z[1])/relative
    nn=aa*norm(mm0+cc/aa,mm1+cc/aa);dd=bb-cc*cc/aa
    assert S.cancel(qnew-qq-sigma*(nn/relative+dd))==0
    SECTORS.partition(PERMUTATIONS)
    return {'graph':'L3E3G35','states':records,'positive_child_determinant_exact':True,
        'host_external_momenta_retained':True,'child_T0_leaves_host_unchanged':True,
        'relative_scale_measure':'d sigma/sigma: sigma^5 Jacobian / sigma^6 determinant power',
        'six_child_Hepp_sectors':'each coefficient 1; v/(1+v+s*v)^3',
        'same_scalar_finite_log_pair_proved':True,'whole_integral_not_computed':True}


def chart(point,permutation):
    assert len(point)==8 and all(type(q) in (Q,int) and 0<q<1 for q in point)
    assert tuple(permutation) in PERMUTATIONS
    s,v,f0,f1,f2,r0,t,z0=point;den=1+s+s*v;weights=[None]*3
    for i,j in enumerate(permutation):weights[j]=(1/den,s/den,s*v/den)[i]
    w0,w1,w2=weights;psi=w0*w1+w0*w2+w1*w2
    ff=(f0,f1,f2);h=[w*f*(1-f) for w,f in zip(weights,ff)];kappa=w0*w1*w2/psi
    x,y=f0-f1,f0-f2
    a=h[0]+h[1]+kappa*x*x;c=h[0]+kappa*x*y
    det=h[0]*h[1]+h[0]*h[2]+h[1]*h[2]+kappa*(h[0]*(y-x)**2+h[1]*y*y+h[2]*x*x)
    d=det/a
    rd=r0*r0+(1-r0)**2;r=r0*r0/rd;dr=2*r0*(1-r0)/rd**2
    host=[r,(1-r)*t,(1-r)*(1-t)]
    q=host[0]*host[1]+host[0]*host[2]+host[1]*host[2]
    n=a*norm(c/a-host[0]-host[2],c/a-host[0]);u=a
    z=z0*z0;y0=1-z;cc=1+u;sigma=y0/(cc*z)
    d1,d2,d3=u+z,1+u*z,1+z*(u+q*cc)
    x0=q*cc*z/d2;x1=y0*(n*cc*z/d1+d)/d3
    ca=u*(d1*d1+d1*cc*z+(cc*z)**2)/d1**3*q*cc/d2
    cb=cc**3*z*z/d1**3*(n*cc*z/d1+d)/d3
    common=(1-r)*dr*2*z0
    assert min(psi,q,u,d,x0,x1)>0 and n>=0
    return {'Q':q,'u':u,'N':n,'d':d,'sigma':sigma,'child_U':psi,
        'alphas':[sigma*weights[i]*ff[i] for i in range(3)]+[sigma*weights[i]*(1-ff[i]) for i in range(3)]+host,
        'mass_sum':1+sigma,'stable_measure':v/(1+v+s*v)**3*common,
        'base_log_argument':x0,'difference_log_argument':x1,'positive_coefficient':ca,'negative_coefficient':cb,
        'original_projective_jacobian':w0*w1*w2*s/den**3*common*sigma**5/(cc*z*z)}


evaluate_chart=SECTORS.evaluate_chart
