"""Exact rank-one reorganization for the original L3E3G6 vertex insertion.

All Schwinger scales and six host sectors are retained. This module is a
candidate finite integrand, not an evaluated or certified whole integral.
"""
from fractions import Fraction as Q
from itertools import combinations,permutations
import sympy as S
import primitive_triangle_factor_v3 as ORIGINAL
PERMUTATIONS=tuple(permutations(range(3)))


def partition(values):
    assert all(len(p)==3 and all(type(v) is int for v in p) for p in values)
    assert len(values)==6 and len(set(map(tuple,values)))==6
    assert set(map(tuple,values))==set(PERMUTATIONS)


def norm(x,y):return x*x+y*y-x*y


def verify(raw):
    assert raw['id']=='L3E3G6'
    gd=raw['gaussian_full_descriptor']
    assert gd['loops']==3 and gd['edges']==9 and len(gd['proper_subgraphs'])==1
    gram=S.Matrix([[S.Rational(x) for x in row] for row in raw['external_gram']])
    assert gram==S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    previous=ORIGINAL.verify(raw);assert previous['triangle_original_edges']==[6,7,8]
    assert len(previous['states'])==2 and len(previous['paired_subtractions'])==1
    h=S.symbols('h0:6',positive=True);r,t,sigma=S.symbols('r t sigma',positive=True)
    x=[r,(1-r)*t,(1-r)*(1-t)]
    a=r*(1-r);c=a*t;b=x[1]*(1-x[1]);d=(1-r)*t*(1-t)
    assert S.expand(a*b-c*c-a*d)==0
    host=None;matrices={};bindings=[]
    for state in previous['states']:
        route=[S.Matrix([[S.Rational(v) for v in row]]) for row in state['original_edge_routing']]
        H=S.zeros(4)
        for e in range(6):
            assert route[e][0]==0
            H+=h[e]*route[e][:,1:].T*route[e][:,1:]
        if host is None:host=H
        else:assert H==host
        V=S.zeros(4)
        shifts=[route[e][:,1:]/route[e][0] for e in (6,7,8)]
        for i,j in combinations(range(3),2):
            delta=shifts[i]-shifts[j];V+=x[i]*x[j]*delta.T*delta
        expected=S.zeros(4)
        if state['word']==[]:
            expected[0,0]=a;expected[0,2]=expected[2,0]=c;expected[2,2]=b
        else:assert state['word']==[0]
        assert (V-expected).applyfunc(S.expand)==S.zeros(4)
        full=S.zeros(5)
        for alpha,row in zip([*h,*[sigma*z for z in x]],route):full+=alpha*row.T*row
        uf=S.expand(full[:3,:3].det(method='domain-ge'))
        ff=S.expand(sum(gram[i,j]*full.extract([0,1,2,3+i],[0,1,2,3+j]).det(method='domain-ge') for i in range(2) for j in range(2)))
        matrices[tuple(state['word'])]=(uf,ff)
        bindings.append({'word':state['word'],'sign':state['sign'],'original_edge_routing':state['original_edge_routing'],
                         'jacobian_determinant':state['jacobian_determinant'],'triangle_variance_exact':True})
    A=host[:2,:2];B=host[:2,2:];inv=A.inv();U=S.factor(A.det())
    q0=S.factor(sum(gram[i,j]*(host[2:,2:]-B.T*inv*B)[i,j] for i in range(2) for j in range(2)))
    mu=(-inv*B)[0,:];variance=S.factor(inv[0,0]);u=a*variance
    shifted_norm=S.factor(norm(mu[0]+t,mu[1]));M=a*shifted_norm+d;K=u*d;R=1+sigma*u
    assert S.cancel(matrices[()][0]-sigma*U*R)==0
    assert S.cancel(matrices[()][1]-sigma*U*(R*q0+sigma*(M+sigma*K)))==0
    assert S.cancel(matrices[(0,)][0]-sigma*U)==0
    assert S.cancel(matrices[(0,)][1]-sigma*U*q0)==0
    w=S.symbols('w0:3',positive=True);f=S.symbols('f0:3',positive=True)
    sub={h[i]:w[i]*(1-f[i]) for i in range(3)}|{h[i+3]:w[i]*f[i] for i in range(3)}
    U0=w[0]*w[1]+w[0]*w[2]+w[1]*w[2]
    P=(f[0]-f[2],f[1]-f[2]);positive_q=sum(w[i]*f[i]*(1-f[i]) for i in range(3))+S.prod(w)*norm(*P)/U0
    positive_mu=(w[1]*w[2]*P[0]/U0-f[0],w[1]*w[2]*P[1]/U0)
    assert S.cancel(q0.subs(sub)-positive_q)==0
    assert all(S.cancel(mu[i].subs(sub)-positive_mu[i])==0 for i in range(2))
    assert S.cancel(variance.subs(sub)-(w[1]+w[2])/U0)==0
    # The five host simplex coordinates factor into two path weights and three
    # edge fractions. The six sectors permute parameter weights, not graph legs.
    W0,W1=S.symbols('W0 W1',positive=True);W=[W0,W1,1-W0-W1]
    host_parameters=[W[0]*(1-f[0]),W[1]*(1-f[1]),W[2]*(1-f[2]),W[0]*f[0],W[1]*f[1]]
    host_jac=S.factor(S.Matrix(host_parameters).jacobian([W0,W1,*f]).det())
    assert S.cancel(host_jac**2-S.prod(W)**2)==0
    s,v=S.symbols('s v',positive=True);den=1+s+s*v;ordered=[1/den,s/den,s*v/den]
    witnesses=[]
    for perm in PERMUTATIONS:
        weights=[None]*3
        for i,j in enumerate(perm):weights[j]=ordered[i]
        jac=S.factor(S.Matrix(weights[:2]).jacobian([s,v]).det())
        assert S.cancel(jac**2-(s/den**3)**2)==0
        psi=weights[0]*weights[1]+weights[0]*weights[2]+weights[1]*weights[2]
        density=S.factor(S.prod(weights)*s/den**3/psi**3)
        assert S.cancel(density-v/(1+v+s*v)**3)==0
        witnesses.append({'parameter_order':perm,'weight':'v/(1+v+s*v)^3','multiplicity':1})
    # Check the rational coefficients of both log1p(x)/x factors after the
    # scale map. No logarithms are expanded or truncated in the integrand.
    uu,qq,nn,dd,y=S.symbols('u Q N d y',positive=True);z=1-y;cc=1+uu
    sig=y/(cc*z);rel=1+sig*uu;mass=1+sig
    base=qq/mass;diff=sig*(nn/rel+dd)/(mass+qq)
    left_a=(1-rel**-3)*base/(y*z);left_b=rel**-3*diff/(y*z)
    den1=uu+z;den2=1+uu*z;den3=1+z*(uu+qq*cc)
    right_a=uu*(den1**2+den1*cc*z+(cc*z)**2)/den1**3 *qq*cc/den2
    right_b=cc**3*z**2/den1**3*(nn*cc*z/den1+dd)/den3
    assert S.cancel(left_a-right_a)==S.cancel(left_b-right_b)==0
    q=S.symbols('q',positive=True);rr=q*q/(q*q+(1-q)**2)
    assert S.cancel(S.diff(rr,q)-2*q*(1-q)/(q*q+(1-q)**2)**2)==0
    partition(PERMUTATIONS)
    return {'graph':'L3E3G6','states':bindings,'rank_one_loop_update_exact':True,
        'original_U_F_rebuilt_after_triangle_scaling':True,'positive_host_formula_exact':True,
        'scale_log_coefficients_exact':True,'host_simplex_jacobian_exact':True,
        'sectors':witnesses,'triangle_simplex_jacobian':'1-r',
        'overall_scale_pair':'-log(1+Q/(1+sigma))',
        'relative_triangle_scale_measure':'d sigma/sigma',
        'external_and_host_momenta_not_projected_with_triangle':True,
        'no_series_truncation':True,'whole_graph_certificate':False}


def chart(point,permutation):
    assert all(type(q) in (Q,int) for q in point)
    assert len(point)==8 and all(Q(0)<q<Q(1) for q in point)
    assert tuple(permutation) in PERMUTATIONS
    s,v,f0,f1,f2,r0,t,z0=point;den=1+s+s*v;weights=[None]*3
    for i,j in enumerate(permutation):weights[j]=(1/den,s/den,s*v/den)[i]
    A,B,C=weights;U=A*B+A*C+B*C;f=(f0,f1,f2)
    p0,p1=f0-f2,f1-f2
    q0=sum(w*h*(1-h) for w,h in zip(weights,f))+A*B*C*norm(p0,p1)/U
    mu0=B*C*p0/U-f0;mu1=B*C*p1/U
    rden=r0*r0+(1-r0)**2;r=r0*r0/rden;dr=2*r0*(1-r0)/rden**2
    a=r*(1-r);d=(1-r)*t*(1-t);N=a*norm(mu0+t,mu1)
    u=a*(B+C)/U;z=z0*z0;y=1-z;cc=1+u
    sigma=y/(cc*z);D1=u+z;D2=1+u*z;D3=1+z*(u+q0*cc)
    x0=q0*cc*z/D2;x1=y*(N*cc*z/D1+d)/D3
    ca=u*(D1*D1+D1*cc*z+(cc*z)**2)/D1**3*q0*cc/D2
    cb=cc**3*z*z/D1**3*(N*cc*z/D1+d)/D3
    assert min(U,q0,u,D1,D2,D3)>0 and min(N,d,x0,x1,ca,cb)>=0
    host_jac=A*B*C*s/den**3
    common=(1-r)*dr*2*z0
    return {'host_weights':weights,'edge_fractions':f,'triangle_weights':[r,(1-r)*t,(1-r)*(1-t)],
        'sigma':sigma,'alphas':[weights[i]*(1-f[i]) for i in range(3)]+[weights[i]*f[i] for i in range(3)]+
                 [sigma*r,sigma*(1-r)*t,sigma*(1-r)*(1-t)],
        'base_log_argument':x0,'difference_log_argument':x1,'positive_coefficient':ca,'negative_coefficient':cb,
        'stable_measure':v/(1+v+s*v)**3*common,
        'original_projective_jacobian':host_jac*common*sigma*sigma/(cc*z*z),
        'mass_sum':1+sigma,'host_U':U,'host_Q':q0}


def evaluate_chart(record,mp):
    def number(q):return mp.mpf(q.numerator)/q.denominator
    def l(x):return mp.log1p(x)/x if x else mp.mpf(1)
    return number(record['stable_measure'])*(number(record['positive_coefficient'])*l(number(record['base_log_argument']))
        -number(record['negative_coefficient'])*l(number(record['difference_log_argument'])))
