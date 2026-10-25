"""Positive effective minors and independent conditional moments for all corners."""
from fractions import Fraction as Q
from itertools import combinations,permutations,product
from math import factorial
import sympy as S


def polynomial(rows,alpha):return S.expand(sum(co*S.prod(alpha[e] for e in range(len(alpha)) if mask&(1<<e)) for mask,co in rows))


def prepare(binding):
    L=binding['host_loops'];host=binding['host_routes'];effective=host+[r for rows in binding['child_pair_differences'] for r in rows]
    assert len(effective)==12;matrix=S.Matrix(effective);U=[];F=[]
    for subset in combinations(range(12),L):
        co=matrix.extract(subset,range(L)).det()**2
        if co:U.append([sum(1<<e for e in subset),int(co)])
    for subset in combinations(range(12),L+1):
        ds=[matrix.extract(subset,[*range(L),L+j]).det() for j in range(2)]
        co=ds[0]**2+ds[1]**2-ds[0]*ds[1];assert co>=0 and co.is_Integer
        if co:F.append([sum(1<<e for e in subset),int(co)])
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    buckets={}
    for e,row in enumerate(host):
        sign=next(v for v in row[:L] if v);assert sign in (-1,1)
        buckets.setdefault(tuple(sign*v for v in row[:L]),[]).append(e)
    groups=list(buckets.values())
    if L==2:
        assert sorted(map(len,groups))==[2,2,2];a=S.symbols('a:3');x=S.symbols('x:3');alphas=[None]*6
        for g,(i,j) in enumerate(groups):alphas[i]=a[g]*x[g];alphas[j]=a[g]*(1-x[g])
        ug=polynomial([r for r in U if r[0]<(1<<6)],alphas)
        assert S.expand(ug-(a[0]*a[1]+a[0]*a[2]+a[1]*a[2]))==0
        for i,j in groups:
            si=next(v for v in host[i][:L] if v);sj=next(v for v in host[j][:L] if v)
            diff=S.Matrix([host[i][L:]])-si*sj*S.Matrix([host[j][L:]])
            assert (diff*gram*diff.T)[0]==1
        gamma={1:2,2:2,4:2,3:1,5:1,6:1};nu={0:0,1:0,2:0,4:0,3:1,5:1,6:1,7:2}
        total=Q()
        for order in permutations(range(3)):
            mask=7;mass=Q(1)
            for e in order[:-1]:mask^=1<<e;mass/=gamma[mask]
            total+=mass
        assert total==3;measure={'kind':'three-paired-host-groups','groups':groups,'gamma':gamma,'nu':nu,'normalization':total,'orders':6}
    else:
        assert L==1 and len(host)==3 and len(groups)==1
        normalized=[S.Matrix([row[1:]])/row[0] for row in host]
        for i,j in combinations(range(3),2):
            d=normalized[i]-normalized[j];assert (d*gram*d.T)[0]==1
        measure={'kind':'host-triangle-simplex','groups':groups,'normalization':Q(1,2),'orders':1}
    return {'id':binding['id'],'host_loops':L,'host_edges':len(host),'effective_routes':effective,
            'effective_U':U,'effective_F':F,'host_measure':measure,'all_coefficients_nonnegative':True}


def verify_point(binding,positive,number,active):
    L=binding['host_loops'];n=L+2;host=binding['host_routes'];count=len(host);groups=positive['host_measure']['groups']
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    if L==2:
        scales=[S.Rational(j+1,number+j+2) for j in range(3)];x=[S.Rational(j+1,number+j+5) for j in range(3)]
        def alpha_at(xs):
            result=[None]*count
            for g,(i,j) in enumerate(groups):result[i]=scales[g]*xs[g];result[j]=scales[g]*(1-xs[g])
            return result
        current=alpha_at(x);center=alpha_at([S.Rational(1,2)]*3)
    else:
        current=[S.Rational(number+1,2*number+7),S.Rational(2,2*number+7),S.Rational(number+4,2*number+7)]
        center=[S.Rational(1,3)]*3
    extra=[]
    for i in range(len(binding['children'])):
        xs=[S.Rational(number+1,2*number+7+i),S.Rational(i+2,2*number+7+i),S.Rational(number+4,2*number+7+i)]
        assert sum(xs)==1;r=S.Rational(number+i+1,i+3) if i in active else S.Integer(0)
        extra.extend(r*xs[j]*xs[k] for j,k in combinations(range(3),2))
    routes=[S.Matrix([row]) for row in positive['effective_routes']]
    def gaussian(alpha):
        matrix=sum((v*r.T*r for v,r in zip(alpha+extra,routes)),S.zeros(n));A=matrix[:L,:L];inverse=A.inv()
        q=S.trace(gram*(matrix[L:,L:]-matrix[L:,:L]*inverse*matrix[:L,L:]))
        return A.det(),S.cancel(q),inverse
    U,Qpoint,inverse=gaussian(current);u=polynomial(positive['effective_U'],current+extra);f=polynomial(positive['effective_F'],current+extra)
    assert U==u and S.cancel(f/u-Qpoint)==0
    _,centerq,_=gaussian(center)
    if L==2:
        mean=centerq-sum(scales[g]**2*(routes[edges[0]][:,:L]*inverse*routes[edges[0]][:,:L].T)[0]/12 for g,edges in enumerate(groups))
        nodes=[(S.Integer(0),S.Rational(1,6)),(S.Rational(1,2),S.Rational(2,3)),(S.Integer(1),S.Rational(1,6))]
        quadrature=sum(S.prod(p[1] for p in row)*gaussian(alpha_at([p[0] for p in row]))[1] for row in product(nodes,repeat=3))
    else:
        mean=centerq-S.Rational(1,12)/U
        quadrature=sum(gaussian(list(p))[1] for p in set(permutations([S.Rational(2,3),S.Rational(1,6),S.Rational(1,6)])))/3
    avgF=0
    for mask,co in positive['effective_F']:
        term=S.Integer(co)
        if L==2:
            for g,edges in enumerate(groups):
                degree=sum(bool(mask&(1<<e)) for e in edges)
                term*=1 if degree==0 else scales[g]/2 if degree==1 else scales[g]**2/6
        else:
            degree=(mask&7).bit_count();term*=S.Rational(2,factorial(degree+2))
        term*=S.prod(extra[e-count] for e in range(count,12) if mask&(1<<e));avgF+=term
    assert S.cancel(mean-quadrature)==0 and S.cancel(avgF/U-mean)==0
    return {'active':list(active),'point':number,'U':Q(U),'Q':Q(Qpoint),'conditional_mean_Q':Q(mean),
            'positive_and_original_equal':True,'three_conditional_mean_constructions_equal':True}
