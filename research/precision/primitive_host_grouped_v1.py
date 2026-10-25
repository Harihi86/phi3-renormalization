"""Merge equal loop coefficients before Hepp sampling; keep external fractions."""
from fractions import Fraction as Q
from itertools import combinations,permutations
import sympy as S
import primitive_host_bubble_factor_v1 as F


def build(binding):
    routes=binding['host_routes'];buckets={}
    for e,row in enumerate(routes):
        sign=next(v for v in row[:3] if v);assert sign in (-1,1)
        key=tuple(v*sign for v in row[:3]);buckets.setdefault(key,[]).append(e)
    groups=list(buckets.values());assert len(groups)==6 and sorted(map(len,groups))==[1,1,1,2,2,2]
    weighted=binding['weighted_host_edge']
    for edges in groups:
        if weighted in edges:edges.sort(key=lambda e:e!=weighted)
    grouped=[S.Matrix([routes[edges[0]][:3]]) for edges in groups];U=[]
    for subset in combinations(range(6),3):
        det=S.Matrix.vstack(*(grouped[g] for g in subset)).det()
        if det:assert abs(det)==1;U.append([sum(1<<g for g in subset),1])
    assert len(U)==16
    a=S.symbols('a:6',positive=True);x=S.symbols('x:3');alpha=[None]*9;fractions=[]
    for g,edges in enumerate(groups):
        if len(edges)==1:alpha[edges[0]]=a[g]
        else:
            j=len(fractions);fractions.append({'group':g,'edges':edges,'index':j,'weighted':weighted in edges})
            alpha[edges[0]]=a[g]*x[j];alpha[edges[1]]=a[g]*(1-x[j])
    assert S.expand(F.polynomial(binding['host_U'],alpha)-F.polynomial(U,a))==0
    fullF=F.polynomial(binding['host_F'],alpha);meanF=fullF;integrated=[]
    for row in fractions:
        if not row['weighted']:meanF=S.integrate(meanF,(x[row['index']],0,1));integrated.append(row['index'])
    # An independent Cauchy--Binet monomial average, using 1/2 and 1/6.
    avg=0
    for mask,co in binding['host_F']:
        term=S.Integer(co)
        for g,edges in enumerate(groups):
            count=sum(bool(mask&(1<<e)) for e in edges)
            if len(edges)==2 and weighted not in edges:
                term*=1 if count==0 else a[g]/2 if count==1 else a[g]**2/6
            else:term*=S.prod(alpha[e] for e in edges if mask&(1<<e))
        avg+=term
    assert S.expand(avg-meanF)==0
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    for row in fractions:
        e0,e1=row['edges'];r0=S.Matrix([routes[e0]]);r1=S.Matrix([routes[e1]])
        sign=next(v for v in routes[e0][:3] if v)*next(v for v in routes[e1][:3] if v)
        delta=r0[:,3:]-sign*r1[:,3:];assert (delta*gram*delta.T)[0]==1
    degrees=list(map(len,groups));masks=[m for m,_ in U];full=63
    nu={s:min((s&m).bit_count() for m in masks) for s in range(64)}
    gamma={s:sum(degrees[i] for i in range(6) if s&(1<<i))-3*nu[s] for s in range(1,full)}
    assert all(v>=1 for v in gamma.values())
    H={0:Q(1)}
    for size in range(1,6):
        for s in range(1,full):
            if s.bit_count()==size:H[s]=sum((H[s^(1<<i)] for i in range(6) if s&(1<<i)),Q())/gamma[s]
    total=sum((H[full^(1<<i)] for i in range(6)),Q());exact=Q()
    for order in permutations(range(6)):
        remaining=full;den=1
        for e in order[:-1]:remaining^=1<<e;den*=gamma[remaining]
        exact+=Q(1,den)
    assert total==exact
    return {'id':binding['id'],'groups':groups,'fractions':fractions,'degrees':degrees,'U':U,
        'gamma':gamma,'nu':nu,'H':H,'normalization':total,'sector_orders':720,
        'exact_host_U_grouping':True,'exact_unweighted_fraction_mean':True,'integrated_fraction_indices':integrated,
        'group_Jacobian':'product a_g for the three split-edge groups; each pair has absolute determinant a_g',
        'weighted_Jacobian':'ratio exponent is sum of group degrees in suffix minus1',
        'original_full_integrand_retained':True,'whole_integral_evaluated':False}
