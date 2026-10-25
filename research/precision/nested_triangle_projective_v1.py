"""Original-edge projective forests with exact paired logarithm subtraction.

The order proposal is a change of variables, not a convergence certificate.
All twelve original parameters and all recursively generated forests remain.
"""
from fractions import Fraction as Q
from itertools import combinations
import math
import exact_rational_sector_transport_v1 as H


def determinant(rows):
    """Fraction-free integer determinant, with explicit exact divisions."""
    a=[list(row) for row in rows];n=len(a)
    assert all(len(row)==n and all(type(v) is int for v in row) for row in a)
    sign=previous=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        value=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=value*a[i][j]-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
        for i in range(k+1,n):a[i][k]=0
        previous=value
    return sign*a[-1][-1]


def minors(routes):
    assert len(routes)==12 and all(len(row)==6 for row in routes)
    u=[];f=[]
    for selected in combinations(range(12),4):
        d=determinant([[routes[e][j] for j in range(4)] for e in selected])
        if d:u.append([sum(1<<e for e in selected),d*d])
    for selected in combinations(range(12),5):
        ds=[determinant([[routes[e][j] for j in [0,1,2,3,4+column]] for e in selected]) for column in range(2)]
        coefficient=ds[0]**2+ds[1]**2-ds[0]*ds[1]
        if coefficient:f.append([sum(1<<e for e in selected),coefficient])
    assert u and f and all(c>0 for m,c in u+f)
    return sorted(u),sorted(f)


def positive(rows):
    out={}
    for row in rows:
        assert len(row) in (3,4)
        mask,slot,channel=row[:3];assert slot==0
        coefficient=Q(row[3] if len(row)==4 else 1)
        assert coefficient.denominator==1 and coefficient>0
        out[mask]=out.get(mask,0)+coefficient.numerator
    return sorted([int(mask),int(co)] for mask,co in out.items())


def compile_original(raw,recompiled):
    gd,pd=raw['gaussian_full_descriptor'],raw['positive_full_descriptor']
    assert (gd['dimension'],gd['loops'],gd['external_legs'],gd['edges'],gd['overall_degree'])==(6,4,3,12,0)
    assert gd['normalization']=='d^6k/pi^3 per loop; no symmetry weight or coupling'
    assert raw['external_gram']==[[1,-.5],[-.5,1]] and raw['positive_invariants']==[1,1,1]
    assert recompiled['terms']==gd['terms'] and recompiled['proper_subgraphs']==gd['proper_subgraphs']
    parts=gd['proper_subgraphs'];assert len(parts) in (2,3)
    assert sorted((p['loops'],p['external_legs'],p['omega'],len(p['edge_ids'])) for p in parts)==[(j,3,0,3*j) for j in range(1,len(parts)+1)]
    indices=sorted(range(len(parts)),key=lambda i:parts[i]['loops'])
    assert all(set(parts[a]['edge_ids'])<set(parts[b]['edge_ids']) for a,b in zip(indices,indices[1:]))
    words={tuple(i for i in indices if mask&(1<<i)) for mask in range(1<<len(parts))}
    assert len(gd['terms'])==len(pd['states'])==len(words)
    assert {tuple(t['subgraphs_inner_first']) for t in gd['terms']}==words
    assert {tuple(s['word']) for s in pd['states']}==words
    states=[]
    for term in gd['terms']:
        word=term['subgraphs_inner_first'];state=next(s for s in pd['states'] if s['word']==word)
        assert term['sign']==state['sign']==(-1)**len(word)
        assert term['jet_subgraphs']==term['jet_orders']==state['variables']==[]
        assert all(v['powers']==[] and type(v['coefficient']) is int for row in term['scaled_routing'] for v in row)
        routes=[[v['coefficient'] for v in row] for row in term['scaled_routing']]
        u,f=minors(routes)
        assert u==positive(state['U']) and f==positive(state['F'])
        assert all(m.bit_count()==4 for m,c in u) and all(m.bit_count()==5 for m,c in f)
        states.append({'word':word,'sign':term['sign'],'routes':routes,'U':u,'F':f})
    triangle=indices[0];pairs=[]
    for index,state in enumerate(states):
        if triangle in state['word']:continue
        target=sorted([*state['word'],triangle],key=lambda i:parts[i]['loops'])
        projected=next(i for i,s in enumerate(states) if s['word']==target)
        assert states[projected]['sign']==-state['sign']
        pairs.append([index,projected])
    assert len(pairs)*2==len(states)
    return {'id':raw['id'],'states':states,'pairs':pairs,'triangle_index':triangle,
        'original_edges':12,'original_loops':4,'physical_dimension':11,'sampling_dimension':12,
        'complete_density_homogeneity':-12,'fresh_graph_BPH_equals_original':True,
        'gaussian_Cauchy_Binet_equals_positive_forests':True,
        'all_overall_and_nested_T0_retained':True,'added_symmetry_or_factorial':False,
        'radial_kernel':'-log1p(F/(M*U))/U^3; M=sum(alpha), m^2=1',
        'paired_kernel':'s*((1/Uz^3-1/Ub^3)*log1p(rz)-log1p((rb-rz)/(1+rz))/Ub^3)'}


def polynomial(terms,alpha):
    return sum((co*math.prod(alpha[e] for e in range(12) if mask&(1<<e)) for mask,co in terms),Q())


def gaussian_uf(state,alpha):
    """Independent exact Gaussian elimination; no monomial evaluation."""
    routes=state['routes']
    matrix=[[sum((a*r[i]*r[j] for a,r in zip(alpha,routes)),Q()) for j in range(6)] for i in range(6)]
    loop=[row[:4] for row in matrix[:4]];rhs=[row[4:] for row in matrix[:4]];u=Q(1)
    for k in range(4):
        pivot=loop[k][k];assert pivot>0;u*=pivot
        for i in range(k+1,4):
            multiplier=loop[i][k]/pivot
            for j in range(k+1,4):loop[i][j]-=multiplier*loop[k][j]
            for j in range(2):rhs[i][j]-=multiplier*rhs[k][j]
    solved=[[Q(),Q()] for _ in range(4)]
    for i in range(3,-1,-1):
        for j in range(2):solved[i][j]=(rhs[i][j]-sum((loop[i][k]*solved[k][j] for k in range(i+1,4)),Q()))/loop[i][i]
    schur=[[matrix[4+i][4+j]-sum((matrix[k][4+i]*solved[k][j] for k in range(4)),Q()) for j in range(2)] for i in range(2)]
    return u,u*(schur[0][0]+schur[1][1]-(schur[0][1]+schur[1][0])/2)


def prepare_pairs(binding,alpha,use_gaussian=False):
    assert len(alpha)==12 and all(type(a) in (int,Q) and a>0 for a in alpha)
    mass=sum(alpha);values=[]
    for state in binding['states']:
        u,f=gaussian_uf(state,alpha) if use_gaussian else (polynomial(state['U'],alpha),polynomial(state['F'],alpha))
        assert u>0 and f>0;values.append((u,f/(mass*u)))
    records=[]
    for bare,zero in binding['pairs']:
        ub,rb=values[bare];uz,rz=values[zero];sign=binding['states'][bare]['sign']
        relative=(rb-rz)/(1+rz);assert relative>-1
        records.append({'base_argument':rz,'base_coefficient':sign*(1/uz**3-1/ub**3),
            'relative_argument':relative,'relative_coefficient':-sign/ub**3})
    return records


def table(binding,power=2):
    gamma={}
    for mask in range(1,4095):
        estimates=[]
        for state in binding['states']:
            nu=min((mask&m).bit_count() for m,c in state['U'])
            nf=min((mask&m).bit_count() for m,c in state['F'])
            assert nf>=nu;estimates.append(mask.bit_count()+nf-4*nu)
        gamma[mask]=min(estimates)
    return H.build(12,gamma,power)


def verify_table(table):
    """Every nonempty subset has positive transitions; induction covers 12! paths."""
    n=table['edges'];full=(1<<n)-1;assert n==12
    assert set(table['weights'])==set(range(full)) and table['weights'][0]==table['scale']
    assert set(table['shapes'])==set(range(1,full)) and min(table['weights'].values())>0
    probabilities={full:Q(1)};counts={full:1};transitions=0
    for size in range(n,0,-1):
        for mask in range(1,full+1):
            if mask.bit_count()!=size:continue
            options=[mask^(1<<e) for e in range(n) if mask&(1<<e)]
            total=sum(table['weights'][s] for s in options)
            if mask==full:assert total==table['normalization']
            else:
                shape=table['shapes'][mask];beta=Q(shape) if shape>0 else Q(1,-shape)
                assert shape>=1 or shape in (-2,-4)
                assert table['weights'][mask]*beta==total
            assert total<2**62
            for child in options:
                probabilities[child]=probabilities.get(child,Q())+probabilities[mask]*Q(table['weights'][child],total)
                counts[child]=counts.get(child,0)+counts[mask];transitions+=1
    assert probabilities[0]==1 and counts[0]==math.factorial(n)
    return {'positive_subset_partitions':full,'transitions':transitions,'complete_order_paths':counts[0],
        'path_mass':probabilities[0],'method':'exact finite directed-acyclic recurrence, not explicit permutation enumeration'}
