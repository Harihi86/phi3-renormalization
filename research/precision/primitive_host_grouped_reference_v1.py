"""Independent literal 720-order measure and original-route MP reference.

No production proposal, positive polynomial, decoder or native runtime imports.
The Gaussian construction shares only the exact original loop-basis helper.
The elementary bubble is checked separately against its hypergeometric integral.
"""
from fractions import Fraction as Q
from itertools import permutations,combinations
import sympy as S
import primitive_triangle_factor_v3 as V


def prepare(raw):
    gd=raw['gaussian_full_descriptor'];selected=gd['proper_subgraphs'][0]['edge_ids']
    _,witnesses=V.basis_all(gd,0,S.symbols('z:1'))
    witness=next(w for t,w in zip(gd['terms'],witnesses) if not t['subgraphs_inner_first'])
    full=witness['routing'];delta=full[selected[0]][:,1:]/full[selected[0]][0]-full[selected[1]][:,1:]/full[selected[1]][0]
    outside=[e for e in range(12) if e not in selected]
    repeated=[e for e in outside if full[e][:,1:] in (delta,-delta)];assert len(repeated)==2
    kept=[e for e in outside if e not in repeated]+[repeated[0]]
    routes=[[int(v) for v in full[e][:,1:]] for e in kept]
    buckets={}
    for e,row in enumerate(routes):
        first=next(v for v in row[:3] if v);key=tuple(first*v for v in row[:3]);buckets.setdefault(key,[]).append(e)
    groups=list(buckets.values());assert len(groups)==6
    for edges in groups:
        if 8 in edges:edges.sort(key=lambda e:e!=8)
    weighted=next(g for g,edges in enumerate(groups) if 8 in edges)
    indices={g:j for j,g in enumerate(g for g,edges in enumerate(groups) if len(edges)==2)}
    bases=[]
    for comb in combinations(range(6),3):
        matrix=S.Matrix([routes[groups[g][0]][:3] for g in comb]);det=matrix.det()
        if det:assert abs(det)==1;bases.append(comb)
    nu={mask:min(sum(bool(mask&(1<<g)) for g in comb) for comb in bases) for mask in range(64)}
    gamma={mask:sum(len(groups[g]) for g in range(6) if mask&(1<<g))-3*nu[mask] for mask in range(1,63)}
    assert min(gamma.values())>0
    literal=[]
    for order in permutations(range(6)):
        mask=63;mass=Q(1);exponents=[];leading=[]
        for e in order[:-1]:mask^=1<<e;mass/=gamma[mask];exponents.append(gamma[mask]);leading.append(nu[mask])
        literal.append({'order':order,'mass':mass,'exponents':exponents,'leading':leading})
    normalization=sum((entry['mass'] for entry in literal),Q());strata=[]
    for first in range(6):
        for rank in range(6):
            if (first==weighted)!=(rank==0):continue
            orders=[entry for entry in literal if entry['order'][0]==first and entry['order'][rank]==weighted]
            mass=sum((entry['mass'] for entry in orders),Q())
            strata.append({'first':first,'position':rank,'probability':mass/normalization,'mass':mass,'orders':orders})
    assert sum((s['probability'] for s in strata),Q())==1
    return {'id':raw['id'],'routes':routes,'groups':groups,'indices':indices,'weighted':weighted,
            'strata':strata,'normalization':normalization,'gamma':gamma,'nu':nu}


def decode(z,stratum):
    z=Q(z);assert 0<=z<1;left=Q()
    for row in stratum['orders']:
        p=row['mass']/stratum['mass']
        if z<left+p:
            return {'order':list(row['order']),'exponents':row['exponents'],'leading':row['leading'],
                    'auxiliary':(z-left)/p,'probability':p}
        left+=p
    raise AssertionError('literal permutation CDF failed')


def bubble(mp,t):
    # Algebraically equal to the hypergeometric renormalized bubble; no native code.
    return (t+4)**mp.mpf('1.5')/(3*mp.sqrt(t))*mp.asinh(mp.sqrt(t)/2)-(t+4)/3-t/9


def evaluate(mp,points,control,prepared,scalar_arguments=None):
    mq=lambda q:mp.mpf(q.numerator)/q.denominator
    p=list(map(Q,points));total=mp.mpf(0);groups=prepared['groups'];routes=prepared['routes']
    gram=mp.matrix([[1,-mp.mpf('.5')],[-mp.mpf('.5'),1]])
    def gaussian(alpha):
        matrix=mp.matrix(5)
        for a,row in zip(alpha,routes):
            for i in range(5):
                for j in range(5):matrix[i,j]+=a*row[i]*row[j]
        A=matrix[:3,:3];U=mp.det(A);inv=A**-1
        schur=matrix[3:,3:]-matrix[3:,:3]*inv*matrix[:3,3:]
        q=sum(gram[i,j]*schur[j,i] for i in range(2) for j in range(2))
        return U,q,inv
    def alphas(a,x):
        alpha=[mp.mpf(0)]*9
        for g,edges in enumerate(groups):
            if len(edges)==1:alpha[edges[0]]=a[g]
            else:
                fraction=x[prepared['indices'][g]];alpha[edges[0]]=a[g]*fraction;alpha[edges[1]]=a[g]*(1-fraction)
        return alpha
    for st in prepared['strata']:
        decoded=decode(p[0],st);a=[mp.mpf(1)]*6;cumulative=tropical=mp.mpf(1)
        for j,e in enumerate(decoded['order'][1:]):
            ratio=mp.root(mq(p[j+1]),decoded['exponents'][j]);cumulative*=ratio;a[e]=cumulative;tropical*=ratio**decoded['leading'][j]
        x=list(map(mq,p[6:]));jac=mp.mpf(1)
        if len(groups[prepared['weighted']])==2:
            j=prepared['indices'][prepared['weighted']];jac=2*x[j];x[j]*=x[j]
        alpha=alphas(a,x);U,q,inv=gaussian(alpha);M=sum(a);d=alpha[8];other=sum(alpha[:8]);u=mq(decoded['auxiliary'])
        value=mp.mpf(0)
        for aux in (u,1-u):
            mass=M+aux*q;arg=(other+aux*q)/d
            if scalar_arguments is not None:scalar_arguments.append(arg)
            value-=d*q*bubble(mp,arg)/(2*mass**2)
        if control:
            center=list(x)
            for g,edges in enumerate(groups):
                if len(edges)==2 and g!=prepared['weighted']:center[prepared['indices'][g]]=mp.mpf('.5')
            _,mean,_=gaussian(alphas(a,center))
            for g,edges in enumerate(groups):
                if len(edges)==2 and g!=prepared['weighted']:
                    l=mp.matrix(routes[edges[0]][:3]);mean-=a[g]**2*(l.T*inv*l)[0]/12
            value+=d*bubble(mp,other/d)/M**2*(q-mean)
        total+=mq(st['probability']*prepared['normalization'])*(tropical/U)**3*jac*value
    return total
