"""Literal original-route matrix reference, independent of native kernels/CDF."""
from fractions import Fraction as Q
from itertools import combinations,permutations
import sympy as S
import primitive_triangle_factor_v3 as V


def prepare(raw):
    gd=raw['gaussian_full_descriptor'];selected=gd['proper_subgraphs'][0]['edge_ids'];assert len(selected)==3
    _,witnesses=V.basis_all(gd,0,S.symbols('z:1'))
    rows=next(w['routing'] for t,w in zip(gd['terms'],witnesses) if not t['subgraphs_inner_first'])
    routes=[[int(v) for v in row[:,1:]] for e,row in enumerate(rows) if e not in selected]
    shifts=[rows[e][:,1:]/rows[e][0] for e in selected]
    delta=[[int(v) for v in shifts[i]-shifts[j]] for i,j in combinations(range(3),2)]
    buckets={}
    for e,row in enumerate(routes):
        sign=next(v for v in row[:3] if v);key=tuple(sign*v for v in row[:3]);buckets.setdefault(key,[]).append(e)
    groups=list(buckets.values());assert len(groups)==6
    indices={g:j for j,g in enumerate(g for g,edges in enumerate(groups) if len(edges)==2)}
    anchor=next(g for g,edges in enumerate(groups) if 8 in edges)
    nu={}
    for mask in range(64):
        rowsleft=[routes[edges[0]][:3] for g,edges in enumerate(groups) if not mask&(1<<g)]
        rank=S.Matrix(rowsleft).rank() if rowsleft else 0;nu[mask]=3-rank
    gamma={mask:sum(len(groups[g]) for g in range(6) if mask&(1<<g))-3*nu[mask] for mask in range(1,63)}
    assert min(gamma.values())>0;literal=[]
    for order in permutations(range(6)):
        mask=63;mass=Q(1);exponents=[];leading=[]
        for e in order[:-1]:mask^=1<<e;mass/=gamma[mask];exponents.append(gamma[mask]);leading.append(nu[mask])
        literal.append({'order':order,'mass':mass,'exponents':exponents,'leading':leading})
    normalization=sum((row['mass'] for row in literal),Q());strata=[]
    for first in range(6):
        for rank in range(6):
            if (first==anchor)!=(rank==0):continue
            orders=[row for row in literal if row['order'][0]==first and row['order'][rank]==anchor]
            mass=sum((row['mass'] for row in orders),Q());strata.append({'first':first,'position':rank,'probability':mass/normalization,'mass':mass,'orders':orders})
    assert sum((row['probability'] for row in strata),Q())==1
    return {'id':raw['id'],'routes':routes,'delta':delta,'groups':groups,'indices':indices,
            'anchor':anchor,'normalization':normalization,'strata':strata,'gamma':gamma,'nu':nu}


def decode(z,stratum):
    z=Q(z);assert 0<=z<1;left=Q()
    for row in stratum['orders']:
        probability=row['mass']/stratum['mass']
        if z<left+probability:return {'order':list(row['order']),'exponents':row['exponents'],'leading':row['leading'],
                                    'auxiliary':(z-left)/probability,'probability':probability}
        left+=probability
    raise AssertionError('literal triangle sector CDF failed')


def evaluate(mp,points,control,prepared):
    p=list(map(Q,points));mq=lambda q:mp.mpf(q.numerator)/q.denominator
    gram=mp.matrix([[1,-mp.mpf('.5')],[-mp.mpf('.5'),1]]);routes=prepared['routes'];delta=prepared['delta'];groups=prepared['groups']
    def schur(matrix):
        A=matrix[:3,:3];inverse=A**-1;U=mp.det(A)
        remaining=matrix[3:,3:]-matrix[3:,:3]*inverse*matrix[:3,3:]
        q=sum(gram[i,j]*remaining[j,i] for i in range(2) for j in range(2))
        return U,q,inverse
    def host(a,x):
        alpha=[mp.mpf(0)]*9
        for g,edges in enumerate(groups):
            if len(edges)==1:alpha[edges[0]]=a[g]
            else:
                z=x[prepared['indices'][g]];alpha[edges[0]]=a[g]*z;alpha[edges[1]]=a[g]*(1-z)
        matrix=mp.matrix(5)
        for aa,row in zip(alpha,routes):
            r=mp.matrix(row);matrix+=aa*r*r.T
        return matrix
    def mean(a,center_q,inverse):
        value=center_q
        for g,edges in enumerate(groups):
            if len(edges)==2:
                l=mp.matrix(routes[edges[0]][:3]);value-=a[g]**2*(l.T*inverse*l)[0]/12
        return value
    def exprel(z):return mp.expm1(z)/z if z else mp.mpf(1)
    total=mp.mpf(0)
    for st in prepared['strata']:
        decoded=decode(p[0],st);a=[mp.mpf(1)]*6;cumulative=tropical=mp.mpf(1)
        for j,e in enumerate(decoded['order'][1:]):
            ratio=mp.root(mq(p[j+1]),decoded['exponents'][j]);cumulative*=ratio;a[e]=cumulative;tropical*=ratio**decoded['leading'][j]
        M=sum(a);H=host(a,list(map(mq,p[6:9])));Hc=host(a,[mp.mpf('.5')]*3)
        U0,q0,inverse=schur(H);mean0=mean(a,schur(Hc)[1],inverse) if control else mp.mpf(0)
        square_root=mp.sqrt(1-mq(p[9]));simplex=[mq(p[9])/(1+square_root),square_root*mq(p[10]),square_root*(1-mq(p[10]))]
        value=mp.mpf(0)
        for child in permutations(simplex):
            variance=mp.matrix(5)
            for (i,j),d in zip(combinations(range(3),2),delta):
                row=mp.matrix(d);variance+=child[i]*child[j]*row*row.T
            vl=variance[:3,:3];ell=sum((inverse*vl)[i,i] for i in range(3));assert ell>0;h=mp.log(M*ell)
            # Original Gaussian residual matrices provide the exact radial endpoint limits.
            shift=mp.matrix(5,2)
            for i in range(3):
                for j in range(2):shift[i,j]=-(inverse*H[:3,3:])[i,j]
            shift[3,0]=shift[4,1]=1
            De=sum(gram[i,j]*(shift.T*variance*shift)[j,i] for i in range(2) for j in range(2))
            meanDe=mp.mpf(0)
            if control:
                center_shift=mp.matrix(shift)
                for i in range(3):
                    for j in range(2):center_shift[i,j]=-(inverse*Hc[:3,3:])[i,j]
                meanDe=sum(gram[i,j]*(center_shift.T*variance*center_shift)[j,i] for i in range(2) for j in range(2))
                for g,edges in enumerate(groups):
                    if len(edges)==2:
                        l=mp.matrix(routes[edges[0]][:3]);meanDe+=a[g]**2*(l.T*inverse*vl*inverse*l)[0]/12
            for uQ in (decoded['auxiliary'],1-decoded['auxiliary']):
                u=mq(uQ)
                if uQ==0:
                    current=M/exprel(h)*(3*ell*mp.log1p(q0/M)-De/(M+q0))
                    if control:current+=M/exprel(h)*(3*ell*(mean0-q0)-(meanDe-De))/M
                elif uQ==1:current=(mean0 if control else q0)/(M*exprel(-h))
                else:
                    radial=M*u/(1-u)*exprel(u*h)/(mp.exp(u*h)*exprel((1-u)*h))
                    jac=exprel(h)/(u*(1-u)*exprel(u*h)*exprel((1-u)*h));mass=M+radial
                    Ur,qr,inverser=schur(H+radial*variance);inv3=(U0/Ur)**3
                    current=mp.log1p(q0/mass)-inv3*mp.log1p(qr/mass)
                    if control:
                        meanr=mean(a,schur(Hc+radial*variance)[1],inverser)
                        current+=((mean0-q0)-inv3*(meanr-qr))/mass
                    current*=jac
                value+=current
        total+=mq(st['probability']*prepared['normalization'])*(tropical/U0)**3*value/24
    return total
