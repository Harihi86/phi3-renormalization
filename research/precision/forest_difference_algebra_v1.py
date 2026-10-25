"""Factored cube differences in x_i²=r_i*x_i; no division by a child scale."""
from itertools import product


def constant(value,k):return [value]+[0]*((1<<k)-1)
def variable(i,k,weight=1):
    result=constant(0,k);result[1<<i]=weight;return result
def add(a,b):return [x+y for x,y in zip(a,b)]
def subtract(a,b):return [x-y for x,y in zip(a,b)]
def scale(a,c):return [c*x for x in a]
def monomial(r,mask):
    out=1
    for i,x in enumerate(r):
        if mask&(1<<i):out*=x
    return out


def multiply(a,b,r):
    out=[0]*len(a);weights=[monomial(r,i) for i in range(len(a))]
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i|j]+=x*y*weights[i&j]
    return out


def vertex(a,r,mask):return sum(v*monomial(r,i) for i,v in enumerate(a) if i&mask==i)
def restrict(a,mask):return [v if i&mask==i else 0 for i,v in enumerate(a)]


def inverse(a,r):
    out=[0]*len(a)
    for mask in sorted(range(len(a)),key=lambda x:(x.bit_count(),x)):
        if mask==0:out[0]=1/a[0];continue
        numerator=0
        for b in range(len(a)):
            if b==mask or b&mask!=b:continue
            for c in range(len(a)):
                if c&mask==c and c|b==mask:numerator+=a[c]*out[b]*monomial(r,c&b)
        out[mask]=-numerator/vertex(a,r,mask)
    return out


def power(a,n,r):
    out=constant(1,len(r))
    for _ in range(n):out=multiply(out,a,r)
    return out


def log_ratio_data(a,r,mask):
    assert mask
    even=[];odd=[]
    for sub in range(len(a)):
        if sub&mask==sub:(even if (mask.bit_count()-sub.bit_count())%2==0 else odd).append(sub)
    products=[];values=[]
    for group in (even,odd):
        p=constant(1,len(r));v=1
        for sub in group:p=multiply(p,restrict(a,sub),r);v*=vertex(a,r,sub)
        products.append(p);values.append(v)
    numerator=subtract(*products)[mask];denominator=values[1]
    t=monomial(r,mask)*numerator/denominator
    return numerator/denominator,t,values[0],values[1]


def logarithm(a,r,log,divided_log):
    out=[log(a[0])]+[0]*(len(a)-1)
    for mask in range(1,len(a)):
        factor,t,positive,negative=log_ratio_data(a,r,mask)
        out[mask]=factor*divided_log(t,positive,negative)
    return out


def forest_kernel(U,F,mass,r,log,divided_log,meanF=None):
    inv=inverse(U,r);P=add(scale(U,mass),F)
    logs=subtract(logarithm(P,r,log,divided_log),logarithm(U,r,log,divided_log))
    logs[0]-=log(mass)
    result=multiply(power(inv,3,r),logs,r)
    if meanF is not None:
        correction=scale(multiply(subtract(meanF,F),power(inv,4,r),r),1/mass)
        result=add(result,correction)
    # Original -Delta_1...Delta_k f / product r_i, with mass fixed across corners.
    return -result[-1]


def polynomial(rows,edge_rings,r):
    total=constant(0,len(r))
    for mask,co in rows:
        term=constant(co,len(r))
        for e,edge in enumerate(edge_rings):
            if mask&(1<<e):term=multiply(term,edge,r)
        total=add(total,term)
    return total
