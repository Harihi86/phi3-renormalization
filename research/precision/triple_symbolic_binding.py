"""Exact all-beta polynomial binding through original Gaussian routing minors.

Cauchy--Binet is applied before floating evaluation. Every Taylor monomial
through order two in each original variable is kept, including odd powers
until their cancellation is verified. No positive U/F rows enter this path.
"""
from fractions import Fraction as Q
from itertools import combinations,permutations
from collections import defaultdict


def add(out,poly,factor=Q(1)):
    for key,value in poly.items():out[key]+=value*factor


def multiply(a,b,k):
    out=defaultdict(Q)
    for p,c in a.items():
        for q,d in b.items():
            exponents=tuple(x+y for x,y in zip(p,q))
            if all(e<=2 for e in exponents):out[exponents]+=c*d
    return {p:c for p,c in out.items()if c}


def minor(route,edges,columns,k):
    result=defaultdict(Q);zero=(0,)*k
    for order in permutations(columns):
        rank=[columns.index(c)for c in order]
        sign=(-1)**sum(rank[i]>rank[j]for i in range(len(rank))for j in range(i+1,len(rank)))
        term={zero:Q(sign)}
        for edge,col in zip(edges,order):
            row=route[edge][col];power=tuple(row['powers'])
            if len(power)!=k or any(type(e)is not int or not 0<=e<=2 for e in power):raise ValueError('Complete original routing Taylor exponents')
            coefficient=Q(row['coefficient'])
            term=multiply(term,{power:coefficient}if coefficient else{},k)
            if not term:break
        add(result,term)
    return {p:c for p,c in result.items()if c}


def gaussian(descriptor,term,group,gram,target_variables):
    L=descriptor['loops'];P=descriptor['external_legs']-1;k=len(term['jet_subgraphs']);route=term['scaled_routing']
    if len(group)!=3 or len(set(group))!=3 or set(term['jet_subgraphs'])!=set(target_variables):raise ValueError('Exact selected edges and original Taylor IDs')
    remaining=tuple(e for e in range(descriptor['edges'])if e not in group)
    active=tuple(sorted((*remaining,group[0])));out=[defaultdict(Q)for _ in range(4)]
    for family,degree in [(0,L),(1,L+1)]:
        for edges in combinations(active,degree):
            if family==0:
                det=minor(route,edges,tuple(range(L)),k);poly=multiply(det,det,k)
            else:
                dets=[minor(route,edges,(*range(L),L+i),k)for i in range(P)];poly=defaultdict(Q)
                for i in range(P):
                    for j in range(P):add(poly,multiply(dets[i],dets[j],k),Q(gram[i][j]))
            mask=sum(1<<e for e in edges if e in remaining);linear=int(group[0]in edges)
            for powers,value in poly.items():
                if not value:continue
                if any(p%2 for p in powers):raise ArithmeticError('Odd Gaussian polynomial failed exact cancellation')
                slot=sum(1<<target_variables.index(v)for v,p in zip(term['jet_subgraphs'],powers)if p==2)
                out[2*family+linear][(mask,slot)]+=value
    return [[[mask,slot,c]for(mask,slot),c in sorted(poly.items())if c]for poly in out]


def positive(state,group,invariants):
    out=[defaultdict(Q)for _ in range(4)];selected=1<<group[0];discarded=sum(1<<e for e in group[1:])
    for family,f in enumerate(('U','F')):
        for row in state[f]:
            mask,slot,ch=row[:3]
            if mask&discarded:continue
            coefficient=Q(row[3]if len(row)==4 else 1)*(Q(invariants[ch])if f=='F'else 1)
            out[2*family+int(bool(mask&selected))][(mask&~selected,slot)]+=coefficient
    return [[[mask,slot,c]for(mask,slot),c in sorted(poly.items())if c]for poly in out]


def evaluate(polynomials,beta,k):
    out=[[Q()for _ in range(1<<k)]for _ in range(4)]
    for target,rows in zip(out,polynomials):
        for mask,slot,c in rows:
            v=c
            for e,x in beta.items():
                if mask&(1<<e):v*=x
            if mask&~sum(1<<e for e in beta):raise ValueError('Unbound original beta')
            target[slot]+=v
    return out
