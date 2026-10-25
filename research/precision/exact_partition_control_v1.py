"""Piecewise Legendre controls with exact rational coefficients and exact means."""
from fractions import Fraction as Q
from itertools import combinations
import math
import exact_partition_transport_v1 as P

def powers(degree):
    assert degree in (0,1,2);result=[(0,)*7]
    if degree:
        for axis in range(7):p=[0]*7;p[axis]=1;result.append(tuple(p))
    if degree==2:
        for axis in range(7):p=[0]*7;p[axis]=2;result.append(tuple(p))
        for a,b in combinations(range(7),2):p=[0]*7;p[a]=p[b]=1;result.append(tuple(p))
    return result

def decode(row):
    assert type(row) is list and len(row)==2 and all(type(x) is int for x in row) and row[1]>0
    assert max(abs(x).bit_length() for x in row)<=2048
    return Q(*row)

def validate(model):
    assert type(model) is dict and set(model)=={'partition','degree','coefficients'}
    rows=P.leaves(model['partition']);basis=powers(model['degree'])
    assert len(rows)==len(model['coefficients'])
    coefficients=[[decode(q) for q in c] for c in model['coefficients']]
    assert all(len(c)==len(basis) for c in coefficients)
    mean=sum((math.prod(b-a for a,b in zip(row['low'],row['high']))*c[0] for row,c in zip(rows,coefficients)),Q())
    return rows,basis,coefficients,mean

def value(model,point):
    rows,basis,coefficients,mean=validate(model)
    assert len(point)==7 and all(type(x) in (int,Q) and 0<x<1 for x in point)
    selected=[i for i,row in enumerate(rows) if all(a<=x<b for a,b,x in zip(row['low'],row['high'],point))]
    assert len(selected)==1;i=selected[0];row=rows[i]
    z=[2*(x-a)/(b-a)-1 for a,b,x in zip(row['low'],row['high'],point)]
    legendre=[[Q(1),x,(3*x*x-1)/2] for x in z]
    full=sum((coefficient*math.prod(legendre[j][degree] for j,degree in enumerate(p)) for p,coefficient in zip(basis,coefficients[i])),Q())
    return full-mean

def exact_checks(model):
    rows,basis,coefficients,mean=validate(model)
    # Integral of P_0,P_1,P_2 over normalized [-1,1] is 1,0,0.
    moments=[Q(1),Q(),(3*Q(1,3)-1)/2];assert moments==[1,0,0]
    integrated=Q();terms=0
    for row,coeffs in zip(rows,coefficients):
        volume=math.prod(b-a for a,b in zip(row['low'],row['high']))
        for p,co in zip(basis,coeffs):integrated+=volume*co*math.prod(moments[k] for k in p);terms+=1
    assert integrated==mean and integrated-mean==0
    return {'leaves':len(rows),'degree':model['degree'],'basis_size':len(basis),'exact_terms_integrated':terms,
            'known_polynomial_mean':mean,'exact_control_integral':Q(),'physics_counterterms_changed':False}
