"""Exact rational inverse-CDF transport for frozen separable histograms."""
from fractions import Fraction as Q
UNITS=1<<16


def validate(model):
 assert type(model) is dict and set(model)=={'bins','counts'}
 bins=model['bins'];rows=model['counts']
 assert type(bins) is int and bins in (1,4,8,16) and type(rows) is list and rows
 assert all(type(row) is list and len(row)==bins and all(type(n) is int and n>0 for n in row) and sum(row)==UNITS for row in rows)
 return bins,rows


def transform(model,point):
 bins,rows=validate(model);assert len(point)==len(rows)
 mapped=[];jac=Q(1)
 for y,counts in zip(point,rows):
  assert type(y) in (Q,int) and 0<=y<=1
  cumulative=0
  for j,n in enumerate(counts):
   if y<Q(cumulative+n,UNITS) or j==bins-1:break
   cumulative+=n
  slope=Q(UNITS,bins*n)
  mapped.append(Q(j,bins)+(y-Q(cumulative,UNITS))*slope);jac*=slope
 assert all(0<=x<=1 for x in mapped) and jac>0
 return mapped,jac


def exact_checks(model):
 bins,rows=validate(model);moments=0;boundaries=0
 for counts in rows:
  cumulative=0
  for j,n in enumerate(counts):
   lo,hi=Q(cumulative,UNITS),Q(cumulative+n,UNITS);slope=Q(UNITS,bins*n)
   assert Q(j,bins)+(hi-lo)*slope==Q(j+1,bins)
   assert (hi-lo)*slope==Q(1,bins)
   cumulative+=n;boundaries+=2
  assert cumulative==UNITS
  for power in range(5):
   integral=sum(((Q(j+1,bins)**(power+1)-Q(j,bins)**(power+1))/Q(power+1) for j in range(bins)),Q())
   assert integral==Q(1,power+1);moments+=1
 return {'exact_polynomial_moments':moments,'exact_cell_endpoints':boundaries,'positive_full_domain_Jacobian':True}
