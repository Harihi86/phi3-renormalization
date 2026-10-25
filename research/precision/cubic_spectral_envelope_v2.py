"""Exact cubic enclosures from complete monotonicity; no sampled extrema."""
from fractions import Fraction as Q
from math import comb
from convex_radial_envelope import Node
from convex_spectral_envelope import SpectralEnvelope
def add(a,b):return tuple((a[i] if i<len(a) else Q())+(b[i] if i<len(b) else Q()) for i in range(max(len(a),len(b))))
def scale(a,c):return tuple(Q(c)*v for v in a)
def mul(a,b):
    out=[Q()]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return tuple(out)
def value(a,x):
    result=Q()
    for c in reversed(a):result=result*x+c
    return result
def bernstein(a,left,right,degree=3):
    assert len(a)<=degree+1 and left<right;h=right-left
    powers=[sum((a[k]*comb(k,j)*left**(k-j)*h**j for k in range(j,len(a))),Q()) for j in range(degree+1)]
    return tuple(sum((powers[j]*Q(comb(i,j),comb(degree,j)) for j in range(i+1)),Q()) for i in range(degree+1))
def interpolate(nodes,indices,left,right,upper):
    mid=(left+right)/2;result=(Q(),);choices=[]
    for j in indices:
        basis=(Q(1),)
        for k in indices:
            if k!=j:basis=scale(mul(basis,(-nodes[k].x,Q(1))),1/(nodes[j].x-nodes[k].x))
        positive=value(basis,mid)>0
        y=nodes[j].upper if positive==upper else nodes[j].lower
        result=add(result,scale(basis,y));choices.append({'node':j,'upper_endpoint_used':positive==upper})
    return result,choices
def dyadic_upper(q,bits=160):
    q=Q(q);assert q>=0;scale=1<<bits
    return Q((q.numerator*scale+q.denominator-1)//q.denominator,scale)
def inspect(nodes):
    SpectralEnvelope(nodes)  # Validate ordering and decreasing compatibility.
    assert len(nodes)>=4;pieces=[];by_interval={}
    for i,(a,b) in enumerate(zip(nodes,nodes[1:])):
        left,right=a.x,b.x;lower_candidates=[((Q(),),[],'positive_zero'),((b.lower,),[{'node':i+1,'upper_endpoint_used':False}],'right_constant')]
        if i>=1:
            prev=nodes[i-1];s=(a.lower-prev.upper)/(a.x-prev.x);lower_candidates.append(((a.lower-s*a.x,s),[{'node':i,'upper_endpoint_used':False},{'node':i-1,'upper_endpoint_used':True}],'left_secant'))
        if i+2<len(nodes):
            after=nodes[i+2];s=(after.upper-b.lower)/(after.x-b.x);lower_candidates.append(((b.lower-s*b.x,s),[{'node':i+1,'upper_endpoint_used':False},{'node':i+2,'upper_endpoint_used':True}],'right_secant'))
        if i>=1 and i+2<len(nodes):
            p,ch=interpolate(nodes,list(range(i-1,i+3)),left,right,False);lower_candidates.append((p,ch,'cubic_middle'))
        upper,choice=interpolate(nodes,[i,i+1],left,right,True);upper_candidates=[(upper,choice,'linear_chord')]
        if i>=2:
            p,ch=interpolate(nodes,list(range(i-2,i+2)),left,right,True);upper_candidates.append((p,ch,'cubic_right'))
        if i+3<len(nodes):
            p,ch=interpolate(nodes,list(range(i,i+4)),left,right,True);upper_candidates.append((p,ch,'cubic_left'))
        candidates=[]
        for L,lc,lt in lower_candidates:
            for U,uc,ut in upper_candidates:
                gap=add(U,scale(L,-1));gapmax=max(bernstein(gap,left,right));assert gapmax>=0
                error=gapmax/2;M=scale(add(L,U),Q(1,2));Mmax=max(map(abs,bernstein(M,left,right)))
                fmax=Mmax+max(left**-4,left**-2)/18;weight=right**5*min(right**4/24,right**2/3)
                propagated=dyadic_upper((right-left)*weight*(2*fmax*error+error*error))
                candidates.append((propagated,lt,ut,L,U,M,error,lc,uc))
        propagated,lt,ut,L,U,M,error,lc,uc=min(candidates,key=lambda x:(x[0],x[1],x[2]))
        by_interval[i]=propagated;pieces.append({'left':left,'right':right,'lower':L,'upper':U,'midpoint':M,'lower_kind':lt,'upper_kind':ut,'lower_choices':lc,'upper_choices':uc,'C_error':error,'density_error_integral_bound':propagated,'node_interval':i})
    return {'pieces':pieces,'by_interval':by_interval,'propagated_error':sum(by_interval.values(),Q())}
