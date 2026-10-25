"""Exact non-quadrature allowances for the radial-nested routing of G380."""
from fractions import Fraction as Q
from math import factorial,isqrt
from mixed_triangle_discrete_bounds import shifted_sum,weights
from exact_exp2_enclosures import log2

def bounds(N=48):
    assert N in (48,64);q=Q(5,12);K=2*N-24
    A=shifted_sum(K,3,q)/288
    coefficient=sum(abs(w)*Q(ell+2,4*(ell+1))*A for ell,w in enumerate(weights(20)))
    highL=49*log2()[1]
    def tail(power):return sum(Q(factorial(power),factorial(power-j))*highL**(power-j) for j in range(power+1))/2**49
    errors={'projection':coefficient*q**(K-6)/32*(Q(1,3)+Q(2,K-6)),
            'harmonics':q**36*shifted_sum(21,1,q*q)/23040,
            'lower_endpoint':Q(1,2**50*5760),
            'upper_endpoint':24*tail(2)+Q(1312,3)*tail(3)}
    errors['total']=sum(errors.values())
    return {'id':'L4E3G380','angular_grid':N,'degrees':20,'domain':(Q(1,1024),Q(2)**49),'errors':errors}

def sqrt_upper(q,bits=160):
    q=Q(q);assert q>=0 and max(q.numerator.bit_length(),q.denominator.bit_length())<1000000
    n=isqrt((q.numerator<<(2*bits))//q.denominator)
    if n*n*q.denominator<q.numerator<<(2*bits):n+=1
    answer=Q(n,1<<bits);assert answer*answer>=q
    return answer
