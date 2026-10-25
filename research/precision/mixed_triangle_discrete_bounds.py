"""Exact, independent error allowances for the finite angular projection."""
from fractions import Fraction as Q
from math import comb,factorial
from exact_exp2_enclosures import log2

def moments(x):
    return [1/(1-x),x/(1-x)**2,x*(1+x)/(1-x)**3,
            x*(1+4*x+x*x)/(1-x)**4,
            x*(1+11*x+11*x*x+x**3)/(1-x)**5]

def shifted_sum(start,power,x):
    return sum(Q(comb(power,j))*(start+3)**(power-j)*v for j,v in enumerate(moments(x)[:power+1]))

def weights(L):
    c=[Q(1),Q(2)]
    for n in range(1,L):c.append(((n+2)*c[n]-(n+3)*c[n-1])/(n+1))
    return [2*c[n]/(n+2) for n in range(L+1)]

def bounds(N=48,L=20,left=-12,right=64):
    assert N in (48,64) and L==20 and left<=0 and right>=2
    K=2*N-L-4;q0=Q(5,12);assert K>=72
    A1=shifted_sum(K,3,q0)/288;A2=shifted_sum(K,4,q0)/864
    e205=e220=Q()
    for ell,w in enumerate(weights(L)):
        A=A1*Q(ell+2,2);B=A2*Q(ell+2,2)
        e205+=abs(w)*(B/10+A/4+A*B*q0**K)
        e220+=abs(w)*(A/5+A*A*q0**K)
    jbare=q0**(K-4)*(Q(1,4)+Q(1,K-4))
    jone=q0**(K-6)/120*(Q(1,3)+Q(2,K-6))
    eps=Q(2)**left;R=Q(2)**right;lnR=right*log2()[1]
    def Itail(k):return sum(Q(factorial(k),factorial(k-j))*lnR**(k-j) for j in range(k+1))/R
    ans={
        'L4E3G205':{'projection':e205*jbare,'harmonics':q0**38*shifted_sum(21,3,q0*q0)/62208,'lower_endpoint':eps**3/162000,'upper_endpoint':8*Itail(2)+Q(4784,27)*Itail(3)},
        'L4E3G220':{'projection':e220*jone,'harmonics':q0**36*shifted_sum(21,2,q0*q0)/2073600,'lower_endpoint':eps**3/648000,'upper_endpoint':Q(16,3)*Itail(2)+Q(920,9)*Itail(3)}}
    for row in ans.values():row['total']=sum(row.values())
    return {'angular_grid':N,'degrees':L,'first_possible_alias':K,'log2_interval':(left,right),'graphs':ans}
