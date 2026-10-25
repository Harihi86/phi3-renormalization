"""Actual subtracted rational INTEGRANDS for three representative two-loop graphs.

No quadrature, integrated values, UV limit or general convergence is asserted.
Coordinates are declared here, separate from a graph's automatically chosen routing.
Generic exact Taylor coefficients use a recurrence, not finite differences.
Fraction inputs keep all arithmetic exact.
"""
from fractions import Fraction


def dot(a,b):
    if len(a)!=len(b):raise ValueError("Vectors must have equal dimensions")
    return sum(x*y for x,y in zip(a,b))
def add(a,b):
    if len(a)!=len(b):raise ValueError("Vectors must have equal dimensions")
    return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)


def delta(q,m2=1):
    if m2<=0:raise ValueError("Positive mass squared is required")
    denominator=dot(q,q)+m2
    return Fraction(1,1)/denominator


def taylor_product(factors,order,m2=1):
    """Return coeffs of prod 1/((a+t*b)^2+m2) through t**order.

    For A=a^2+m2, B=2a.b, C=b^2, each reciprocal satisfies
    u_0=1/A; u_n=-(B*u_(n-1)+C*u_(n-2))/A.
    Evaluate the Taylor polynomial at t=1 by summing the coefficients.
    """
    if order<0:raise ValueError("Taylor order must be nonnegative")
    if m2<=0:raise ValueError("Positive mass squared is required")
    result=[Fraction(1)]+[Fraction()]*order
    for a,b in factors:
        A=dot(a,a)+m2;B=2*dot(a,b);C=dot(b,b)
        u=[Fraction(1)/A]
        for n in range(1,order+1):
            u.append(-(B*u[n-1]+(C*u[n-2] if n>1 else 0))/A)
        result=[sum(result[j]*u[n-j] for j in range(n+1)) for n in range(order+1)]
    return tuple(result)


def bubble_taylor(q,r,m2=1):
    d=delta(q,m2);qr=dot(q,r);r2=dot(r,r)
    return d*d-2*qr*d**3-r2*d**3+4*qr**2*d**4


def self_energy_insertion(k,q,p,m2=1):
    """Nested two-point graph. One bubble subgraph; four forest terms.

    Bare denominator product: D(k)^2 D(k+p) D(q) D(q+k).
    Internal Taylor subtracts q-bubble in input k through degree 2;
    overall subtraction in p is also degree 2 and leaves k,q fixed.
    """
    dk=delta(k,m2);dq=delta(q,m2)
    h=dq*delta(add(q,k),m2);ct=bubble_taylor(q,k,m2)
    a=dk**2*delta(add(k,p),m2)
    ta=dk**2*sum(taylor_product([(k,p)],2,m2))
    terms=[a*h,-a*ct,-ta*h,ta*ct]
    return {"terms":terms,"sum":sum(terms),"factorized":(a-ta)*(h-ct),
            "operators":["1","-t_gamma^2","-t_G^2","+t_G^2 t_gamma^2"]}


def triangle_insertion(k,q,p1,p2,m2=1):
    """Nested three-point graph. Subgraph and overall Taylor orders are zero.

    a=D(k)D(k+p2)D(k-p1); h=D(q)D(q+k)D(q+k-p1).
    t_gamma sets the independent subgraph inputs (k,p1) to zero, not q.
    t_G sets physical (p1,p2) to zero, retaining outer loop k.
    """
    dk=delta(k,m2);dq=delta(q,m2);dqk=delta(add(q,k),m2)
    a=dk*delta(add(k,p2),m2)*delta(add(k,neg(p1)),m2)
    h=dq*dqk*delta(add(add(q,k),neg(p1)),m2)
    terms=[a*h,-a*dq**3,-dk**3*dq*dqk**2,dk**3*dq**3]
    return {"terms":terms,"sum":sum(terms),
            "factorized":a*(h-dq**3)-dk**3*(dq*dqk**2-dq**3),
            "operators":["1","-t_gamma^0","-t_G^0","+t_G^0 t_gamma^0"]}


def overlapping_self_energy(k,q,p,m2=1):
    """Diamond graph: two overlapping triangle subgraphs; SIX forest terms.

    I=D(k)D(p-k)D(q)D(p-q)D(k-q).
    Subtractions: D(k)^3 D(q)D(p-q), and its k<->q partner.
    The forbidden t_gamma1 t_gamma2 term is NOT constructed.
    The remaining overall operation is Taylor degree 2 in physical p.
    """
    zero=(0,)*len(k);dk=delta(k,m2);dq=delta(q,m2);dkq=delta(add(k,neg(q)),m2)
    dpk=delta(add(p,neg(k)),m2);dpq=delta(add(p,neg(q)),m2)
    raw=dk*dpk*dq*dpq*dkq
    c1=dk**3*dq*dpq;c2=dq**3*dk*dpk
    overall_raw=dk*dq*dkq*sum(taylor_product([(neg(k),p),(neg(q),p)],2,m2))
    overall_c1=dk**3*dq*sum(taylor_product([(neg(q),p)],2,m2))
    overall_c2=dq**3*dk*sum(taylor_product([(neg(k),p)],2,m2))
    terms=[raw,-c1,-c2,-overall_raw,overall_c1,overall_c2]
    return {"terms":terms,"sum":sum(terms),
            "operators":["1","-t_gamma1^0","-t_gamma2^0","-t_G^2",
                         "+t_G^2 t_gamma1^0","+t_G^2 t_gamma2^0"]}


def demo():
    k=(Fraction(1),Fraction(2),0,0,0,0)
    q=(Fraction(-1),0,Fraction(1),0,0,0)
    p=(Fraction(1,3),Fraction(-1,4),Fraction(1,5),0,0,0)
    p2=(0,Fraction(1,7),0,Fraction(1,6),0,0)
    result={}
    for name,r in [("nested_two_point",self_energy_insertion(k,q,p)),
                   ("nested_three_point",triangle_insertion(k,q,p,p2)),
                   ("overlapping_two_point",overlapping_self_energy(k,q,p))]:
        result[name]={"terms":[str(x) for x in r["terms"]],"sum":str(r["sum"]),
                      "operators":r["operators"],"integrated":False}
    return result
