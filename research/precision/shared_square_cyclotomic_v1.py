"""Real rational six-periodic series for the complete squared-bubble primitive.

The Hermite rational part and zero-pole logarithms are kept exactly. All
nonzero-pole Li1/Li2/Li3 contributions are combined before evaluation.
"""
from fractions import Fraction as Q

COS2=(2,1,-1,-2,-1,1)
SIN2ROOT3=(0,3,3,0,-3,-3)


def periodic(formula):
    rows=[]
    for k in range(3):
        poles=formula['simple_poles'][str(k)];a,b=poles['rho_plus']
        rows.append([poles['real_poles']['1']+(-1)**n*poles['real_poles']['-1']+
                     a*COS2[n]+b*SIN2ROOT3[n] for n in range(6)])
    return rows


def coefficient(rows,n):
    assert type(n) is int and n>=1
    a,b,c=(row[n%6] for row in rows)
    return [-a/n+b/n**2-2*c/n**3,-b/n+2*c/n**2,-c/n]


def number(ctx,q):return ctx.mpf(q.numerator)/q.denominator if isinstance(q,Q) else ctx.mpf(q)


def rational_part(ctx,formula,w):
    log=ctx.log(w);value=ctx.mpf(0)
    def polynomial(coeffs):
        total=ctx.mpf(0)
        for c in coeffs:total=total*w+number(ctx,c)
        return total
    for k in range(3):
        r=formula['rational_parts'][str(k)];zero=formula['simple_poles'][str(k)]['real_poles']['0']
        value+=polynomial(r['numerator'])/polynomial(r['denominator'])*log**k+number(ctx,zero)*log**(k+1)/(k+1)
    return value


def partial(ctx,formula,w,count):
    w=ctx.mpf(w);assert 0<w<1 and type(count) is int and count>=1
    rows=periodic(formula);log=ctx.log(w);value=ctx.mpf(0)
    for n in range(count,0,-1):
        a,b,c=map(lambda q:number(ctx,q),coefficient(rows,n))
        value=value*w+a+log*b+log**2*c
    return rational_part(ctx,formula,w)+w*value


def tail(ctx,rows,w,count):
    w=ctx.mpf(w);assert 0<w<1 and type(count) is int and count>=1
    m0,m1,m2=(number(ctx,max(map(abs,row))) for row in rows)
    n=ctx.mpf(count+1);l=abs(ctx.log(w))
    return w**(count+1)/(1-w)*(m0/n+m1*(l/n+1/n**2)+m2*(l*l/n+2*l/n**2+2/n**3))
