"""Directed integer Horner bounds, independent of native Arb polynomial code."""
from pathlib import Path
from fractions import Fraction as Q
import theta_two_radial_majorant_audit_v4 as A

BITS=1536
SCALE=1<<BITS


def scaled(q):
    q=A.bounded(Q(q));assert q>=0
    n=q.numerator*SCALE;d=q.denominator
    return n//d,-(-n//d)


def mul(a,b):
    assert max(x.bit_length() for x in a)+max(x.bit_length() for x in b)<=A.CAP
    return a[0]*b[0]//SCALE,-(-(a[1]*b[1])//SCALE)


class Polynomial:
    def __init__(self,path):
        with Path(path).open() as stream:
            self.order,count=map(int,stream.readline().split());assert 1<=self.order<=22 and count<=100000
            terms={}
            for line in stream:
                i,j,n,d=map(int,line.split());assert 0<=i and 0<=j and i+j<=12*self.order+16 and d>0
                assert (i,j) not in terms;terms[i,j]=scaled(abs(Q(n,d)))
        assert len(terms)==count;self.count=count
        self.rows=[]
        for i in range(max(i for i,j in terms)+1):
            js=[j for ii,j in terms if ii==i]
            self.rows.append([terms.get((i,j),(0,0)) for j in range(max(js)+1)] if js else [])

    def verify(self,proof):
        p={k:A.interval(v) for k,v in proof.items()}
        r,w=scaled(p['r_mod'][1]),scaled(p['w_mod'][1]);total=(0,0)
        for row in reversed(self.rows):
            inner=(0,0)
            for coefficient in reversed(row):
                inner=mul(inner,w);inner=inner[0]+coefficient[0],inner[1]+coefficient[1]
            total=mul(total,r);total=total[0]+inner[0],total[1]+inner[1]
        upper=A.bounded(Q(total[1],SCALE))
        assert p['numerator_absolute'][1]>=upper,'Independent positive polynomial upper bound'
        return {'all_coefficients_replayed':self.count,'directed_integer_bits':BITS,'integer_cap':A.CAP}
