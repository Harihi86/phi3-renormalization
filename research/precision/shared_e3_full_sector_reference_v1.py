"""Unreduced720-order reference, exact Gaussian elimination and independent primitive evaluation."""
from fractions import Fraction as Q
from itertools import permutations
import mpmath as mp
import shared_bubble_square_primitive_v1 as S


def solve(matrix,rhs):
    n=len(rhs);a=[list(row)+[rhs[i]] for i,row in enumerate(matrix)]
    for i in range(n):
        pivot=next(k for k in range(i,n) if a[k][i]);a[i],a[pivot]=a[pivot],a[i]
        value=a[i][i];a[i]=[v/value for v in a[i]]
        for j in range(n):
            if j!=i:
                scale=a[j][i];a[j]=[a[j][k]-scale*a[i][k] for k in range(n+1)]
    return [a[i][-1] for i in range(n)]


def schur(state,alpha):
    matrix=[[sum((a*route[i]*route[j] for a,route in zip(alpha,state['routes'])),Q()) for j in range(4)] for i in range(4)]
    block=[row[:2] for row in matrix[:2]];cross=[row[2:] for row in matrix[:2]]
    solved=[solve(block,[cross[i][j] for i in range(2)]) for j in range(2)]
    remainder=[[matrix[2+i][2+j]-sum((cross[k][i]*solved[j][k] for k in range(2)),Q()) for j in range(2)] for i in range(2)]
    U=block[0][0]*block[1][1]-block[0][1]*block[1][0]
    z=remainder[0][0]+remainder[1][1]-(remainder[0][1]+remainder[1][0])/2
    assert U>0 and z>0;return U,z


def coefficients(point,binding,table):
    assert len(point)==5 and all(type(v) is Q and 0<v<1 for v in point)
    assert binding['weighted_kind']=='shared';out={};weights=binding['paired_states'][0]['weighted_edges']
    def add(t,c):out[t]=out.get(t,Q())+c
    count=0
    for order in permutations(range(6)):
        remaining=set(range(6));remaining.remove(order[0]);alpha=[Q()]*6;alpha[order[0]]=Q(1);scale=Q(1);jac=Q(1)
        for j,u in enumerate(point,1):
            mask=sum(2**edge for edge in remaining);shape=table['shapes'][mask]
            if shape<0:
                r=u**(-shape);derivative=(-shape)*u**(-shape-1)
            else:
                r=shape*u/(1-u+shape*u);derivative=Q(shape)/(1-u+shape*u)**2
            scale*=r;alpha[order[j]]=scale;jac*=derivative*r**(5-j);remaining.remove(order[j])
        assert not remaining;mass=sum(alpha,Q());d=alpha[weights[0]];t0=(mass-d)/d
        for state in binding['paired_states']:
            assert state['weighted_edges']==weights;U,z=schur(state,alpha);c=state['sign']*jac/U**3
            add(t0,c);add(t0+z/d,-c)
        count+=1
    out={t:c for t,c in out.items() if c};assert count==720 and sum(out.values(),Q())==0 and all(t>0 for t in out)
    return out


def evaluate_coefficients(exact,formulas,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda q:ctx.mpf(q.numerator)/q.denominator;values=[]
    assert sum(exact.values(),Q())==0
    for t,c in sorted(exact.items()):
        value=number(t);w=4/(ctx.sqrt(value+4)+ctx.sqrt(value))**2
        values.append(number(c)*S.evaluate_w(ctx,formulas['0'],w))
    return ctx,ctx.fsum(values)


def evaluate(point,binding,table,formulas,digits=640):
    return evaluate_coefficients(coefficients(point,binding,table),formulas,digits)
