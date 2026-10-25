"""Independent unpaired original-forest reference via mpmath LU solves."""
import mpmath as mp


def value(binding,alpha,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator
    a=list(map(number,alpha));mass=ctx.fsum(a);pieces=[]
    for state in binding['states']:
        rows=state['routes']
        matrix=ctx.matrix([[ctx.fsum(ae*r[i]*r[j] for ae,r in zip(a,rows)) for j in range(5)] for i in range(5)])
        loop=matrix[:4,:4];rhs=matrix[:4,4];solved=ctx.lu_solve(loop,rhs)
        u=ctx.det(loop);schur=matrix[4,4]-ctx.fsum(matrix[4,k]*solved[k] for k in range(4))
        ratio=schur/mass
        # The original forest sum, with no paired coefficients or relative logs.
        pieces.append(state['sign']*mass*((1+ratio)*ctx.log1p(ratio)-ratio)/u**3)
    return ctx.fsum(pieces)
