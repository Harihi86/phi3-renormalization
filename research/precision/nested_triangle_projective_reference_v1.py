"""Unpaired original-forest logarithms from independent mpmath Schur solves."""
import mpmath as mp
import exact_rational_sector_transport_v1 as A


def evaluate(point,binding,table,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda q:ctx.mpf(q.numerator)/q.denominator
    total=ctx.mpf(0)
    for chart in A.map_orbit(table,point):
        alpha=list(map(number,chart['alpha']));mass=ctx.fsum(alpha);subtotal=ctx.mpf(0)
        for state in binding['states']:
            routes=state['routes']
            matrix=ctx.matrix([[ctx.fsum(a*r[i]*r[j] for a,r in zip(alpha,routes)) for j in range(6)] for i in range(6)])
            loop=matrix[:4,:4];cross=matrix[:4,4:];solved=ctx.matrix(4,2)
            for j in range(2):
                column=ctx.lu_solve(loop,cross[:,j])
                for i in range(4):solved[i,j]=column[i]
            schur=matrix[4:,4:]-cross.T*solved
            u=ctx.det(loop);h=schur[0,0]+schur[1,1]-(schur[0,1]+schur[1,0])/2
            assert u>0 and h>0
            subtotal-=state['sign']*ctx.log1p(h/mass)/u**3
        total+=number(chart['weight'])*subtotal/12
    return ctx,total
