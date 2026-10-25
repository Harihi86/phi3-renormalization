"""Exact proposal coordinates, then independent original Gaussian Schur solves."""
from fractions import Fraction as Q
import math
import mpmath as mp
import exact_prefix_sector_transport_v1 as A
import mixed_single_weighted_e2_reference_v1 as B
J,X=B.J,B.X

def map_orbit(table,point,depth):
    return A.mapped_rows(table,point,depth)


def evaluate(points,binding,table,digits=640,depth=1):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda x:ctx.mpf(x.numerator)/x.denominator
    total=ctx.mpf(0)
    for row in map_orbit(table,points,depth):
        alpha=list(map(number,row['alpha']));mass=ctx.fsum(alpha);d=alpha[binding['weighted_edge']];t0=(mass-d)/d
        a0=J.M.bubble_primitive(ctx,t0);c0=J.weighted_moment(ctx,t0);subtotal=ctx.mpf(0)
        for state in binding['paired_states']:
            routes=state['routes'];matrix=ctx.matrix(4)
            for i in range(4):
                for j in range(4):matrix[i,j]=ctx.fsum(a*r[i]*r[j] for a,r in zip(alpha,routes))
            block=ctx.matrix([[matrix[i,j] for j in range(3)] for i in range(3)])
            U=ctx.det(block);external=ctx.matrix([matrix[i,3] for i in range(3)])
            q=matrix[3,3]-(external.T*ctx.lu_solve(block,external))[0];assert U>0 and q>0
            t1=t0+q/d;primitive=d*(t1*(J.M.bubble_primitive(ctx,t1)-a0)-(J.weighted_moment(ctx,t1)-c0))
            subtotal+=state['sign']*primitive/U**3
        total+=number(row['weight'])*subtotal
    return ctx,total
