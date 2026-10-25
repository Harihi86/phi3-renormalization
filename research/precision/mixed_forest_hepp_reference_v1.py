"""Exact proposal coordinates, then independent original Gaussian Schur solves."""
from fractions import Fraction as Q
import math
import mpmath as mp
import mixed_forest_hepp_v1 as A
import mixed_single_weighted_e2_reference_v1 as B
J,X=B.J,B.X

def map_orbit(table,points):
    assert len(points)==8 and all(type(x) in (int,Q) and 0<x<1 for x in points)
    result=[]
    for orbit in range(8):
        selector=(points[0]+Q(orbit,8))%1;order,powers,probability=A.decode(table,selector)
        alpha=[Q()]*8;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,power) in enumerate(zip(points[1:],powers),1):
            scale*=u**power;alpha[order[j]]=scale;jacobian*=power*u**(power*(8-j)-1)
        result.append({'alpha':alpha,'order':order,'powers':powers,'probability':probability,'weight':jacobian/probability})
    return result

def evaluate(points,binding,table,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda x:ctx.mpf(x.numerator)/x.denominator
    total=ctx.mpf(0)
    for row in map_orbit(table,points):
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
        total+=number(row['weight'])*subtotal/8
    return ctx,total
