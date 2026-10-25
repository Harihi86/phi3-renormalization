"""Independent nine-sector maps, original Gaussian Schur solve and T0 primitive."""
from fractions import Fraction as Q
import mpmath as mp
import exact_rational_sector_transport_v1 as A
import nested_analytic_marginals_v1 as M
import shared_weighted_e2_scalar_v1 as B


def map_orbit(table,point):
    assert table['edges']==9 and len(point)==9
    assert all(type(x) in (int,Q) and 0<x<1 for x in point)
    rows=[]
    for orbit in range(9):
        order,shapes,probability=A.decode(table,(point[0]+Q(orbit,9))%1)
        alpha=[Q()]*9;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,shape) in enumerate(zip(point[1:],shapes),1):
            if shape<0:
                ratio=u**(-shape);derivative=(-shape)*u**(-shape-1)
            else:
                denominator=1+(shape-1)*u
                ratio=shape*u/denominator;derivative=Q(shape)/denominator**2
            scale*=ratio;alpha[order[j]]=scale
            jacobian*=derivative*ratio**(8-j)
        rows.append({'alpha':alpha,'order':order,'shapes':shapes,
                     'probability':probability,'weight':jacobian/probability})
    return rows


def kernel(ctx,mass,edge,z,quadrature=False):
    assert mass>edge>0 and z>0
    t0=(mass-edge)/edge;delta=z/edge
    if quadrature:
        # An integral of the independent hypergeometric bubble, with no dilogarithms.
        def f(v):
            t=t0+delta*v
            return -delta*B.hyper_bubble(ctx,t)/(1+t)**2
        return ctx.quad(f,[0,ctx.mpf('.5'),1])
    return M.bubble_primitive(ctx,t0)-M.bubble_primitive(ctx,t0+delta)


def evaluate(point,binding,table,digits=640,quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator
    total=ctx.mpf(0)
    for row in map_orbit(table,point):
        alpha=list(map(number,row['alpha']));mass=ctx.fsum(alpha);subtotal=ctx.mpf(0)
        for state in binding['paired_states']:
            matrix=ctx.matrix(5)
            for i in range(5):
                for j in range(5):
                    matrix[i,j]=ctx.fsum(a*route[i]*route[j] for a,route in zip(alpha,state['routes']))
            block=ctx.matrix([[matrix[i,j] for j in range(3)] for i in range(3)])
            cross=ctx.matrix([[matrix[i,3+j] for j in range(2)] for i in range(3)])
            solved=ctx.matrix(3,2)
            for j in range(2):
                column=ctx.lu_solve(block,ctx.matrix([cross[i,j] for i in range(3)]))
                for i in range(3):solved[i,j]=column[i]
            schur=ctx.matrix([[matrix[3+i,3+j] for j in range(2)] for i in range(2)])-cross.T*solved
            U=ctx.det(block);z=schur[0,0]+schur[1,1]-(schur[0,1]+schur[1,0])/2
            assert U>0 and z>0
            subtotal+=state['sign']*kernel(ctx,mass,alpha[state['weighted_edge']],z,quadrature)/U**3
        total+=number(row['weight'])*subtotal/9
    return ctx,total
