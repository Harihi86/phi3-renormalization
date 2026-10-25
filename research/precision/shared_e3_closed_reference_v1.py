"""Original Gaussian Schur solve and optional direct integral of the entire shared auxiliary coordinate."""
from fractions import Fraction as Q
import mpmath as mp
import exact_rational_sector_transport_v1 as A
import shared_bubble_square_primitive_v1 as S
import shared_weighted_e2_scalar_v1 as B


def map_orbit(table,point):
    assert table['edges']==6 and len(point)==6
    assert all(type(x) in (int,Q) and 0<x<1 for x in point)
    rows=[]
    for orbit in range(6):
        order,shapes,probability=A.decode(table,(point[0]+Q(orbit,6))%1)
        alpha=[Q()]*6;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,shape) in enumerate(zip(point[1:],shapes),1):
            if shape<0:
                ratio=u**(-shape);derivative=(-shape)*u**(-shape-1)
            else:
                denominator=1+(shape-1)*u;ratio=shape*u/denominator;derivative=Q(shape)/denominator**2
            scale*=ratio;alpha[order[j]]=scale;jacobian*=derivative*ratio**(5-j)
        rows.append({'alpha':alpha,'order':order,'shapes':shapes,'probability':probability,'weight':jacobian/probability})
    return rows


def evaluate(point,binding,table,formulas,digits=640,quadrature=False):
    assert binding['weighted_kind']=='shared'
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda x:ctx.mpf(x.numerator)/x.denominator
    def unnormalized(t):
        w=4/(ctx.sqrt(t+4)+ctx.sqrt(t))**2
        return S.evaluate_w(ctx,formulas['0'],w)
    total=ctx.mpf(0)
    for row in map_orbit(table,point):
        alpha=list(map(number,row['alpha']));mass=ctx.fsum(alpha);subtotal=ctx.mpf(0)
        weights=binding['paired_states'][0]['weighted_edges'];assert len(weights)==1
        d=alpha[weights[0]];t0=(mass-d)/d
        common=None if quadrature else unnormalized(t0)
        for state in binding['paired_states']:
            assert state['weighted_edges']==weights;matrix=ctx.matrix(4)
            for i in range(4):
                for j in range(4):matrix[i,j]=ctx.fsum(a*r[i]*r[j] for a,r in zip(alpha,state['routes']))
            block=ctx.matrix([[matrix[i,j] for j in range(2)] for i in range(2)])
            cross=ctx.matrix([[matrix[i,2+j] for j in range(2)] for i in range(2)]);solved=ctx.matrix(2)
            for j in range(2):
                col=ctx.lu_solve(block,ctx.matrix([cross[i,j] for i in range(2)]))
                for i in range(2):solved[i,j]=col[i]
            schur=ctx.matrix([[matrix[2+i,2+j] for j in range(2)] for i in range(2)])-cross.T*solved
            U=ctx.det(block);z=schur[0,0]+schur[1,1]-(schur[0,1]+schur[1,0])/2
            assert U>0 and z>0;delta=z/d
            if quadrature:
                def f(v):
                    t=t0+delta*v
                    return -delta*B.hyper_bubble(ctx,t)**2/(1+t)**3
                value=ctx.quad(f,[0,ctx.mpf('.5'),1])
            else:value=common-unnormalized(t0+delta)
            subtotal+=state['sign']*value/U**3
        total+=number(row['weight'])*subtotal/6
    return ctx,total
