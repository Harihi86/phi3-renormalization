"""Six exact sector charts and the original two-external-momentum Schur solve."""
from fractions import Fraction as Q
import mpmath as mp
import exact_rational_sector_transport_v1 as A
import weighted_e3_scalar_v1 as S


def map_orbit(table,point):
    assert table['edges']==6 and len(point)==7
    assert all(type(x) in (int,Q) and 0<x<1 for x in point)
    rows=[]
    for orbit in range(6):
        order,shapes,probability=A.decode(table,(point[0]+Q(orbit,6))%1)
        alpha=[Q()]*6;alpha[order[0]]=Q(1);scale=Q(1);jacobian=Q(1)
        for j,(u,shape) in enumerate(zip(point[1:6],shapes),1):
            if shape<0:
                ratio=u**(-shape);derivative=(-shape)*u**(-shape-1)
            else:
                denominator=1+(shape-1)*u
                ratio=shape*u/denominator;derivative=Q(shape)/denominator**2
            scale*=ratio;alpha[order[j]]=scale
            jacobian*=derivative*ratio**(5-j)
        rows.append({'alpha':alpha,'order':order,'shapes':shapes,
                     'probability':probability,'weight':jacobian/probability})
    return rows


def evaluate(point,binding,table,digits=640,quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda x:ctx.mpf(x.numerator)/x.denominator
    shared=binding['weighted_kind']=='shared';v=number(point[6])
    r=1-v*v;tau=(1-r*r)/4;spectral=ctx.mpf(1) if shared else v*r**4/6
    total=ctx.mpf(0)
    for row in map_orbit(table,point):
        alpha=list(map(number,row['alpha']));mass=ctx.fsum(alpha);subtotal=ctx.mpf(0)
        for state in binding['paired_states']:
            matrix=ctx.matrix(4)
            for i in range(4):
                for j in range(4):
                    matrix[i,j]=ctx.fsum(a*route[i]*route[j] for a,route in zip(alpha,state['routes']))
            block=ctx.matrix([[matrix[i,j] for j in range(2)] for i in range(2)])
            cross=ctx.matrix([[matrix[i,2+j] for j in range(2)] for i in range(2)])
            solved=ctx.matrix(2)
            for j in range(2):
                column=ctx.lu_solve(block,ctx.matrix([cross[i,j] for i in range(2)]))
                for i in range(2):solved[i,j]=column[i]
            schur=ctx.matrix([[matrix[2+i,2+j] for j in range(2)] for i in range(2)])-cross.T*solved
            U=ctx.det(block);z=schur[0,0]+schur[1,1]-(schur[0,1]+schur[1,0])/2
            assert U>0 and z>0
            weights=[alpha[i] for i in state['weighted_edges']]
            if shared:
                assert len(weights)==1;d=weights[0];a=mass+v*z
                bubble=S.H.direct_bubble if quadrature else S.H.hyper_bubble
                value=-d*d*z*bubble(ctx,(a-d)/d)**2/a**3
            else:
                assert len(weights)==2
                kernel=S.nine_term_reference if quadrature else S.distinct_kernel
                value=kernel(ctx,mass,*weights,z,tau)
            subtotal+=state['sign']*value/U**3
        total+=number(row['weight'])*subtotal*spectral/6
    return ctx,total
