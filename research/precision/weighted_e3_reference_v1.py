"""Original routed four-by-four Gaussian Schur reference for complete E3 forests."""
import mpmath as mp
import weighted_e3_scalar_v1 as S

def evaluate(points,binding,digits=640,chart='weighted-e3-linear',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    p=list(map(number,points));assert len(p)==6 and all(0<v<1 for v in p)
    assert chart in ('weighted-e3-linear','weighted-e3-square')
    square=chart=='weighted-e3-square';values=[x*x if square else x for x in p[:5]]
    jac=ctx.fprod(2*x for x in p[:5]) if square else ctx.mpf(1);shared=binding['weighted_kind']=='shared'
    r=1-p[5]**2;tau=(1-r*r)/4;weight=ctx.mpf(1) if shared else p[5]*r**4/6;total=ctx.mpf(0)
    for pivot in range(6):
        alpha=values.copy();alpha.insert(pivot,ctx.mpf(1));mass=ctx.fsum(alpha)
        for state in binding['paired_states']:
            matrix=ctx.matrix(4)
            for i in range(4):
                for j in range(4):matrix[i,j]=ctx.fsum(a*v[i]*v[j] for a,v in zip(alpha,state['routes']))
            block=ctx.matrix([[matrix[i,j] for j in range(2)] for i in range(2)])
            cross=ctx.matrix([[matrix[i,2+j] for j in range(2)] for i in range(2)])
            solved=ctx.matrix(2)
            for j in range(2):
                vector=ctx.lu_solve(block,ctx.matrix([cross[i,j] for i in range(2)]))
                for i in range(2):solved[i,j]=vector[i]
            schur=ctx.matrix([[matrix[2+i,2+j] for j in range(2)] for i in range(2)])-cross.T*solved
            U=ctx.det(block);z=schur[0,0]+schur[1,1]-(schur[0,1]+schur[1,0])/2
            assert U>0 and z>=0
            weights=[alpha[i] for i in state['weighted_edges']]
            if shared:
                assert len(weights)==1;d=weights[0];a=mass+p[5]*z
                bubble=S.H.direct_bubble if quadrature else S.H.hyper_bubble
                value=-d*d*z*bubble(ctx,(a-d)/d)**2/a**3
            else:
                assert len(weights)==2;kernel=S.nine_term_reference if quadrature else S.distinct_kernel
                value=kernel(ctx,mass,*weights,z,tau)
            total+=state['sign']*value/U**3
    return ctx,total*jac*weight
