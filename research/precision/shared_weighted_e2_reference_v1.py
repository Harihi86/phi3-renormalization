"""Independent Gaussian Schur solves and hypergeometric/Feynman bubble integrals."""
import mpmath as mp
import shared_weighted_e2_scalar_v1 as S

def evaluate(points,binding,digits=640,chart='shared-weighted-linear',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    p=list(map(number,points));assert len(p)==5 and all(0<v<1 for v in p)
    assert chart in ('shared-weighted-linear','shared-weighted-square')
    square=chart=='shared-weighted-square';values=[x*x if square else x for x in p[:4]]
    jac=ctx.fprod(2*x for x in p[:4]) if square else ctx.mpf(1)
    total=ctx.mpf(0);bubble=S.direct_bubble if quadrature else S.hyper_bubble
    for pivot in range(5):
        alpha=values.copy();alpha.insert(pivot,ctx.mpf(1));mass=ctx.fsum(alpha)
        for state in binding['paired_states']:
            routes=state['routes'];matrix=ctx.matrix(3)
            for i in range(3):
                for j in range(3):matrix[i,j]=ctx.fsum(a*v[i]*v[j] for a,v in zip(alpha,routes))
            block=ctx.matrix([[matrix[i,j] for j in range(2)] for i in range(2)])
            U=ctx.det(block);external=ctx.matrix([matrix[i,2] for i in range(2)])
            z=matrix[2,2]-(external.T*ctx.lu_solve(block,external))[0]
            assert U>0 and z>=0 and len(state['weighted_edges'])==1
            d=alpha[state['weighted_edges'][0]];a=mass+p[4]*z
            # This original-route evaluator independently assembles the scalar
            # density instead of calling the production scalar density helper.
            total+=state['sign']*d*d*z*z*(1-p[4])*bubble(ctx,(a-d)/d)**2/(a**3*U**3)
    return ctx,total*jac
