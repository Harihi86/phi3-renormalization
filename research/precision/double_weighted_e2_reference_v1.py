"""Original effective Gaussian Schur solves and independent full nine-term integrals."""
import mpmath as mp
import double_weighted_e2_scalar_v1 as S
J=S.F
X=S.B

def evaluate(points,binding,digits=640,chart='double-weighted-linear',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    p=list(map(number,points));assert len(p)==5 and all(0<v<1 for v in p)
    assert chart in ('double-weighted-linear','double-weighted-square')
    square=chart=='double-weighted-square';values=[x*x if square else x for x in p[:4]]
    jac=ctx.fprod(2*x for x in p[:4]) if square else ctx.mpf(1)
    r=1-p[4]**2;tau=(1-r*r)/4;weight=p[4]*r**4/6;total=ctx.mpf(0)
    for pivot in range(5):
        alpha=values.copy();alpha.insert(pivot,ctx.mpf(1));mass=ctx.fsum(alpha)
        for state in binding['paired_states']:
            routes=state['routes'];matrix=ctx.matrix(3)
            for i in range(3):
                for j in range(3):matrix[i,j]=ctx.fsum(a*v[i]*v[j] for a,v in zip(alpha,routes))
            block=ctx.matrix([[matrix[i,j] for j in range(2)] for i in range(2)])
            U=ctx.det(block);external=ctx.matrix([matrix[i,2] for i in range(2)])
            z=matrix[2,2]-(external.T*ctx.lu_solve(block,external))[0]
            assert U>0 and z>=0;d,e=[alpha[i] for i in state['weighted_edges']]
            kernel=S.nine_term_reference if quadrature else S.spectral_kernel
            total+=state['sign']*kernel(ctx,mass,d,e,z,tau)/U**3
    return ctx,total*jac*weight
