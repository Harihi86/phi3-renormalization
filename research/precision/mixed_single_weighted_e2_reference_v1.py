"""Original effective Gaussian matrices, high-precision closed J or full remainder quadrature."""
import mpmath as mp
import weighted_bubble_e2_closed_v1 as J
import weighted_bubble_e2_integral_reference_v1 as X

def evaluate(points,binding,digits=640,chart='weighted-linear',quadrature=False):
    ctx=mp.mp.clone();ctx.dps=digits
    def number(x):return ctx.mpf(x.numerator)/x.denominator if hasattr(x,'numerator') else ctx.mpf(x)
    p=list(map(number,points));assert len(p)==7 and all(0<x<1 for x in p)
    assert chart in ('weighted-linear','weighted-square')
    square=chart=='weighted-square';values=[x*x if square else x for x in p]
    jacobian=ctx.fprod(2*x for x in p) if square else ctx.mpf(1);total=ctx.mpf(0)
    for pivot in range(8):
        alpha=values.copy();alpha.insert(pivot,ctx.mpf(1));mass=ctx.fsum(alpha)
        d=alpha[binding['weighted_edge']];t0=(mass-d)/d
        if not quadrature:a0=J.M.bubble_primitive(ctx,t0);c0=J.weighted_moment(ctx,t0)
        for state in binding['paired_states']:
            routes=state['routes'];matrix=ctx.matrix(4)
            for i in range(4):
                for j in range(4):matrix[i,j]=ctx.fsum(a*r[i]*r[j] for a,r in zip(alpha,routes))
            block=ctx.matrix([[matrix[i,j] for j in range(3)] for i in range(3)])
            U=ctx.det(block);external=ctx.matrix([matrix[i,3] for i in range(3)])
            # Independent Schur solve, not the positive minor list or F determinant.
            Q=matrix[3,3]-(external.T*ctx.lu_solve(block,external))[0]
            assert U>0 and Q>0;d=alpha[state['weighted_edge']]
            if quadrature:primitive=X.reference_J(ctx,mass,d,Q)
            else:
                t1=t0+Q/d
                primitive=d*(t1*(J.M.bubble_primitive(ctx,t1)-a0)-(J.weighted_moment(ctx,t1)-c0))
            total+=state['sign']*primitive/U**3
    return ctx,total*jacobian
