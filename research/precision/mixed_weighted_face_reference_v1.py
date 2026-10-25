"""Uncombined original Schur reference with independently expanded primitive A."""
import mpmath as mp
import exact_rational_sector_transport_v1 as A
import bubble_primitive_local_series_reference_v1 as SMALL
import bubble_primitive_large_series_reference_v2 as LARGE
import weighted_bubble_e2_closed_v1 as J
import weighted_bubble_e2_integral_reference_v1 as QUAD


def primitive(ctx,t):
    return (SMALL if t<16 else LARGE).evaluate(ctx,t)[0]


def value(ctx,binding,alpha):
    """No exact-coefficient collector, positive minors, or face matching used."""
    mass=ctx.fsum(alpha);d=alpha[binding['weighted_edge']];t0=(mass-d)/d
    a0=primitive(ctx,t0);c0=J.weighted_moment(ctx,t0);values=[]
    for state in binding['paired_states']:
        routes=state['routes'];matrix=ctx.matrix(4)
        for i in range(4):
            for j in range(4):matrix[i,j]=ctx.fsum(a*r[i]*r[j] for a,r in zip(alpha,routes))
        block=ctx.matrix([[matrix[i,j] for j in range(3)] for i in range(3)])
        u=ctx.det(block);external=ctx.matrix([matrix[i,3] for i in range(3)])
        h=matrix[3,3]-(external.T*ctx.lu_solve(block,external))[0]
        assert u>0 and h>0
        t1=t0+h/d
        j=d*(t1*(primitive(ctx,t1)-a0)-(J.weighted_moment(ctx,t1)-c0))
        values.append(state['sign']*j/u**3)
    return ctx.fsum(values)


def evaluate(point,binding,table,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda q:ctx.mpf(q.numerator)/q.denominator
    return ctx,ctx.fsum(number(row['weight'])*value(ctx,binding,list(map(number,row['alpha'])))/8
                       for row in A.map_orbit(table,point))
