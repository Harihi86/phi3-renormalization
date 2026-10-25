"""Unpaired LU reference on the complete probabilistic projective measure."""
import mpmath as mp
import exact_rational_sector_transport_v1 as A
import overlap_e2_projective_reference_v1 as U


def evaluate(point,binding,table,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits
    number=lambda q:ctx.mpf(q.numerator)/q.denominator
    pieces=[]
    for chart in A.map_orbit(table,point):
        pieces.append(number(chart['weight'])*ctx.mpf(U.value(binding,chart['alpha'],digits))/11)
    return ctx,ctx.fsum(pieces)
