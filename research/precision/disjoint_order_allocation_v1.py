"""Exact minimum-node four-axis Gauss orders under a shared remainder budget.

This only selects orders. Actual central values, isolated Gauss rules and
complex majorants must still be evaluated and independently audited.
"""
from bisect import bisect_right
from fractions import Fraction as Q
from itertools import product
from math import prod
import gauss_cauchy_resolvent as RES
A=RES.A


def inputs(volume,moduli,ratios,orders,allowance):
    volume=A.bounded(Q(volume));allowance=A.bounded(Q(allowance))
    moduli=tuple(A.bounded(Q(x)) for x in moduli)
    ratios=tuple(A.bounded(Q(x)) for x in ratios);orders=tuple(orders)
    assert volume>0 and allowance>0 and len(moduli)==len(ratios)==4
    assert all(m>=0 for m in moduli) and all(0<q<1 for q in ratios)
    assert orders and tuple(sorted(set(orders)))==orders
    assert all(type(n) is int and 2<=n<=64 for n in orders)
    return volume,moduli,ratios,orders,allowance


def choose(volume,moduli,ratios,orders,allowance):
    volume,moduli,ratios,orders,allowance=inputs(volume,moduli,ratios,orders,allowance)
    limit=A.bounded(3*allowance/4)
    factors={(q,n):RES.factor(n,q) for q in set(ratios) for n in orders}
    errors=[{n:A.bounded(volume*m*factors[q,n]) for n in orders} for m,q in zip(moduli,ratios)]
    def pairs(i,j):
        return [(A.bounded(errors[i][a]+errors[j][b]),a*b,(a,b)) for a,b in product(orders,repeat=2)]
    right=sorted(pairs(2,3));right_errors=[r[0] for r in right];prefix=[];best=None
    for row in right:
        if best is None or (row[1],row[2])<(best[1],best[2]):best=row
        prefix.append(best)
    winner=None
    for left in pairs(0,1):
        available=A.bounded(limit-left[0]);index=bisect_right(right_errors,available)-1
        if index<0:continue
        partner=prefix[index];candidate=(left[1]*partner[1],left[2]+partner[2])
        if winner is None or candidate<winner:winner=candidate
    if winner is None:return None
    nodes,chosen=winner;error,parts=RES.error(volume,moduli,ratios,chosen)
    assert error<=limit and nodes==prod(chosen)
    return {'orders':chosen,'tensor_nodes':nodes,'axis_errors':parts,'remainder_bound':error,
            'remainder_limit':limit,'cell_allowance':allowance,
            'reserved_for_central_enclosure':A.bounded(allowance-error),
            'requires_new_central_evaluation_and_full_cell_audit':True}


def audit_choice(result,volume,moduli,ratios,orders,allowance):
    volume,moduli,ratios,orders,allowance=inputs(volume,moduli,ratios,orders,allowance)
    chosen=tuple(result['orders']);assert len(chosen)==4 and all(n in orders for n in chosen)
    error,parts=RES.error(volume,moduli,ratios,chosen);limit=A.bounded(3*allowance/4)
    assert result['tensor_nodes']==prod(chosen)
    assert tuple(result['axis_errors'])==tuple(parts) and result['remainder_bound']==error<=limit
    assert result['remainder_limit']==limit and result['cell_allowance']==allowance
    assert result['reserved_for_central_enclosure']==allowance-error>=allowance/4
    assert result['requires_new_central_evaluation_and_full_cell_audit'] is True
    return error
