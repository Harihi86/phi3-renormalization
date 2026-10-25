"""Exact propagation audit for the documented FLINT certified integral API."""
from fractions import Fraction as Q
from fourier_spectral_audit import interval,exp_negative64
def verify(data):
    assert data['success'] and data['calc_exit']==0 and data['bits']>=128
    assert data['calls']>=data['analytic_calls']>0
    rl,rh=interval(data['r']);assert Q(1,262144)<=rl<=rh<=12
    ql,qh=interval(data['quadrature']);il,ih=interval(data['quadrature_imag']);assert il<=0<=ih
    el,eh=exp_negative64();tl,th=interval(data['spectral_tail']);assert tl<=Q(8,27)*396*el and th>=Q(8,27)*396*eh
    cl,ch=interval(data['C_scaled']);assert cl<=max(Q(),ql) and ch>=qh+th
    assert (ch-cl)/2<(1+abs((ch+cl)/2))/10**12
    return {'r_interval':(rl,rh),'C_scaled':(cl,ch),'quadrature':(ql,qh),'spectral_tail_upper':th,'calls':data['calls'],'analytic_calls':data['analytic_calls'],'FLINT_integrator_is_bound_authority':True}
