"""Add exact replay of signed range enclosures; retain the previous guard."""
from theta_two_radial_majorant_audit import *


def range_enclosure(data,N,area,polynomial=None):
    assert data['finite'] and data['method']=='range' and area>0
    assert len(data['majorant_methods'])==len(data['propagator_majorants'])==len(data['coefficient_majorants'])==1
    method=data['majorant_methods'][0]
    if method=='original_propagators':bound=propagators(data['propagator_majorants'][0],N)
    else:
        assert method=='absolute_coefficients'
        bound=coefficients(data['coefficient_majorants'][0],N,polynomial)
    upper=interval(data['upper_density'])[1];assert upper>=bound
    lo,hi=interval(data['real']);assert lo<=-upper*area and hi>=upper*area
    return {'interval':(lo,hi),'range_error_upper':upper*area,'methods':[method]}
