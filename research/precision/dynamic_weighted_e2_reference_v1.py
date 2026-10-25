"""Independent original-route references selected by verified shared-edge structure."""
import double_weighted_e2_reference_v1 as HD
import shared_weighted_e2_reference_v1 as HS
J,X,S=HD.J,HD.X,HD.S

def evaluate(points,binding,digits=640,chart='dynamic-weighted-linear',quadrature=False):
    assert chart in ('dynamic-weighted-linear','dynamic-weighted-square')
    assert binding['parent_dependent_boundaries_retained']
    if binding['weighted_kind']=='shared':
        assert binding['id']=='L4E2G1'
        return HS.evaluate(points,binding,digits,chart.replace('dynamic-weighted','shared-weighted'),quadrature)
    assert binding['weighted_kind']=='double' and binding['id'] in ('L4E2G8','L4E2G17')
    return HD.evaluate(points,binding,digits,chart.replace('dynamic-weighted','double-weighted'),quadrature)
