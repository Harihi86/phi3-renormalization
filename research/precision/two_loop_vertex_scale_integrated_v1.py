"""The independently bound G35 parameters in the universal exact K3 formula."""
from fractions import Fraction as Q
import two_loop_vertex_rank_one_v1 as F
import triangle_vertex_scale_integrated_v1 as SCALAR
K=SCALAR.K
scale_integral=SCALAR.scale_integral


def parameters(point,permutation):
    assert len(point)==7
    record=F.chart(tuple(point)+(Q(1,2),),permutation)
    return {**{k:record[k] for k in ('Q','u','N','d')},'jacobian':record['stable_measure']}


def evaluate(data,mp):return SCALAR.evaluate(data,mp)
