"""Exact type boundary to the frozen integer-series transcendental routines."""
from fractions import Fraction
import exact_exp2_enclosures as BASE_EXP
import exact_radial_node_bounds_v2 as BASE_BUB
import exact_range_enclosures_v2 as BASE_SQRT


def plain(q):
    return Fraction(int(q.numerator),int(q.denominator))


def exp2(q):return BASE_EXP.exp2(plain(q))
def log2():return BASE_EXP.log2()
def log_bounds(q):return BASE_BUB.log_bounds(plain(q))
def bubble_bounds(q):return BASE_BUB.bubble_bounds(plain(q))
def sqrt_bounds(q):return BASE_SQRT.sqrt_bounds(plain(q))


def clear_caches():
    BASE_EXP.exp2.cache_clear();BASE_BUB.log_bounds.cache_clear()
