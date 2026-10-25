"""Convention-aware leading bare dimensional-regularization pole benchmarks.

J_G(d) = integral d^d k / pi^(d/2) prod_e 1/(q_e^2 + m2), d=D-2*eps.
There are no diagram symmetry factors, couplings, mu factors, or BPH operators
in J.  A positive common mass makes the one-loop Euclidean parameter integral
infrared finite.  Its simple UV pole comes from Gamma(I-d/2).
"""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod

from .graphs import generate, validate_scope
from .parametric import evaluate_symanzik, symanzik


def simplex_moment(powers):
    """Integral over x_i>=0, sum x_i=1 of prod x_i**powers_i.

    Uses delta(1-sum x_i) prod dx_i; its unit-volume integral is 1/(N-1)!,
    not one. This convention is essential for the triangle's residue 1/2.
    """
    powers = tuple(powers)
    if not powers or any(type(n) is not int or n < 0 for n in powers):
        raise ValueError("A nonempty sequence of nonnegative integer powers is required")
    return Fraction(prod(factorial(n) for n in powers),
                    factorial(sum(powers) + len(powers) - 1))


def _multiply(a, b):
    result = defaultdict(Fraction)
    for powers_a, value_a in a.items():
        for powers_b, value_b in b.items():
            result[tuple(x + y for x, y in zip(powers_a, powers_b))] += value_a * value_b
    return dict(result)


def _euclidean_f(polys, m2, vectors):
    result = defaultdict(Fraction)
    for term in polys['F_mass']:
        result[tuple(term['powers'])] += m2 * term['coefficient']
    dimension = len(vectors[0]) if vectors else 0
    for term in polys['F_kinematic']:
        q = tuple(sum((c * v[j] for c, v in zip(term['momentum'], vectors)), Fraction(0))
                  for j in range(dimension))
        result[tuple(term['powers'])] += sum((x * x for x in q), Fraction(0))
    return dict(result)


def one_loop_pole(g, dim=6, m2=1, external_vectors=None):
    """Return an exact JSON-ready 1/eps coefficient for one bare one-loop graph.

    Scope is the atlas's bare trivalent connected 1PI one-loop E=1,2,3 graphs,
    with D=4 or 6 and d=D-2*eps. external_vectors supplies E-1 independent
    Euclidean momenta in sorted-label order; the last external leg is minus
    their sum. If omitted, fixed rational sample momenta are recorded.

    For n=D/2-I>=0, residue=(-1)^n/n! * integral_simplex F^n. For n<0,
    the massive one-loop integral has no UV pole. There are no finite parts.
    """
    validate_scope(g)
    if g.loops != 1 or type(dim) is not int or dim not in (4, 6):
        raise ValueError("Pole benchmarks support one-loop atlas graphs in D=4 or D=6")
    polys = symanzik(g)
    if external_vectors is None:
        external_vectors = [(Fraction(i + 1, 2), Fraction(1, i + 2))
                            for i in range(len(polys['external_basis']))]
    vectors = tuple(tuple(Fraction(x) for x in v) for v in external_vectors)
    m2 = Fraction(m2)
    # Reuse the public parameter evaluator's kinematic and positive-mass checks.
    evaluate_symanzik(polys, [1] * len(g.edges), vectors, m2)
    n = dim // 2 - len(g.edges)
    if n >= 0:
        polynomial = {(0,) * len(g.edges): Fraction(1)}
        f = _euclidean_f(polys, m2, vectors)
        for _ in range(n):
            polynomial = _multiply(polynomial, f)
        moment = sum((value * simplex_moment(powers)
                      for powers, value in polynomial.items()), Fraction(0))
        gamma_factor = Fraction((-1) ** n, factorial(n))
        residue = gamma_factor * moment
    else:
        polynomial = {}
        moment = None
        gamma_factor = Fraction(0)
        residue = Fraction(0)
    legs = len(g.external)
    square = sum((x * x for x in vectors[0]), Fraction(0)) if vectors else Fraction(0)
    if dim == 6:
        known = {1: m2 ** 2 / 2, 2: -m2 - square / 6, 3: Fraction(1, 2)}[legs]
        expression = {1: 'm2^2/2',
                      2: f'-m2-{polys["external_basis"][0]}_E^2/6' if legs == 2 else '',
                      3: '1/2'}[legs]
    else:
        known = {1: -m2, 2: Fraction(1), 3: Fraction(0)}[legs]
        expression = {1: '-m2', 2: '1', 3: '0 (UV finite)'}[legs]
    return {
        'dimension': f'{dim}-2*eps',
        'loops': 1,
        'external_legs': legs,
        'm2': str(m2),
        'external_basis': polys['external_basis'],
        'external_vectors': [[str(x) for x in v] for v in vectors],
        'gamma_index_n': n,
        'gamma_simple_pole_factor': str(gamma_factor),
        'simplex_f_power': n if n >= 0 else None,
        'simplex_integral': str(moment) if moment is not None else None,
        'integrated_polynomial': [{'powers': list(p), 'coefficient': str(c)}
                                  for p, c in sorted(polynomial.items())],
        'simple_pole_coefficient': str(residue),
        'known_closed_form': expression,
        'known_closed_form_value': str(known),
        'checks': {'simplex_matches_closed_form': residue == known},
        'normalization': ('Euclidean J = integral d^d k/pi^(d/2) product_e '
                          '1/(q_e^2+m2), d=D-2*eps; unit propagator powers; '
                          'no 1/S, coupling, mu factor, or BPH subtraction'),
        'limits': ('Leading bare DR 1/eps coefficient only. This is neither a '
                   'finite integral evaluation nor a BPH-renormalized amplitude; '
                   'it is not directly comparable to a different loop-measure '
                   'or d=D-eps convention without converting factors.'),
    }


def one_loop_pole_benchmarks(dim=6):
    """Three reproducible tadpole/bubble/triangle normalization benchmarks."""
    return [{'id': f'L1E{legs}G1',
             **one_loop_pole(generate(1, legs)[0], dim=dim)}
            for legs in (1, 2, 3)]
