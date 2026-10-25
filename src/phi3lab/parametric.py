"""Exact Euclidean Symanzik polynomials from edge-labelled multigraphs.

For unit propagator powers and a common positive mass squared, the convention is
U = sum_T prod_(e not in T) x_e and
F = m2 U sum_e x_e + sum_(spanning 2-forest A|B) p_A**2 prod_(e not in A|B) x_e.
Each unordered two-forest is counted once.  External momenta are all incoming;
the last (sorted) external label is minus the sum of the others.  This module
does not perform integration or apply any BPH subtraction operator.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import prod

from .graphs import components, connected


def _powers(size, edge_ids):
    selected = set(edge_ids)
    return tuple(int(i not in selected) for i in range(size))


def _determinant(matrix):
    """Exact Gaussian elimination; det of the empty matrix is one."""
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        value = a[i][i]
        result *= value
        for j in range(i + 1, len(a)):
            factor = a[j][i] / value
            for k in range(i + 1, len(a)):
                a[j][k] -= factor * a[i][k]
            a[j][i] = Fraction(0)
    return result


def matrix_tree_u(g, alpha):
    """Independent U evaluation using the weighted Laplacian minor.

    Self-loops contribute to prod(alpha), but not to the Laplacian. Parallel
    edge conductances add.  This algorithm does not enumerate spanning trees.
    """
    alpha = tuple(Fraction(x) for x in alpha)
    if len(alpha) != len(g.edges) or any(x <= 0 for x in alpha):
        raise ValueError("One positive Feynman parameter is required per edge")
    laplacian = [[Fraction(0) for _ in range(g.n)] for _ in range(g.n)]
    for x, (a, b) in zip(alpha, g.edges):
        if a == b:
            continue
        value = 1 / x
        laplacian[a][a] += value
        laplacian[b][b] += value
        laplacian[a][b] -= value
        laplacian[b][a] -= value
    minor = [row[:-1] for row in laplacian[:-1]]
    return prod(alpha, start=Fraction(1)) * _determinant(minor)


def _monomial(powers, alpha):
    return prod((x ** n for x, n in zip(alpha, powers)), start=Fraction(1))


def evaluate_symanzik(polys, alpha, external_vectors, m2=1):
    """Evaluate U,F exactly at positive parameters and Euclidean kinematics.

    external_vectors is a sequence in polys['external_basis'] order, i.e. E-1
    vectors; the last external momentum is dependent and must not be supplied.
    Fractions/integers preserve rational arithmetic. Positive m2 is required.
    These samples are values of polynomials, never values of loop integrals.
    """
    alpha = tuple(Fraction(x) for x in alpha)
    m2 = Fraction(m2)
    if len(alpha) != len(polys['variables']) or any(x <= 0 for x in alpha):
        raise ValueError("One positive Feynman parameter is required per edge")
    if m2 <= 0:
        raise ValueError("A positive common mass squared is required")
    vectors = tuple(tuple(Fraction(x) for x in v) for v in external_vectors)
    if len(vectors) != len(polys['external_basis']):
        raise ValueError("Supply exactly E-1 independent external vectors")
    if vectors and (not vectors[0] or any(len(v) != len(vectors[0]) for v in vectors)):
        raise ValueError("External vectors must have the same positive dimension")
    u = sum((term['coefficient'] * _monomial(term['powers'], alpha)
             for term in polys['U']), Fraction(0))
    f = m2 * sum((term['coefficient'] * _monomial(term['powers'], alpha)
                  for term in polys['F_mass']), Fraction(0))
    dimension = len(vectors[0]) if vectors else 0
    for term in polys['F_kinematic']:
        momentum = tuple(sum((c * v[j] for c, v in zip(term['momentum'], vectors)),
                             Fraction(0)) for j in range(dimension))
        square = sum((x * x for x in momentum), Fraction(0))
        f += square * _monomial(term['powers'], alpha)
    return u, f


def symanzik(g):
    """Return JSON-serializable sparse U,F and reproducible local checks.

    Powers refer to the original internal edge indices, including duplicate
    endpoints and self-loops.  F_mass multiplies common m2; each F_kinematic
    record multiplies the squared Euclidean momentum represented by its
    integer coefficient vector.  Graphs must be connected and have an edge.
    """
    if not connected(range(g.n), g.edges) or not g.edges:
        raise ValueError("Symanzik construction needs a connected graph with an edge")
    size = len(g.edges)
    vertices = range(g.n)
    u = []
    for retained in combinations(range(size), g.n - 1):
        # A connected graph with n-1 edges is a spanning tree, also for
        # multigraphs: a self-loop or a parallel cycle cannot pass this test.
        if connected(vertices, tuple(g.edges[i] for i in retained)):
            u.append({'coefficient': 1, 'powers': list(_powers(size, retained))})
    labels = [label for label, _ in g.external]
    independent = labels[:-1]
    label_vectors = {label: tuple(int(i == j) for i in range(len(independent)))
                     for j, label in enumerate(independent)}
    if labels:
        label_vectors[labels[-1]] = (-1,) * len(independent)
    f_kinematic = []
    two_forests_count = 0
    if g.n >= 2:
        for retained in combinations(range(size), g.n - 2):
            blocks = components(vertices, tuple(g.edges[i] for i in retained))
            if len(blocks) != 2:
                continue
            two_forests_count += 1
            # components() chooses the component containing vertex 0 first.
            # The complement would give the negative momentum and is not
            # enumerated a second time.
            block = blocks[0]
            momentum = tuple(sum(label_vectors[label][i]
                                 for label, v in g.external if v in block)
                             for i in range(len(independent)))
            # A consistent sign gives stable records under changes of which
            # side is named first; the squared momentum is unchanged.
            first = next((x for x in momentum if x), 0)
            if first < 0:
                momentum = tuple(-x for x in momentum)
            if any(momentum):
                f_kinematic.append({'powers': list(_powers(size, retained)),
                                    'momentum': list(momentum)})
    mass_terms = Counter()
    for term in u:
        for i in range(size):
            powers = list(term['powers'])
            powers[i] += 1
            mass_terms[tuple(powers)] += term['coefficient']
    f_mass = [{'coefficient': coefficient, 'powers': list(powers)}
              for powers, coefficient in sorted(mass_terms.items())]
    result = {
        'convention': 'Euclidean F = m2*U*sum(x) + sum_2forests(p_A^2*prod_complement(x))',
        'variables': [f'x{i}' for i in range(size)],
        'external_labels': labels,
        'external_basis': [f'p{label}' for label in independent],
        'dependent_external': f'p{labels[-1]} = -sum(independent)' if labels else None,
        'loops': g.loops,
        'U': u,
        'F_mass': f_mass,
        'F_kinematic': f_kinematic,
        'spanning_trees_count': len(u),
        'spanning_two_forests_count': two_forests_count,
    }
    samples = []
    external_vectors = [(Fraction(i + 1, 2), Fraction(1, i + 2))
                        for i in range(len(independent))]
    for alpha in ((Fraction(1),) * size,
                  tuple(Fraction(i + 2, i + 1) for i in range(size))):
        value_u, value_f = evaluate_symanzik(result, alpha, external_vectors)
        determinant_u = matrix_tree_u(g, alpha)
        samples.append({'alpha': [str(x) for x in alpha],
                        'external_vectors': [[str(x) for x in v] for v in external_vectors],
                        'm2': '1', 'U': str(value_u), 'F': str(value_f),
                        'matrix_tree_U': str(determinant_u),
                        'matrix_tree_matches': value_u == determinant_u,
                        'euclidean_positive': value_u > 0 and value_f > 0})
    result['checks'] = {
        'u_homogeneous': all(sum(t['powers']) == g.loops for t in u),
        'f_homogeneous': all(sum(t['powers']) == g.loops + 1
                             for t in f_mass + f_kinematic),
        'matrix_tree_identity_samples': all(s['matrix_tree_matches'] for s in samples),
        'euclidean_positive_samples': all(s['euclidean_positive'] for s in samples),
    }
    result['samples'] = samples
    result['limits'] = ('Exact graph polynomials; the determinant equality and positivity '
                        'are checked at the recorded rational samples. No integration, '
                        'UV subtraction, or general analytic convergence proof is performed.')
    return result
