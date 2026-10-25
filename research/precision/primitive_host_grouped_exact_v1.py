"""Rational conditional decoder; no floating point enters sector selection."""
from fractions import Fraction as Q


def decode(value, table, proposal, stratum):
    value = Q(value)
    assert 0 <= value < 1
    gamma = {int(k): v for k, v in table['gamma'].items()}
    nu = {int(k): v for k, v in table['nu'].items()}
    weights = proposal['weights']
    edge = proposal['edge']
    position, first = stratum['position'], stratum['first']
    remaining, order, exponents, leading, probability = 63, [], [], [], Q(1)
    for depth in range(6):
        k = position-depth if position != -1 and remaining & (1 << edge) else -1
        options = [e for e in range(6) if remaining & (1 << e)
                   and (k == -1 or ((e == edge) == (k == 0)))
                   and (depth != 0 or first == -1 or e == first)]
        children = [remaining ^ (1 << e) for e in options]
        masses = [weights[f'{s}:{-1 if k == -1 or e == edge else k-1}'] / gamma.get(s, 1)
                  for e, s in zip(options, children)]
        total = sum(masses, Q())
        left = Q()
        for e, child, mass in zip(options, children, masses):
            p = mass/total
            if value < left+p:
                value = (value-left)/p
                probability *= p
                order.append(e)
                remaining = child
                if depth < 5:
                    exponents.append(gamma[child]); leading.append(nu[child])
                break
            left += p
        else:
            raise AssertionError('conditional CDF did not cover the point')
    assert remaining == 0 and sorted(order) == list(range(6))
    assert first == -1 or order[0] == first
    assert position == -1 or order[position] == edge
    return {'order': order, 'exponents': exponents, 'leading': leading,
            'auxiliary': value, 'probability': probability}
