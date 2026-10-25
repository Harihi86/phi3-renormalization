"""Exhaustive fixed diagnostic partition; no numerical selection or evaluator."""
GRAPHS = ('L4E3G39', 'L4E3G394')
VARIANTS = {f'sector{power}-smooth{order}': (f'sector{power}-positive', order)
            for power in (2, 4) for order in (1, 2)}
SEEDS = [12410244001 + 239 * repeat for repeat in range(4)]


def fixed_plan():
    return [dict(id=gid, variant=variant, repeat=repeat, seed=seed,
                 mode=mode, order=order, nodes=512)
            for gid in GRAPHS for variant, (mode, order) in VARIANTS.items()
            for repeat, seed in enumerate(SEEDS)]


def key(row):
    return row['id'], row['variant'], row['repeat']


def validate(reused, missing):
    assert type(reused) is list and type(missing) is list and reused and missing
    expected = {key(row): row for row in fixed_plan()}
    seen = set()
    for origin, rows in [('reused', reused), ('missing', missing)]:
        for row in rows:
            identity = key(row)
            assert identity in expected and identity not in seen
            seen.add(identity)
            assert all(type(row[k]) is type(v) and row[k] == v
                       for k, v in expected[identity].items())
            if origin == 'reused':
                summary = row['verified_summary']
                assert all(type(summary[k]) is type(row[k]) and summary[k] == row[k]
                           for k in ('id', 'variant', 'repeat', 'seed', 'nodes'))
                assert set(summary['full']) == set(summary['half']) == {'midpoint', 'radius'}
                assert all(summary[name]['radius'] >= 0 for name in ('full', 'half'))
                assert len(row['native_receipts']) == 2
                assert len({receipt['folder'] for receipt in row['native_receipts']}) == 2
            else:
                assert row['reason'] == 'old group did not finish all512 nodes'
    assert seen == set(expected)
    return fixed_plan()
