"""Strict normalization of inline and separately stored independent tile lists."""
from pathlib import Path


def validate(job,cfg,summary,role,records,snapshot=None):
    job=Path(job);assert role in ('central','geometry')
    assert type(records) is list and type(summary['tiles']) is int
    count=cfg['snapshot_tiles']
    if type(count) is list:
        assert role=='geometry' and 'tile_records' in summary
        assert all(type(v) is str and Path(v).is_absolute() and Path(v).name=='summary.json' and Path(v).parents[1].name=='tiles' for v in count)
        assert len(set(count))==len(count)
        count=len(count)
    assert type(count) is int and summary['tiles']==len(records)==count
    assert summary['snapshot_only'] is True and cfg['whole_graph_certificates_added']==0
    seen=set()
    for row in records:
        assert row['side'] in ('aligned','opposite') and type(row['tile']) is list and len(row['tile'])==3
        assert all(type(v) is int for v in row['tile'])
        key=(row['side'],tuple(row['tile']));assert key not in seen;seen.add(key)
        name='-'.join(map(str,row['tile']));path=Path(row['file'])
        assert path in (job/row['side']/name/'summary.json',job/row['side']/(name+'.json'))
        assert role!='central' or path.name=='summary.json'
        assert isinstance(row['sha256'],str) and len(row['sha256'])==64 and all(c in '0123456789abcdef' for c in row['sha256'])
    if 'tile_records' in summary:
        assert records==summary['tile_records'] and snapshot is None
    else:
        assert role=='geometry' and summary['exact_arithmetic_backend']=='gmpy2.mpq'
        assert summary['native_arithmetic_used'] is True and summary['no_production_native_evaluator_calls'] is True
        assert type(snapshot) is list and len(snapshot)==len(records)
        expected=set()
        for target in snapshot:
            assert target['side'] in ('aligned','opposite')
            path=Path(target['path']);assert path.name=='summary.json' and path.parents[1].name=='tiles'
            assert isinstance(target['sha256'],str) and len(target['sha256'])==64 and all(c in '0123456789abcdef' for c in target['sha256'])
            key=(target['side'],str(job/target['side']/path.parent.name/'summary.json'))
            assert key not in expected;expected.add(key)
        assert {(row['side'],row['file']) for row in records}==expected
    return records


def read(job,cfg,summary,role,evidence):
    if 'tile_records' in summary:
        records=validate(job,cfg,summary,role,summary['tile_records'])
        if type(cfg['snapshot_tiles']) is list:
            actual=[]
            for row in records:
                saved=evidence.read(row['file'],row['sha256'])
                actual.append(saved['input_summary'])
                evidence.bind(saved['input_summary'],saved['sha256'])
            assert len(set(actual))==len(actual) and set(actual)==set(cfg['snapshot_tiles'])
        return records
    path=Path(job)/'snapshot.json'
    snapshot=evidence.read(path,cfg['sources'][str(path)])
    records=evidence.read(Path(job)/'tile-records.json')
    validate(job,cfg,summary,role,records,snapshot)
    targets={(s['side'],Path(s['path']).parent.name):s for s in snapshot}
    for row in records:
        target=targets[(row['side'],'-'.join(map(str,row['tile'])))]
        saved=evidence.read(row['file'],row['sha256'])
        assert saved['input_summary']==target['path'] and saved['input_sha256']==target['sha256']
        assert saved['side']==row['side'] and saved['tile']==row['tile']
    return records
