"""Read-only, lossless ledger selection even when numbered cohorts finish out of order."""
from pathlib import Path
import json,hashlib,re


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'))


def inventory(summary):
    strict=summary['strict_analytic_ids'];rows=summary['empirical_integrals'];pending=summary['unresolved_ids']
    ids=[row['id'] for row in rows];universe=strict+ids+pending
    assert len(universe)==len(set(universe))==162,'duplicate or missing original graph'
    assert summary['strict_analytic_count']==len(strict)==24
    assert summary['empirical_dual_kernel_count']==len(rows)
    assert summary['distinct_resolved_original_ids']==len(strict)+len(rows)
    assert summary['unresolved_count']==len(pending)
    assert summary['complete162_objective_achieved']==(not pending)
    for row in rows:
        assert row['empirical_gate_passed'] is True and type(row['actual_audit_exit']) is int and row['actual_audit_exit']==0
        assert row['rigorous_quadrature_error_bound'] is False
    return {'universe':frozenset(universe),'strict':frozenset(strict),
            'strict_index':canonical(summary['strict_analytic_index']),
            'empirical':{row['id']:canonical(row) for row in rows}}


def select(records):
    assert records,'at least one actually completed ledger required'
    decoded=[(record,inventory(record['data'])) for record in records]
    reference=decoded[0][1]
    for _,item in decoded:
        for field in ('universe','strict','strict_index'):assert item[field]==reference[field],'ledger universe or strict evidence changed'
    candidates=[]
    for candidate,head in decoded:
        if all(all(gid in head['empirical'] and head['empirical'][gid]==row for gid,row in prior['empirical'].items()) for _,prior in decoded):
            candidates.append(candidate)
    assert candidates,'completed ledgers have incomparable or conflicting evidence; explicit reconciliation required'
    # A later cohort number is not a logical parent. Equal inventories are harmless;
    # otherwise the selected record must contain every prior accepted row verbatim.
    chosen=max(candidates,key=lambda r:(r['completed_mtime_ns'],r['number']))
    return {key:chosen[key] for key in ('index','summary','supervisor')}


def resolve(root,floor=59):
    root=Path(root);records=[]
    for folder in sorted((root/'jobs').glob('precision-outcomes-index-*')):
        match=re.fullmatch(r'precision-outcomes-index-(\d+)',folder.name)
        if not match or int(match[1])<floor or not folder.is_dir():continue
        number=int(match[1]);sup=root/'jobs'/('precision-outcomes-index-supervision-'+str(number))/'status.json'
        status=folder/'status.json';summary=folder/'summary.json'
        if sup.is_file():
            actual=json.loads(sup.read_text())
            if actual.get('state')=='terminal' and type(actual.get('exit_code')) is int and actual['exit_code']!=0:continue
        else:raise AssertionError('uncertain ledger attempt must be inspected: '+folder.name)
        assert actual.get('state')=='terminal','wait for the active ledger writer'
        assert type(actual.get('exit_code')) is int and actual['exit_code']==0 and actual['sources_unchanged'] is True
        worker=json.loads(status.read_text());data=json.loads(summary.read_text())
        assert worker['state']=='terminal' and type(worker['exit_code']) is int and worker['exit_code']==0 and data['source_unchanged'] is True
        records.append({'number':number,'index':folder.name,'summary':sha(summary),'supervisor':sha(sup),
                        'completed_mtime_ns':sup.stat().st_mtime_ns,'data':data})
    return select(records)
