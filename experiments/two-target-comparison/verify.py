#!/usr/bin/env python3
"""Verify Holm boundary/tie cases, exact replay and the recorded audit families."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib,json,math,subprocess,sys
from audit import P,holm,sha
fixtures=[([.01,.04,.03],[.03,.06,.06]),([.02,.02,.9],[.06,.06,.9]),([0.,1.,1.],[0.,1.,1.]),([1.,1.],[1.,1.])]
for values,expected in fixtures:
    assert all(math.isclose(a,b,abs_tol=1e-15,rel_tol=0) for a,b in zip(holm(values),expected))
try:holm([float('nan')]);raise AssertionError('Non-finite p was accepted.')
except ValueError:pass
with TemporaryDirectory() as tmp:
    subprocess.run([sys.executable,str(P/'audit.py'),'--output',tmp],check=True,stdout=subprocess.DEVNULL)
    for name in ['audit.json','primary.csv','paired.csv']:assert (P/'results'/name).read_bytes()==(Path(tmp)/name).read_bytes()
r=json.loads((P/'results/audit.json').read_text());assert len(r['primary'])==3 and len(r['paired'])==6
for family in [r['primary'],r['paired']]:
    for row in family:
        assert 0<=row['raw_p']<=row['holm_adjusted_p']<=1
        if row['target']=='MYH7':
            assert row['raw_p']==row['holm_adjusted_p']==1 and row['status']=='not_run_conservative_correction_slot'
            assert row.get('effect',row.get('mean_primary_minus_control_percentile')) is None
assert not any(x['below_0p05_in_this_retrospective_family'] for x in r['paired'])
kr=next(x for x in r['primary'] if x['target']=='KRAS');assert math.isclose(kr['holm_adjusted_p'],3/10001,abs_tol=1e-15)
ah=next(x for x in r['paired'] if x['target']=='ABL' and x['comparator']=='harmonic');assert math.isclose(ah['holm_adjusted_p'],6*ah['raw_p'],abs_tol=1e-15)
receipt={'status':'passed','holm_fixtures':len(fixtures),'nonfinite_input_rejected':True,'three_output_files_identical_on_fresh_process_replay':True,'unrun_slots_have_no_fabricated_effect':True,'paired_family_has_no_adjusted_p_at_or_below_0p05':True,'code_sha256':sha(P/'audit.py'),'verification_code_sha256':sha(__file__),'results_sha256':{n:sha(P/'results'/n) for n in ['audit.json','primary.csv','paired.csv']},'scope':'Arithmetic, family construction and faithful input selection. This does not validate exchangeability or biological efficacy.'}
(P/'verification-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
