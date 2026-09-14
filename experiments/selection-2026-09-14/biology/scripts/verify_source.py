"""Independent workbook-to-derived-site audit; external table is optional for public replays."""
from pathlib import Path
import json,sys,hashlib,math,collections,argparse
from decimal import Decimal
import openpyxl
P=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--source-xlsx',type=Path,default=P/'local-only/supplementary-table-7.xlsx');args=ap.parse_args()
expected='58abcbd1a08ad9fed8778e3b3e8672bbca84ac5a576f1b0530ddd9f5164fbce2';assert hashlib.sha256(args.source_xlsx.read_bytes()).hexdigest()==expected
w=openpyxl.load_workbook(args.source_xlsx,read_only=True,data_only=True);it=iter(w['TableS7'].values);head=next(it);sums=collections.defaultdict(lambda:[0,Decimal(0),Decimal(0)]);valid=0
for row in it:
 d=dict(zip(head,row))
 if d['protein']!='GRB2-SH3' or d['mut_order']!=1:continue
 vals=[]
 try:vals=[float(d[k]) for k in ['b_ddg_pred','f_ddg_pred','b_ddg_pred_sd','f_ddg_pred_sd']]
 except (ValueError,TypeError):continue
 if not all(math.isfinite(x) for x in vals) or vals[2]<=0 or vals[3]<=0:continue
 key=int(d['Pos'])+158;sums[key][0]+=1;sums[key][1]+=abs(Decimal(str(vals[0])));sums[key][2]+=abs(Decimal(str(vals[1])));valid+=1
import csv
with (P/'results/site-aggregates-and-predictions.csv').open() as f:derived=list(csv.DictReader(f))
errors=[]
for r in derived:
 key=int(r['canonical']);n,b,fo=sums.get(key,[0,0,0])
 if n:
  assert n==int(float(r['mutation_count']));errors.extend([abs(float(b/n)-float(r['binding_magnitude'])),abs(float(fo/n)-float(r['folding_magnitude']))])
assert valid==756 and max(errors)<1e-14
print(json.dumps({'status':'PASS','source_xlsx_sha256':expected,'valid_common_substitutions':valid,'independent_Decimal_site_aggregate_max_abs_error':max(errors),'scope':'Workbook parsed independently of production CSV and pandas aggregation. Publisher workbook is local-only and not packaged for redistribution.'},indent=2))
