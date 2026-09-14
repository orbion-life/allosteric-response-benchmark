"""Create local assay inputs from an independently acquired, hash-pinned workbook."""
from pathlib import Path
import argparse, hashlib, shutil
import pandas as pd
P = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser()
ap.add_argument('--source-xlsx', type=Path, required=True)
args = ap.parse_args()
expected = '58abcbd1a08ad9fed8778e3b3e8672bbca84ac5a576f1b0530ddd9f5164fbce2'
assert hashlib.sha256(args.source_xlsx.read_bytes()).hexdigest() == expected
out = P / 'local-only'
out.mkdir(exist_ok=True)
target = out / 'supplementary-table-7.xlsx'
if args.source_xlsx.resolve() != target.resolve():
    shutil.copy2(args.source_xlsx, target)
df = pd.read_excel(target, sheet_name='TableS7')
assert len(df) == 3691
df.to_csv(out / 'faure-table7.csv', index=False)
print('Verified TableS7 and prepared local-only assay inputs. Do not redistribute these files.')
