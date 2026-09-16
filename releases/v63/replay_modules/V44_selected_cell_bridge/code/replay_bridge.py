"""Portable numerical replay; never overwrites the released reference results."""
from pathlib import Path
import argparse,shutil,sys,json
import pandas as pd,numpy as np
import audit_inputs
parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
source=Path(__file__).resolve().parents[1];dest=Path(args.out).resolve()
assert not dest.exists(),'Choose a new output directory'
dest.mkdir(parents=True)
for d in ['inputs','audit']:shutil.copytree(source/d,dest/d)
for d in ['results','qa','figures']: (dest/d).mkdir()
for f in ['FROZEN_PLAN.md','FREEZE_RECEIPT.json']:shutil.copy2(source/f,dest/f)
audit_inputs.P=dest;audit_inputs.B=source/'reference/B';audit_inputs.I=audit_inputs.B/'inputs';audit_inputs.OTHER=source/'reference/OTHER'
import analyze_bridge
analyze_bridge.main()
checks=[]
for name in ['sample_estimates.csv','donor_estimates.csv','paired_definition_change.csv','summary_estimates.csv','fixed_four_test_family.csv']:
 a=pd.read_csv(source/'results'/name);b=pd.read_csv(dest/'results'/name);assert a.shape==b.shape and a.columns.equals(b.columns)
 maxerr=0.
 for c in a:
  if pd.api.types.is_numeric_dtype(a[c]):
   assert np.allclose(a[c],b[c],atol=1e-12,rtol=0,equal_nan=True);maxerr=max(maxerr,float(np.nanmax(abs(a[c]-b[c]))) if a[c].notna().any() else 0)
  else:assert a[c].fillna('').equals(b[c].fillna(''))
 checks.append(dict(file=name,max_numeric_error=maxerr))
(dest/'PORTABLE_REPLAY_QA.json').write_text(json.dumps(dict(status='PASS',checks=checks),indent=2));print('PORTABLE REPLAY PASS')
