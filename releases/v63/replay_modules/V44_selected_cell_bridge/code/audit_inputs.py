from pathlib import Path
import pandas as pd, numpy as np, gzip, json, re, hashlib
P=Path(__file__).resolve().parents[1]
B=P.parent/'V42_Measurement_Response/base_V41'; I=B/'inputs'
OTHER=Path('E:/codex/2026-09-04/koa-piezo1-oa-mechanoadaptive-state-mechanopathological')
OLD=Path('E:/codex/2026-08-31/yu/outputs/PIEZO1_OA_Round1_Revision_20260908')
M5='COL2A1 SOX9 ACAN CHAD HAPLN1'.split()
M14='COL2A1 ACAN HAPLN1 COL9A1 COL9A2 COL9A3 COL11A1 COL11A2 COMP MATN3 PRG4 CILP CILP2 OGN'.split()
GENES=list(dict.fromkeys(['PIEZO1','TRPV4']+M5+M14))
def read(p,**kw): return pd.read_csv(p,sep='\t' if '.tsv' in str(p) else ',',**kw)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 for d in ['audit','inputs','results','qa','figures','manuscript']: (P/d).mkdir(parents=True,exist_ok=True)
 q=read(OTHER/'work/phase5/bridge_cell_QC_v2.csv.gz');q=q[q.ID.str.match(r'OA[1-8]-[12]$')].copy()
 tail=q.original_barcode.str.extract(r'([ACGT]{16}-\d+)$',expand=False);assert tail.notna().all()
 q['cell_id']=q.ID+'::'+tail;assert q.cell_id.is_unique and q.barcode.is_unique
 q['E2_selected']=q.phase5_pass_qc_v2.astype(str).str.lower().eq('true')
 meta=read(I/'audit/p07_prefreeze/candidate_design_16x9.tsv');rows=[];members=[];sources=[]
 for m in meta.itertuples():
  x=read(I/f'audit/p07_execution/E3_cell_input_{m.legacy_id}.tsv.gz');c=read(I/f'audit/p09_exploration/channel_counts_{m.legacy_id}.tsv.gz')
  assert x.cell_id.equals(c.cell_id) and np.array_equal(x.PIEZO1,c.PIEZO1)
  sub=q[q.ID.eq(m.legacy_id)].set_index('cell_id');assert sub.new_ID.eq(m.source_id).all() and sub.donor_ID.eq(m.legacy_id.split('-')[0]).all() and sub.condition.eq(m.region).all()
  assert set(x.cell_id)<=set(sub.index)
  sel=sub.loc[x.cell_id];f=OTHER/f'work/counts_v2/{m.legacy_id}_target149_bridge.csv.gz'
  a=read(f,index_col=0);assert a.index.is_unique and set(GENES)<=set(a.index) and set(sel.barcode)<=set(a.columns)
  counts=a.loc[GENES,sel.barcode].to_numpy().T
  assert np.array_equal(counts[:,[GENES.index(g) for g in M5]],x[M5].to_numpy())
  assert np.array_equal(counts[:,:2],c[['PIEZO1','TRPV4']].to_numpy())
  assert np.array_equal(sel.counts_v2_total.to_numpy(),x.total_counts.to_numpy())
  intersect=sel.E2_selected.to_numpy();expr=np.log1p(1e4*counts/x.total_counts.to_numpy()[:,None])
  np.savez_compressed(P/f'inputs/{m.legacy_id}.npz',counts=counts,expression=expr,genes=GENES,barcodes=x.cell_id.to_numpy(dtype=str),full_UMI=x.total_counts.to_numpy(),intersection=intersect)
  sub['V42_selected']=sub.index.isin(x.cell_id);sub['sample']=m.legacy_id;sub['region']=m.region;sub['donor']=m.donor_id
  members.append(sub.reset_index()[['cell_id','barcode','sample','donor','region','V42_selected','E2_selected']])
  rows.append(dict(cohort='discovery',sample=m.legacy_id,source_sample=m.source_id,GSM=m.GSM,donor=m.donor_id,E2_legacy_donor=m.legacy_id.split('-')[0],region=m.region,V42_cells=len(x),E2_cells=int(sub.E2_selected.sum()),intersection=int(intersect.sum()),V42_only=int((~intersect).sum()),E2_only=int((sub.E2_selected&~sub.V42_selected).sum()),all_target_counts_and_totals_equal=True))
  sources.append(dict(path=str(f),sha256=sha(f)));print('AUDIT',m.legacy_id,len(x),int(intersect.sum()),flush=True)
 pd.concat(members).to_csv(P/'audit/discovery_cell_membership.csv.gz',index=False)
 pd.DataFrame(rows).to_csv(P/'audit/discovery_sample_audit.csv',index=False)
 assert sum(r['V42_cells'] for r in rows)==101702 and sum(r['E2_cells'] for r in rows)==119207
 external=[]
 for rec in json.loads((I/'audit/p21_external_effects/input_manifest.json').read_text()):
  gsm=rec['GSM'];f=next((OLD/'audit/p17_external').glob(gsm+'*features.tsv.gz'));bc=next((OLD/'audit/p17_external').glob(gsm+'*barcodes.tsv.gz'))
  feat=read(f,header=None);ids={g:np.flatnonzero(feat[1].eq(g)) for g in GENES};assert all(len(v)==1 for v in ids.values())
  cell=read(I/f'audit/p21_external_effects/{gsm}_normalized_targets_and_masks.tsv.gz');cell=cell[cell.primary_P19.astype(str).str.lower().eq('true')]
  bars=gzip.open(bc,'rt').read().splitlines();assert len(bars)==len(set(bars)) and set(cell.barcode)<=set(bars)
  with gzip.open(rec['matrix_path'],'rt') as h:
   line=h.readline()
   while line.startswith('%'):line=h.readline()
  external.append(dict(GSM=gsm,primary_cells=len(cell),dimensions=line.strip(),matrix_path=rec['matrix_path'],expected_sha256=rec['matrix_sha256'],feature_path=str(f),barcode_path=str(bc),all_target_identities_unique=True))
  print('EXTERNAL AUDIT',external[-1],flush=True)
 assert sum(r['primary_cells'] for r in external)==35389
 (P/'audit/external_input_audit.json').write_text(json.dumps(external,indent=2))
 (P/'audit/discovery_sources.json').write_text(json.dumps(sources,indent=2))
 print('INPUT IDENTITY AUDIT PASS; no associations calculated',flush=True)
if __name__=='__main__':main()
