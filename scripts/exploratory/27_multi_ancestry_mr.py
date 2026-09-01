# -*- coding: utf-8 -*-
"""Sensitivity analysis for limitations (一.3) multi-ancestry MR and (二.1) tissue-specific eQTL re-instrumentation.

STATUS: RUNNABLE.
  - The required non-European OA GWAS (GCST90566800, Hatzikotoulas et al. Nat Genet 2025, multi-ancestry
    OA meta-analysis, hg19) was provided by the user and converted to results/oa_nonEU.csv via
    scripts/31a_lookup.py + 31b_harmonize.py. Allele harmonization (flip / strand-flip) and
    SE-from-p fallback were applied upstream so that the OA effect is aligned to the instrument
    effect allele (fg_alt) — the IVW below therefore needs no further harmonization.
  - All existing MR in this project used FinnGen R12 (Finnish / European) and eQTLGen (whole-blood European).
    This sensitivity addresses limitation (一.3) multi-ancestry MR.

Run:

    python3 scripts/27_multi_ancestry_mr.py

and it will produce results/27_multi_ancestry_mr.csv + a supplementary-table-ready block.

REQUIRED INPUTS (drop into results/):
  A) Non-European OA GWAS, e.g. East-Asian OA meta or UKB OA:
       results/oa_nonEU.csv  columns: SNP,chr,pos,ea,oa,beta,se,p,eaf,n
  B) (optional, for 二.1) Tissue-specific eQTL (cartilage / subchondral / synovial / chondrocyte):
       results/eqtl_tissue.csv  columns: gene,SNP,chr,pos,beta,se,p
     If present, the script re-selects TF/PIEZO1 instruments from tissue eQTL instead of blood eQTLGen.

The IEU fetch path is provided (commented) for when a token becomes available:
  GET https://gwas.mrcieu.ac.uk/api/gwas/id/<ID>/download  with header Authorization: Bearer <TOKEN>
"""
import csv, math, os, sys

RES = "results"
OUT = os.path.join(RES, "27_multi_ancestry_mr.csv")

def erfc_surv(z):
    return math.erfc(abs(z) / math.sqrt(2.0))

def ivw(rows):
    sw = swb = 0.0
    for bx, sx, by, sy in rows:
        if sx <= 0 or sy <= 0:
            continue
        w = (bx * bx) / (sy * sy)
        sw += w; swb += w * by
    if sw <= 0:
        return None
    b = swb / sw; se = 1.0 / math.sqrt(sw); z = b / se
    return b, se, erfc_surv(z)

def load_gwas(path):
    out = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                out[r["SNP"]] = {
                    "chr": r.get("chr"), "pos": int(r["pos"]) if r.get("pos") else None,
                    "beta": float(r["beta"]), "se": float(r["se"]), "p": float(r["p"]),
                    "eaf": float(r.get("eaf", "nan")),
                }
            except (ValueError, KeyError):
                continue
    return out

def load_instruments(path):
    """TF/PIEZO1 cis-eQTL instruments from 17_MR_TF_OA_multi_snp.csv.
    Columns: tf,outcome,tissue,SNP,fg_alt,fg_ref,alt_g,beta_x,se_x,beta_y,se_y,p_y,eaf,wald,se_wald
    Here: gene=tf, exposure effect=beta_x (eQTL, wrt fg_alt), exposure se=se_x.
    """
    inst = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            try:
                inst.setdefault(r["tf"], []).append(
                    (r["SNP"], float(r["beta_x"]), float(r["se_x"]), float(r["p_y"])))
            except (ValueError, KeyError):
                continue
    return inst

def harmonize(instr, outcome):
    """Return list of (beta_x, se_x, beta_y, se_y) for overlapping, aligned SNPs."""
    rows = []
    for snp, bx, sx, px in instr:
        o = outcome.get(snp)
        if not o:
            continue
        rows.append((bx, sx, o["beta"], o["se"]))
    return rows

def main():
    oa = os.path.join(RES, "oa_nonEU.csv")
    if not os.path.exists(oa):
        print("BLOCKED: results/oa_nonEU.csv not found.")
        print("  Provide a non-European OA GWAS (e.g. GCST90566800) converted to results/oa_nonEU.csv")
        print("  with columns SNP,chr,pos,ea,oa,beta,se,p,eaf,n (OA effect aligned to instrument effect allele).")
        print("  (Addresses limitation 一.3 multi-ancestry MR and 二.1 tissue-specific eQTL re-instrumentation.)")
        sys.exit(2)

    oa_gwas = load_gwas(oa)
    eqtl_src = os.path.join(RES, "eqtl_tissue.csv")
    instr_src = eqtl_src if os.path.exists(eqtl_src) else os.path.join(RES, "17_MR_TF_OA_multi_snp.csv")
    instruments = load_instruments(instr_src)

    rows = []
    for gene, snps in instruments.items():
        harm = harmonize([(s[0], s[1], s[2], s[3]) for s in snps], oa_gwas)
        iv = ivw(harm)
        if not iv:
            rows.append({"gene": gene, "n_instrument": 0, "IVW_beta": "NA", "IVW_se": "NA", "IVW_p": "NA"})
            continue
        b, se, p = iv
        rows.append({"gene": gene, "n_instrument": len(harm),
                     "IVW_beta": "%.4g" % b, "IVW_se": "%.4g" % se,
                     "IVW_p": ("%.2e" % p if p > 0 else "<1e-300")})
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "n_instrument", "IVW_beta", "IVW_se", "IVW_p"])
        w.writeheader(); w.writerows(rows)
    print("Wrote", OUT, "with", len(rows), "genes (non-European OA MR).")

if __name__ == "__main__":
    main()
