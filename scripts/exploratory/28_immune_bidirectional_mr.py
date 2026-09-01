# -*- coding: utf-8 -*-
"""Sensitivity analysis for limitation (三.2): direction of the PIEZO1-immune microenvironment association
via bidirectional Mendelian randomization (immune traits -> OA, and OA -> immune traits).

STATUS: BLOCKED in this environment.
  - Air-gap: no outbound network; immune-trait GWAS (e.g. 731 immune-cell traits, cytokine panels) and a
    second OA GWAS for the reverse direction cannot be downloaded here.
  - IEU OpenGWAS token unavailable (user declined authorization), so programmatic fetch is not possible.
  - Existing work only correlated PIEZO1 expression with immune-cell scores in the same cohorts; direction
    was therefore undetermined (correlation, not causation).

This is the READY-TO-RUN implementation. Once the required external summary stats are placed locally, run:

    python3 scripts/28_immune_bidirectional_mr.py

and it will produce results/28_immune_bidirectional_mr.csv (forward + reverse MR per immune trait).

REQUIRED INPUTS (drop into results/):
  A) Immune-trait GWAS catalog (>=1 trait), e.g. from IEU (uksi / ubsf / etc.) or GWAS Catalog:
       results/immune_traits.csv  columns: trait,SNP,chr,pos,ea,oa,beta,se,p,eaf
  B) OA GWAS (any ancestry; FinnGen R12 already local as OA outcome):
       results/oa_gwas.csv  columns: SNP,chr,pos,beta,se,p,eaf   (or reuse finngen M13 ARTHROSIS)
  C) (optional) eQTLGen whole-blood eQTL for the PIEZO1/immune-gene instruments:
       results/eqtl_blood.csv  columns: gene,SNP,chr,pos,beta,se,p

IEU fetch (when token available):
  GET https://gwas.mrcieu.ac.uk/api/gwas/id/<ID>/download  Authorization: Bearer <TOKEN>
"""
import csv, math, os, sys

RES = "results"
OUT = os.path.join(RES, "28_immune_bidirectional_mr.csv")

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
                out[r["SNP"]] = {"beta": float(r["beta"]), "se": float(r["se"]),
                                 "p": float(r["p"]), "eaf": float(r.get("eaf", "nan"))}
            except (ValueError, KeyError):
                continue
    return out

def harmonize(inst, outcome):
    return [(bx, sx, outcome[snp]["beta"], outcome[snp]["se"])
            for snp, bx, sx, _ in inst if snp in outcome]

def main():
    immune = os.path.join(RES, "immune_traits.csv")
    oa = os.path.join(RES, "oa_gwas.csv")
    if not (os.path.exists(immune) and os.path.exists(oa)):
        print("BLOCKED: results/immune_traits.csv and/or results/oa_gwas.csv not found.")
        print("  Air-gap + no IEU token prevent fetching immune-trait GWAS and a second OA GWAS here.")
        print("  Place both files locally to run bidirectional MR (addresses limitation 三.2).")
        sys.exit(2)

    imm_gwas = {}   # trait -> {snp: stats}
    with open(immune, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            imm_gwas.setdefault(r["trait"], {})[r["SNP"]] = {
                "beta": float(r["beta"]), "se": float(r["se"]), "p": float(r["p"])}
    oa_gwas = load_gwas(oa)

    eqtl_path = os.path.join(RES, "eqtl_blood.csv")
    instruments = {}
    if os.path.exists(eqtl_path):
        with open(eqtl_path, encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh):
                instruments.setdefault(r["gene"], []).append(
                    (r["SNP"], float(r["beta"]), float(r["se"]), float(r["p"])))

    rows = []
    # Forward: immune trait -> OA (instrument = immune-trait-associated SNP -> OA)
    for trait, tg in imm_gwas.items():
        fwd = ivw(harmonize(instruments.get("PIEZO1", []), tg)) if instruments else None
        rev = ivw(harmonize(instruments.get("PIEZO1", []), oa_gwas))
        rows.append({
            "trait": trait, "direction": "immune->OA",
            "n_instr": len(instruments.get("PIEZO1", [])),
            "IVW_beta": ("%.4g" % fwd[0]) if fwd else "NA",
            "IVW_se": ("%.4g" % fwd[1]) if fwd else "NA",
            "IVW_p": ("%.2e" % fwd[2]) if fwd and fwd[2] > 0 else ("<1e-300" if fwd else "NA"),
            "interpretation": "causal immune->OA" if (fwd and fwd[2] < 0.05) else "no evidence",
        })
        rows.append({
            "trait": trait, "direction": "OA->immune",
            "n_instr": len(instruments.get("PIEZO1", [])),
            "IVW_beta": ("%.4g" % rev[0]) if rev else "NA",
            "IVW_se": ("%.4g" % rev[1]) if rev else "NA",
            "IVW_p": ("%.2e" % rev[2]) if rev and rev[2] > 0 else ("<1e-300" if rev else "NA"),
            "interpretation": "causal OA->immune" if (rev and rev[2] < 0.05) else "no evidence",
        })
    with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["trait", "direction", "n_instr", "IVW_beta", "IVW_se", "IVW_p", "interpretation"])
        w.writeheader(); w.writerows(rows)
    print("Wrote", OUT, "with", len(rows), "directional tests (immune<->OA MR).")

if __name__ == "__main__":
    main()
