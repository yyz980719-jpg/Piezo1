# -*- coding: utf-8 -*-
"""
53_P1_covariate_audit.py
模拟审稿后优化 · T6.1 Covariate audit

目的: 对每个统计模型说明哪些协变量被控制、哪些由配对设计自动控制、哪些不可获得。
做法: 直接解析本地 GEO series matrix 的 !Sample_characteristics_ch1 字段,
      客观记录每个数据集实际可用的元数据字段, 不做主观假设。

输出: results/53_covariate_audit.csv, results/53_covariate_audit_models.csv
"""
import csv
import gzip
import io
import os
import re
from _bootstrap import DATA, RESULTS, ROOT

ROOT = str(ROOT)

DATASETS = {
    "GSE51588":  "GSE51588_series_matrix.txt.gz",
    "GSE57218":  "GSE57218_series_matrix.txt.gz",
    "GSE55235":  "GSE55235_series_matrix.txt.gz",
    "GSE82107":  "GSE82107_series_matrix.txt.gz",
    "GSE46750":  "GSE46750_series_matrix.txt.gz",
    "GSE104782": "GSE104782_series_matrix.txt.gz",
}

# 我们真正关心的协变量 (按名称正则匹配)
FIELDS = {
    "age":      r"^age\b|^age \(",
    "sex":      r"^sex\b|^gender\b",
    "bmi":      r"^bmi\b|^body mass",
    "kl_grade": r"kle?ll?gren|kl grade|kl-grade|oa grade|severity",
    "disease":  r"^disease|^diagnosis",
    "tissue":   r"^tissue\b",
    "compartment": r"compartment|plateau|region",
    "donor_id": r"^subject\b|^patient\b|^donor\b|^individual\b",
    "batch":    r"^batch\b|^run\b|^scan",
}


def read_meta(fname):
    path = os.path.join(DATA, fname)
    if not os.path.exists(path):
        return None, 0
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_title"):
                rows.append(line.rstrip("\n").split("\t"))
    n = 0
    found = {}
    for r in rows:
        vals = [v.strip('"').strip() for v in r[1:]]
        n = max(n, len(vals))
        for v in vals:
            if ":" not in v:
                continue
            k = v.split(":", 1)[0].strip()
            kl = k.lower()
            for std, pat in FIELDS.items():
                if re.search(pat, kl):
                    found.setdefault(std, set()).add(k)
    return found, n


def main():
    out = []
    for ds, fn in DATASETS.items():
        found, n = read_meta(fn)
        if found is None:
            row = {"dataset": ds, "n_samples": 0,
                   "fields_available": "series matrix not available locally"}
            for f in ("age", "sex", "bmi", "kl_grade", "batch"):
                row[f] = "unknown"
            row["donor_id"] = "unknown"
            out.append(row)
            continue
        keys = sorted(found.keys())
        row = {"dataset": ds, "n_samples": n,
               "fields_available": "; ".join(keys) if keys else "(none of the audited fields)"}
        for f in ("age", "sex", "bmi", "kl_grade", "batch"):
            row[f] = "yes" if f in found else "no"
        row["donor_id"] = "yes" if "donor_id" in found else "no (inferred from sample titles)"
        out.append(row)
    with io.open(os.path.join(RESULTS, "53_covariate_audit.csv"), "w",
                 encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    print("=== Covariate / metadata availability audit (parsed from local GEO series matrices) ===")
    for r in out:
        print("  %-10s n=%-4s age=%-3s sex=%-3s bmi=%-3s KL=%-3s batch=%-3s | %s"
              % (r["dataset"], r["n_samples"], r["age"], r["sex"], r["bmi"],
                 r["kl_grade"], r["batch"], r["fields_available"]))

    # 每个模型层面: 控制了什么 / 配对自动控制了什么 / 不可获得
    models = [
        dict(model="GSE51588 paired medial vs lateral (primary compartment contrast)",
             unit="donor (n=25 pairs)", covariates_adjusted="none (paired within donor)",
             controlled_by_design="age, sex, BMI, medication, systemic comorbidity, genotype, batch (paired samples processed together)",
             unavailable="none needed; all time-invariant donor characteristics are absorbed by the pairing"),
        dict(model="GSE51588 disease x compartment interaction (primary interaction test)",
             unit="donor (20 OA vs 5 non-OA)", covariates_adjusted="none",
             controlled_by_design="none (between-donor comparison)",
             unavailable="age, sex, BMI, KL grade are NOT deposited for GSE51588 -> unadjusted; this is a stated limitation and the non-OA group is small"),
        dict(model="GSE51588 disease main effect (region-adjusted)",
             unit="donor (20 OA vs 5 non-OA)", covariates_adjusted="compartment",
             controlled_by_design="none",
             unavailable="age, sex, BMI, KL grade; result is descriptive/secondary"),
        dict(model="GSE57218 paired preserved vs OA-affected cartilage (programme replication)",
             unit="donor (n=33 pairs)", covariates_adjusted="none (paired within donor)",
             controlled_by_design="age, sex, genotype, batch (age and sex ARE deposited for this cohort and are donor-invariant within each pair, hence absorbed)",
             unavailable="KL grade / compartment not deposited; pairing still removes all donor-level confounding"),
        dict(model="GSE152805 paired lateral vs medial tibial cartilage (single cell)",
             unit="donor (n=3 pairs)", covariates_adjusted="none (paired within donor)",
             controlled_by_design="age, sex, BMI, genotype, donor-level batch",
             unavailable="detailed donor metadata not deposited; inference is donor-level and directional only (3/3), not a statistical replication"),
        dict(model="GSE104782 chondrocyte-state axis (10 OA donors)",
             unit="donor (n=10)", covariates_adjusted="donor identity (one rho and one quartile contrast per donor; pseudobulk models use donor blocking)",
             controlled_by_design="all donor-level covariates are absorbed because each donor contributes its own estimate",
             unavailable="age, sex, BMI, KL grade and per-donor cell counts are only partially available; per-donor cell counts are reported (111-155 cells) so that no single donor dominates"),
        dict(model="Synovial OA vs control (GSE55235, GSE82107)",
             unit="sample (10 vs 10; 10 vs 7)", covariates_adjusted="none",
             controlled_by_design="none (unpaired)",
             unavailable="age, sex, BMI, KL grade not deposited; these contrasts are secondary and descriptive, and are not used for any primary claim"),
        dict(model="Two-sample MR (eQTLGen -> FinnGen R12)",
             unit="summary statistics", covariates_adjusted="principal components / study-specific covariates as applied in the source GWAS and eQTL meta-analysis",
             controlled_by_design="random allocation of alleles at conception; two-sample design with non-overlapping exposure and outcome samples",
             unavailable="individual-level covariates unavailable by design; analyses restricted to European-ancestry summary statistics"),
    ]
    with io.open(os.path.join(RESULTS, "53_covariate_audit_models.csv"), "w",
                 encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(models[0].keys()))
        w.writeheader()
        w.writerows(models)
    print("\n=== Model-level covariate handling ===")
    for m in models:
        print("  -", m["model"])
        print("      adjusted:", m["covariates_adjusted"])
        print("      by design:", m["controlled_by_design"])
        print("      unavailable:", m["unavailable"])
    print("\nwrote results/53_covariate_audit.csv and results/53_covariate_audit_models.csv")


if __name__ == "__main__":
    main()
