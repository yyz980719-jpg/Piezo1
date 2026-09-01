# -*- coding: utf-8 -*-
"""Sensitivity analyses for manuscript limitation (二.2): TF->OA MR LD-clumping / multiple-testing.

S18: LD-clump-INDEPENDENT lead single cis-eQTL IVW per TF x OA definition + Benjamini-Hochberg (BH) q.
S19: multi-SNP IVW WITHOUT LD clumping (numeric p) + BH q, contrasted with lead -> demonstrates P-inflation
     when instruments are not clumped. Plus a PIEZO1 distance-clump (chr/pos available) worked example.

All data local; no external download. Distance clumping uses 1 Mb window (PIEZO1 only, which has chr/pos).
TF multi-SNP lacks chr/pos offline, so lead-instrument is used as the clump-equivalent (one tag SNP per locus).
"""
import csv, math, os
from collections import defaultdict

RES = "results"
OUT_ZH = "（解读见正文局限 §二.2）"

def erfc_surv(z):
    return math.erfc(abs(z) / math.sqrt(2.0))

def bh_q(pvals):
    """Return list of BH q-values aligned to input order."""
    ps = [p for p in pvals if p is not None and p == p]
    order = sorted(range(len(pvals)), key=lambda i: (pvals[i] if pvals[i] is not None else 9))
    qs = [None] * len(pvals)
    m = sum(1 for p in pvals if p is not None)
    if m == 0:
        return qs
    # ranks among non-None, ascending
    ranked = [i for i in order if pvals[i] is not None]
    prev = 1.0
    for rank, idx in enumerate(reversed(ranked), 1):  # rank from largest p
        r = m - rank + 1
        q = pvals[idx] * m / r
        q = min(q, prev)
        prev = q
        qs[idx] = q
    return qs

def ivw(rows):
    """Inverse-variance weighted MR from per-SNP (beta_x, se_x, beta_y, se_y)."""
    sw = 0.0; swb = 0.0
    for bx, sx, by, sy in rows:
        if sx <= 0 or sy <= 0:
            continue
        w = (bx * bx) / (sy * sy)
        sw += w; swb += w * by
    if sw <= 0:
        return None
    b = swb / sw
    se = 1.0 / math.sqrt(sw)
    z = b / se
    p = erfc_surv(z)
    return b, se, p

# ---------- load ----------
sup = list(csv.DictReader(open(os.path.join(RES, "17_MR_TF_OA_supplementary_table.csv"), encoding="utf-8-sig")))
multi = list(csv.DictReader(open(os.path.join(RES, "17_MR_TF_OA_multi_snp.csv"), encoding="utf-8-sig")))

OUTMAP = {"ART": "All-site OA", "KNEE": "Knee OA", "COX": "Hip OA"}
# sup-table Outcome labels -> 17_multi outcome codes
OMAP = {"Overall OA": "ART", "Knee OA": "KNEE", "Hip OA": "COX"}
def ocode(s):
    return OMAP.get(s["Outcome"], s["Outcome"])

# group per-SNP by (tf, outcome)
per_tf = defaultdict(list)
for r in multi:
    try:
        bx = float(r["beta_x"]); sx = float(r["se_x"]); by = float(r["beta_y"]); sy = float(r["se_y"])
    except ValueError:
        continue
    per_tf[(r["tf"], r["outcome"])].append((bx, sx, by, sy))

# ---------- S18: lead + BH ----------
lead_ps = []
for s in sup:
    try:
        lead_ps.append(float(s["Lead_P"]))
    except ValueError:
        lead_ps.append(None)
q_lead = bh_q(lead_ps)

s18 = []
for i, s in enumerate(sup):
    s18.append({
        "TF": s["TF"],
        "OA_definition": OUTMAP.get(s["Outcome"], s["Outcome"]),
        "Lead_cis_eQTL_SNP": s["Lead_cis_eQTL_SNP"],
        "eQTL_tissue": s["Best_GTEx_tissue"],
        "Lead_IVW_beta_logOR": s["Lead_IVW_beta_logOR"],
        "Lead_95CI": s["Lead_95pct_CI"],
        "Lead_P": s["Lead_P"],
        "BH_q_24tests": ("%.4g" % q_lead[i]) if q_lead[i] is not None else "NA",
        "Significant_BH": "yes" if (q_lead[i] is not None and q_lead[i] < 0.05) else "no",
    })

# ---------- S19: multi (unclumped) + BH, vs lead ----------
multi_res = {}
for (tf, outc), rows in per_tf.items():
    iv = ivw(rows)
    multi_res[(tf, outc)] = iv  # (b, se, p) or None

multi_ps = []
for s in sup:
    key = (s["TF"], ocode(s))
    iv = multi_res.get(key)
    multi_ps.append(iv[2] if iv else None)
q_multi = bh_q(multi_ps)

s19 = []
for i, s in enumerate(sup):
    key = (s["TF"], ocode(s))
    iv = multi_res.get(key)
    if iv:
        b, se, p = iv
        z = b / se
        ci = "[%.4g, %.4g]" % (b - 1.96 * se, b + 1.96 * se)
        multi_beta = "%.4g" % b; multi_ci = ci; multi_p = "%.3g" % p
    else:
        multi_beta = "NA"; multi_ci = "NA"; multi_p = "NA (no cis-eQTL instruments)"
    s19.append({
        "TF": s["TF"],
        "OA_definition": OUTMAP.get(s["Outcome"], s["Outcome"]),
        "n_SNP_unclumped": s["Multi_n_SNP"],
        "Multi_IVW_beta_logOR": multi_beta,
        "Multi_95CI": multi_ci,
        "Multi_P": multi_p,
        "BH_q_24tests": ("%.4g" % q_multi[i]) if q_multi[i] is not None else "NA",
        "Lead_P": s["Lead_P"],
        "Verdict": ("significant but UNCLUMPED -> hypothesis-generating only"
                    if (q_multi[i] is not None and q_multi[i] < 0.05)
                    else "not significant"),
    })

# ---------- PIEZO1 distance-clump worked example ----------
pie = list(csv.DictReader(open(os.path.join(RES, "11_MR_snp_level.csv"), encoding="utf-8-sig")))
pie_rows = []
for r in pie:
    try:
        pie_rows.append({
            "chr": r["chr"], "pos": int(r["pos"]),
            "SNP": r["SNP"],
            "bx": float(r["beta_x"]), "sx": float(r["se_x"]),
            "by": float(r["beta_y"]), "sy": float(r["se_y"]),
            "py": float(r["p_y"]),
        })
    except (ValueError, KeyError):
        continue

def ivw_pie(rows):
    sw = 0.0; swb = 0.0
    for d in rows:
        w = (d["bx"] ** 2) / (d["sy"] ** 2)
        sw += w; swb += w * d["by"]
    if sw <= 0:
        return None
    b = swb / sw; se = 1.0 / math.sqrt(sw); z = b / se; p = erfc_surv(z)
    return b, se, p

full = ivw_pie(pie_rows)
# distance clump 1 Mb on same chr, keep most significant per cluster
order = sorted(pie_rows, key=lambda d: d["py"])
kept = []
used = []
for d in order:
    if any(d["chr"] == u["chr"] and abs(d["pos"] - u["pos"]) < 1_000_000 for u in used):
        continue
    kept.append(d); used.append(d)
pruned = ivw_pie(kept)

piezo1 = {
    "n_full": len(pie_rows),
    "n_pruned": len(kept),
    "full": full, "pruned": pruned,
}

# ---------- write CSVs ----------
def write_csv(name, rows):
    with open(os.path.join(RES, name), "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

write_csv("26_s18_lead_bh.csv", s18)
write_csv("26_s19_multi_bh.csv", s19)

with open(os.path.join(RES, "26_piezo1_clump.csv"), "w", encoding="utf-8-sig", newline="") as fh:
    fh.write("metric,value\n")
    fh.write("n_SNP_full,%d\n" % piezo1["n_full"])
    fh.write("n_SNP_after_1Mb_distance_clump,%d\n" % piezo1["n_pruned"])
    if full: fh.write("IVW_beta_full,%.5g\nIVW_se_full,%.5g\nIVW_P_full,%.4g\n" % (full[0], full[1], full[2]))
    if pruned: fh.write("IVW_beta_clumped,%.5g\nIVW_se_clumped,%.5g\nIVW_P_clumped,%.4g\n" % (pruned[0], pruned[1], pruned[2]))

# ---------- console summary ----------
print("=== S18 lead-instrument + BH (24 tests) ===")
nsig = sum(1 for x in s18 if x["Significant_BH"] == "yes")
print("TF x OA tests:", len(s18), "| significant after BH:", nsig)
print("min BH q:", min(float(x["BH_q_24tests"]) for x in s18 if x["BH_q_24tests"] != "NA"))
print("\n=== S19 multi (unclumped) + BH ===")
nsig_m = sum(1 for x in s19 if x["Verdict"].startswith("significant"))
print("multi tests with p:", sum(1 for x in s19 if "NA" not in x["Multi_P"]), "| 'significant but unclumped':", nsig_m)
print("\n=== PIEZO1 distance-clump example ===")
print("full n=%d beta=%.4g p=%.3g | clumped n=%d beta=%.4g p=%.3g"
      % (piezo1["n_full"], full[0], full[2], piezo1["n_pruned"], pruned[0], pruned[2]))
print("\nWrote: 26_s18_lead_bh.csv, 26_s19_multi_bh.csv, 26_piezo1_clump.csv")
