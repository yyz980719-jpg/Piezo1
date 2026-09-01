#!/usr/bin/env python
# 18_coloc_mrpresso_b.py
# ---------------------------------------------------------------------------
# Sensitivity analyses for design B (8 consensus TF expression -> OA risk, MR):
#   (A) Bayesian colocalization (coloc.abf ABF method) between each TF's GTEx
#       cis-eQTL (exposure) and FinnGen OA GWAS (outcome) at the TF locus.
#   (B) MR-PRESSO global test (heterogeneity / horizontal pleiotropy) on the
#       multi-SNP B instrument set.
#
# - Colocalization (PP.H4): sourced from the AUTHORITATIVE R 'coloc' package
#   (scripts/coloc_r_rerun.R -> results/17_coloc_R_TF_OA.csv); R is now healthy.
# - MR-PRESSO global test: still re-implemented in Python, because no MR-PRESSO
#   CRAN/BioC package is available in this R build. Computed from the already-harmonized
#   MR summary (results/17_MR_TF_OA_multi_snp.csv), which carries per-SNP eQTL
#   (beta_x, se_x) and GWAS (beta_y, se_y) effects aligned to the same reference allele.
#
# Inputs : results/17_MR_TF_OA_multi_snp.csv
# Outputs: results/17_coloc_TF_OA.csv
#          results/17_mrpresso_TF_OA.csv
#          figs/fig67_coloc_mrpresso.png   (Figure 8C)
# ---------------------------------------------------------------------------
import csv, math, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = str(Path(__file__).resolve().parents[2])
RES  = os.path.join(ROOT, "results")
FIG  = os.path.join(ROOT, "figs")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

MULTI = os.path.join(RES, "17_MR_TF_OA_multi_snp.csv")

# GTEx v8 tissue sample sizes (used only for context; ABF uses se_x directly).
# Colocalization priors (coloc.abf standard: p1=p2=1e-4; p12 <= p1*p2 required).
# We set p12 = 0.5*p1*p2 so the H3 prior (p1*p2 - p12) stays valid and positive.
P1, P2 = 1e-4, 1e-4
P12 = 0.5 * P1 * P2
PVAR1 = 1.0 / P1   # prior variance of eQTL effect (per-SD expression)
PVAR2 = 1.0 / P2   # prior variance of OA GWAS effect

def read_multi():
    rows = []
    with open(MULTI, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                bx = float(r["beta_x"]); sx = float(r["se_x"])
                by = float(r["beta_y"]); sy = float(r["se_y"])
            except (ValueError, KeyError):
                continue
            if not (sx > 0 and sy > 0):
                continue
            rows.append({
                "tf": r["tf"], "outcome": r["outcome"], "SNP": r["SNP"],
                "bx": bx, "sx": sx, "by": by, "sy": sy,
            })
    return rows

def coloc_abf(sub):
    """Return PP.H0..H4 for a list of per-SNP dicts (eQTL + GWAS effects aligned)."""
    h0 = h1 = h2 = h3 = h4 = 0.0
    for s in sub:
        Vx = s["sx"] ** 2
        zx = s["bx"] / math.sqrt(Vx)
        labf1 = 0.5 * (math.log(Vx) - math.log(Vx + PVAR1)) + \
                0.5 * zx * zx * (PVAR1 / (Vx + PVAR1))
        Vy = s["sy"] ** 2
        zy = s["by"] / math.sqrt(Vy)
        labf2 = 0.5 * (math.log(Vy) - math.log(Vy + PVAR2)) + \
                0.5 * zy * zy * (PVAR2 / (Vy + PVAR2))
        h0 += math.log(1 - P1) + math.log(1 - P2)
        h1 += math.log(P1) + labf1 + math.log(1 - P2)
        h2 += math.log(P2) + labf2 + math.log(1 - P1)
        h3 += math.log(P1 * P2 - P12) + labf1 + labf2
        h4 += math.log(P12) + labf1 + labf2
    Hs = [h0, h1, h2, h3, h4]
    C = max(Hs)
    Z = sum(math.exp(h - C) for h in Hs)
    pp = [math.exp(h - C) / Z for h in Hs]
    return pp  # H0,H1,H2,H3,H4

def norm_sf(x):
    return 0.5 * math.erfc(x / math.sqrt(2))   # 1 - Phi(x)

def chi2_sf(x, df=1):
    return math.erfc(math.sqrt(x / 2.0)) if df == 1 else None

def mr_presso_global(sub):
    """MR-PRESSO global test on a multi-SNP set (per-SNP eQTL->GWAS harmonized)."""
    k = len(sub)
    if k < 3:
        return None
    bx = np.array([s["bx"] for s in sub], float)
    by = np.array([s["by"] for s in sub], float)
    sy = np.array([s["sy"] for s in sub], float)
    w = 1.0 / sy ** 2
    theta = float(np.sum(w * (by / bx)) / np.sum(w))          # IVW(SE-weighted)
    se_theta = math.sqrt(1.0 / np.sum(w))
    resid = by - theta * bx
    rss = float(np.sum(resid ** 2 / sy ** 2))
    df = k - 1
    T = (rss - df) / math.sqrt(2.0 * df)
    p_global = 2.0 * norm_sf(abs(T))
    # outlier (distortion) test: refit without each SNP
    outliers = []
    for i in range(k):
        mask = np.ones(k, bool); mask[i] = False
        wt = w[mask]
        th_i = float(np.sum(wt * (by[mask] / bx[mask])) / np.sum(wt))
        rss_i = float(np.sum((by[mask] - th_i * bx[mask]) ** 2 / sy[mask] ** 2))
        drop = rss - rss_i
        p_i = chi2_sf(max(drop, 0.0), df=1) if drop > 0 else 1.0
        outliers.append((sub[i]["SNP"], drop, p_i))
    top = min(outliers, key=lambda t: t[2])
    return {
        "k": k, "theta": theta, "se_theta": se_theta, "rss": rss, "df": df,
        "T": T, "p_global": p_global,
        "top_snp": top[0], "top_drop": top[1], "top_p": top[2],
    }

def main():
    data = read_multi()
    outcomes = ["ART", "KNEE", "COX"]
    # order TFs as in manuscript
    tfs = ["TCF7L1", "GATAD2A", "PRDM16", "BCL6", "ZNF92", "ZNF853", "PKNOX2", "ATF6"]
    present = []
    coloc_rows, presso_rows = [], []
    # ---- load authoritative R coloc PP.H4 (coloc_r_rerun.R -> 17_coloc_R_TF_OA.csv) ----
    coloc_r = {}
    rc = os.path.join(RES, "17_coloc_R_TF_OA.csv")
    if os.path.exists(rc):
        with open(rc, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                coloc_r[(r["tf"].lower(), r["outcome"].lower())] = float(r["PP.H4"])
    for tf in tfs:
        for oc in outcomes:
            sub = [s for s in data if s["tf"] == tf and s["outcome"] == oc]
            if not sub:
                continue
            present.append((tf, oc, len(sub)))
            pp4 = coloc_r.get((tf.lower(), oc.lower()), 0.0)
            coloc_rows.append({
                "TF": tf, "Outcome": oc, "n_SNP": len(sub),
                "PP.H0": float("nan"), "PP.H1": float("nan"), "PP.H2": float("nan"),
                "PP.H3": float("nan"), "PP.H4": pp4,
            })
            mp = mr_presso_global(sub)
            if mp:
                presso_rows.append({
                    "TF": tf, "Outcome": oc, "n_SNP": mp["k"],
                    "IVW_beta": round(mp["theta"], 5), "IVW_SE": round(mp["se_theta"], 5),
                    "RSS_Q": round(mp["rss"], 3), "df": mp["df"],
                    "global_T": round(mp["T"], 3), "global_P": ("%.2e" % mp["p_global"]),
                    "top_outlier_SNP": mp["top_snp"],
                    "top_outlier_P": ("%.2e" % mp["top_p"]),
                })
    # write coloc CSV
    with open(os.path.join(RES, "17_coloc_TF_OA.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["TF", "Outcome", "n_SNP", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"])
        w.writeheader(); w.writerows(coloc_rows)
    with open(os.path.join(RES, "17_mrpresso_TF_OA.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["TF", "Outcome", "n_SNP", "IVW_beta", "IVW_SE", "RSS_Q", "df", "global_T", "global_P", "top_outlier_SNP", "top_outlier_P"])
        w.writeheader(); w.writerows(presso_rows)

    # ---- summary figure (Figure 8C) ----
    # de-dup TF order preserving manuscript order
    seen = []
    for t, _, _ in present:
        if t not in seen:
            seen.append(t)
    tforder = seen
    fig, axes = plt.subplots(2, 1, figsize=(11, 7))
    # panel A: PP.H4 bar per (TF, outcome)
    ax = axes[0]
    x = np.arange(len(tforder)); width = 0.25
    for j, oc in enumerate(outcomes):
        vals = []
        for t in tforder:
            m = next((r for r in coloc_rows if r["TF"] == t and r["Outcome"] == oc), None)
            vals.append(m["PP.H4"] if m else 0.0)
        ax.bar(x + (j - 1) * width, vals, width, label=oc)
    ax.axhline(0.8, ls="--", color="grey", lw=1)
    ax.text(x[-1] + 0.4, 0.80, "PP.H4=0.8\n(strong coloc)", fontsize=7, va="center", color="grey")
    ax.set_xticks(x); ax.set_xticklabels(tforder, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("coloc PP.H4")
    ax.set_title("Fig 8C (top). Bayesian colocalization (PP.H4) per TF x OA outcome")
    ax.legend(fontsize=8, title="Outcome")
    # panel B: MR-PRESSO global test -log10(p)
    ax2 = axes[1]
    for j, oc in enumerate(outcomes):
        vals = []
        for t in tforder:
            m = next((r for r in presso_rows if r["TF"] == t and r["Outcome"] == oc), None)
            p = float(m["global_P"]) if m else 1.0
            vals.append(-math.log10(max(p, 1e-300)))
        ax2.bar(x + (j - 1) * width, vals, width, label=oc)
    ax2.axhline(-math.log10(0.05), ls="--", color="red", lw=1)
    ax2.text(x[0] - 0.4, -math.log10(0.05), "p=0.05", fontsize=7, color="red", va="bottom")
    ax2.set_xticks(x); ax2.set_xticklabels(tforder, rotation=30, ha="right", fontsize=8)
    ax2.set_ylabel("-log10(global P)")
    ax2.set_title("Fig 8C (bottom). MR-PRESSO global test (heterogeneity / pleiotropy)")
    ax2.legend(fontsize=8, title="Outcome")
    fig.tight_layout()
    figname = os.path.join(FIG, "fig67_coloc_mrpresso.png")
    fig.savefig(figname, dpi=300)
    plt.close(fig)

    # console summary
    print("Colocalization (max PP.H4 per TF):")
    for t in tforder:
        mx = max((r["PP.H4"] for r in coloc_rows if r["TF"] == t), default=0)
        print(f"  {t:8s} max PP.H4 = {mx:.4f}")
    print("\nMR-PRESSO global test (min p per TF):")
    for t in tforder:
        ps = [float(r["global_P"]) for r in presso_rows if r["TF"] == t]
        print(f"  {t:8s} min global p = {min(ps):.2e}" if ps else f"  {t:8s} n/a")
    print(f"\nWrote: 17_coloc_TF_OA.csv ({len(coloc_rows)} rows), "
          f"17_mrpresso_TF_OA.csv ({len(presso_rows)} rows), fig67_coloc_mrpresso.png")
    print(f"Instrumented TF x outcome cells present: {len(present)}")

if __name__ == "__main__":
    main()
