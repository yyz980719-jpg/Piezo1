#!/usr/bin/env python3
# 23_finemap_abf.py
# Fine-mapping of the PIEZO1 locus and the 6 instrumented (consensus) TF loci
# using an ABF-based 95% credible set under a single-causal-variant assumption
# (Wakefield 2009 approximate Bayes factor). Pure stdlib; no LD reference required.
#
# Input : FinnGen R12 M13_ARTHROSIS summary stats (gz)
# Output: results/23_finemap_summary.csv
#         results/23_finemap_<GENE>.tsv   (per-SNP PP, in_credible_set)
#
# NOTE: This is a conservative fine-mapping that does NOT model LD. An LD-aware
#       method (SuSiE / POLYFUN+SuSiE / CAVIAR) with a 1000G EUR reference panel
#       would refine the credible sets; it is listed as a future direction.

import gzip, csv, math, os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OA = str(ROOT / "data" / "mr" / "finngen_R12_M13_ARTHROSIS.gz")
OUTDIR = str(ROOT / "results")
os.makedirs(OUTDIR, exist_ok=True)

PRIOR_W = 0.2          # prior variance on log(OR) for a causal variant
WINDOW = 500_000       # +/- bp around lead eQTL SNP

# gene -> (chrom, center_pos); center = lead cis-eQTL SNP position
LOCUS = {
    "PIEZO1":  ("16", 88_746_077),
    "TCF7L1":  ("2",  85_264_624),
    "GATAD2A": ("19", 19_529_715),
    "BCL6":    ("3", 187_983_973),
    "ZNF853":  ("7",  6_499_492),
    "PKNOX2":  ("11", 125_090_281),
    "ATF6":    ("1", 161_942_933),
}

# gene -> list of eQTL instrument SNP rsIDs (from MR pipeline) to test membership
EQTL_SNP = {
    "PIEZO1":  ["rs137932643", "rs149194020", "rs56158123"],  # PIEZO1 cis-eQTL IVs (11_MR_snp_level)
    "TCF7L1":  ["rs883650"],
    "GATAD2A": ["rs4808967"],
    "BCL6":    ["rs115132457"],
    "ZNF853":  ["rs147888555"],
    "PKNOX2":  ["rs7116245"],
    "ATF6":    ["rs12048680"],
}

# build per-chrom window list
windows = {}  # chrom -> list of (low, high, gene)
for g, (c, center) in LOCUS.items():
    windows.setdefault(c, []).append((center - WINDOW, center + WINDOW, g))

# buffers
buff = {g: [] for g in LOCUS}

print("Streaming OA GWAS (%.1f M rows)..." % 21.3, flush=True)
n = 0
with gzip.open(OA, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
    # expected: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta ...
    for line in f:
        n += 1
        if n % 5_000_000 == 0:
            print("  ...%d M rows" % (n // 1_000_000), flush=True)
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        chrom = parts[0]
        wl = windows.get(chrom)
        if not wl:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        for (low, high, g) in wl:
            if low <= pos <= high:
                try:
                    pval = float(parts[6])
                    beta = float(parts[8])
                    se = float(parts[9])
                except ValueError:
                    continue
                if se <= 0 or not math.isfinite(beta) or not math.isfinite(se):
                    continue
                buff[g].append((pos, parts[2], parts[3], parts[4], pval, beta, se))
                break  # one gene per row

print("Streaming done. Rows per locus:")
for g in LOCUS:
    print("  %-8s %d" % (g, len(buff[g])))


def abf_credible(snps, W=PRIOR_W):
    """snps: list of (pos,ref,alt,rsids,pval,beta,se). Returns (PP list, summary)."""
    abfs = []
    for (pos, ref, alt, rsids, pval, beta, se) in snps:
        V = se * se
        z = beta / se
        a = math.sqrt(V / (V + W)) * math.exp((z * z / 2.0) * (W / (V + W)))
        abfs.append(a)
    s = sum(abfs)
    if s <= 0:
        return [0.0] * len(snps)
    return [a / s for a in abfs]


rows_summary = []
detail_writers = {}
for g in LOCUS:
    snps = buff[g]
    if not snps:
        rows_summary.append([g, LOCUS[g][0], LOCUS[g][1], 0, "", "", 0, 0, 0, "NO_DATA", ""])
        continue
    pp = abf_credible(snps)
    # attach PP
    tagged = [(s + (p,)) for s, p in zip(snps, pp)]
    tagged.sort(key=lambda x: x[-1], reverse=True)
    cum = 0.0
    cs = []
    for t in tagged:
        cum += t[-1]
        cs.append(t)
        if cum >= 0.95:
            break
    cs_set = set(id(t) for t in cs)
    # lead by pval
    best = min(snps, key=lambda x: x[4])
    max_nlp = max(-math.log10(max(x[4], 1e-300)) for x in snps)
    lead_p = best[4]
    # eQTL instrument SNP: report its p-value and rank within the locus window
    eqtl_hits = []
    snps_sorted_p = sorted(snps, key=lambda x: x[4])  # ascending pval
    rank_map = {id(s): i + 1 for i, s in enumerate(snps_sorted_p)}
    for rs in EQTL_SNP.get(g, []):
        for s in snps_sorted_p:
            if rs in (s[3] or ""):
                eqtl_hits.append((rs, s[4], round(-math.log10(max(s[4], 1e-300)), 2), rank_map[id(s)]))
                break
    cs_fraction = len(cs) / len(snps)
    top_cs = ";".join("%s(%.3f)" % (t[3].split(";")[0], t[-1]) for t in cs[:8])
    # signal-strength tier (genome-wide sig threshold: -log10(p)=7.30)
    if max_nlp >= 7.3:
        resolv = "GWAS_significant"
    elif max_nlp >= 5.0:
        resolv = "suggestive"
    else:
        resolv = "weak"
    rows_summary.append([
        g, LOCUS[g][0], LOCUS[g][1], len(snps),
        "%s:%d-%s" % (best[2], best[0], best[3]) if best[3] else best[0],
        "%.2e" % lead_p, "%.1f" % max_nlp, len(cs), "%.3f" % cs_fraction,
        resolv,
        ("eqtl_rank:" + ";".join("%s p=%.2e nlp=%.1f rank=%d/%d" % (r, p, nlp, rk, len(snps)) for r, p, nlp, rk in eqtl_hits)) if eqtl_hits else "eqtl_not_in_window",
        top_cs,
    ])
    # write detail TSV
    with open(os.path.join(OUTDIR, "23_finemap_%s.tsv" % g), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "PP", "in_credible_set"])
        for t in tagged:
            w.writerow([t[0], t[1], t[2], t[3], t[4], t[5], t[6], "%.6f" % t[7], 1 if id(t) in cs_set else 0])

with open(os.path.join(OUTDIR, "23_finemap_summary.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "chrom", "center", "n_snps", "lead_snp", "lead_p", "max_neglog10p",
                "credible_set_size", "cs_fraction", "resolution", "eqtl_membership", "top_cs_snps"])
    w.writerows(rows_summary)

print("\n=== FINE-MAPPING SUMMARY (prior W=%.2f, window +/-%d kb) ===" % (PRIOR_W, WINDOW // 1000))
for r in rows_summary:
    print("  %-8s chr%s center=%d  n=%d  lead_p=%s  max_nlp=%s  CS_size=%s (%.3f)  %s  %s"
          % (r[0], r[1], r[2], r[3], r[5], r[6], r[7], float(r[8]), r[9], r[10]))
print("\nWrote results/23_finemap_summary.csv + per-locus TSVs.")
