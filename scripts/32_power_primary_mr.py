# -*- coding: utf-8 -*-
"""Primary PIEZO1->OA MR: statistical power and minimum detectable effect (MDE).

Addresses the key reviewer challenge to a null MR: is the null real, or simply
underpowered? Uses only the already-computed MR standard errors.

MDE: |beta|_min = (z_{1-alpha/2} + z_{1-beta}) * se,  alpha=0.05, power=80%.
Power for a given OR: Phi(|log(OR)|/se - z_{1-alpha/2}).
"""
import csv, math, os
from _bootstrap import RESULTS

RES = str(RESULTS)
SRC = os.path.join(RES, "11_MR_main_results.csv")
OUT = os.path.join(RES, "32_power_primary_mr.csv")
ALPHA = 0.05
Z_A = 1.959963985
Z80 = 0.841621234

def Phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))

def power_for_beta(se, beta_true):
    return Phi(abs(beta_true) / se - Z_A) + Phi(-abs(beta_true) / se - Z_A)

rows = []
with open(SRC, encoding="utf-8-sig", newline="") as fh:
    for r in csv.DictReader(fh):
        if r["method"] != "IVW (RE)":
            continue
        se = float(r["se"]); est = float(r["estimate"]); n = int(r["n_snp"])
        mde = (Z_A + Z80) * se
        rows.append({
            "outcome": r["outcome"], "n_instrument": n,
            "ivw_beta": round(est, 5), "ivw_se": round(se, 5),
            "mde_beta_80pct": round(mde, 4),
            "mde_OR_80pct": round(math.exp(mde), 3),
            "power_OR_1.10": round(power_for_beta(se, math.log(1.10)), 4),
            "power_OR_1.20": round(power_for_beta(se, math.log(1.20)), 4),
            "power_OR_1.50": round(power_for_beta(se, math.log(1.50)), 4),
            "power_OR_2.00": round(power_for_beta(se, math.log(2.00)), 4),
        })

with open(OUT, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print(f"{'outcome':<34}{'OR_MDE':>8}{'P(OR1.1)':>10}{'P(OR1.2)':>10}{'P(OR1.5)':>10}{'P(OR2.0)':>10}")
for r in rows:
    print(f"{r['outcome']:<34}{r['mde_OR_80pct']:>8}{r['power_OR_1.10']:>10}{r['power_OR_1.20']:>10}{r['power_OR_1.50']:>10}{r['power_OR_2.00']:>10}")
print("\nSaved ->", OUT)
