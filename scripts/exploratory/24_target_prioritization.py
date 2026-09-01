#!/usr/bin/env python3
# 24_target_prioritization.py
# Drug-target prioritization / repurposing matrix for PIEZO1 and candidate TFs.
#
# LEFT half  = genetic evidence  -> pulled from THIS study's result CSVs (data-driven)
# RIGHT half = druggability       -> CURATED FROM LITERATURE KNOWLEDGE (no network in this
#                                    environment; Open Targets / DrugBank / ChEMBL APIs are
#                                    unreachable). MUST be re-verified when connectivity returns.
#
# Re-verify with (when network available):
#   Open Targets: https://platform-api.opentargets.org/api/v4/target/<ensemblId>
#   DrugBank, ChEMBL public downloads.

import csv, os

OUT = "results/24_target_prioritization.csv"
os.makedirs("results", exist_ok=True)

# gene, study_role, genetic_evidence, druggability_class, approved_drugs,
# clinical_development, repurposing_priority, evidence_source
ROWS = [
    ["PIEZO1",
     "Central mechanosensitive candidate (gene of focus)",
     "TWAS/SMR p_fdr=3.55e-4 (significant); MR of blood expression null (IVW); coloc PP.H4~7.2e-4 (null) -> tissue-specific mechanism",
     "Mechanosensitive PIEZO1 cation channel (ion channel)",
     "None approved",
     "Discovery-stage modulators only (Yoda1 agonist, Dooku1 antagonist; pharma PIEZO1 programs)",
     "HIGH (novel OA target; no approved drug yet)",
     "Genetic: study CSVs (20_twas_smr, 11_* , 17_*). Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["GATAD2A",
     "Consensus TF (MR + TWAS)",
     "MR beta=-0.0184 (ART, protective direction); TWAS p_fdr=1.59e-4 (significant); coloc S15",
     "Chromatin reader / NuRD-complex co-repressor",
     "None approved",
     "None reported",
     "LOW-MODERATE (epigenetic, historically hard)",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["TCF7L1",
     "Consensus TF (MR)",
     "MR beta=+0.0358 (ART, risk direction, p<1e-34 multi-SNP); TWAS nominal",
     "Wnt-pathway transcription factor (TCF/LEF family)",
     "None targeting TCF7L1 directly",
     "Wnt pathway modulated upstream (porcupine inhibitors) - not TCF7L1 directly",
     "LOW (indirect only)",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["BCL6",
     "Consensus TF (MR, weak)",
     "MR beta=+0.0447 (ART, n=2, weak); TWAS not significant in this study",
     "Transcriptional repressor; validated oncoprotein (DLBCL)",
     "None BCL6-specific approved",
     "Multiple BCL6 inhibitors in preclinical/early clinical (fx1, RI-BPI, small molecules)",
     "MODERATE-HIGH (validated target, tool compounds exist; possible immune-mechanism repurposing)",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["ZNF853",
     "Consensus TF (MR + TWAS)",
     "MR beta=+0.0255 (ART); TWAS p_fdr=3.55e-4 (significant)",
     "Zinc-finger transcription factor (poorly characterized)",
     "None",
     "None",
     "LOW",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["PKNOX2",
     "Consensus TF (MR)",
     "MR beta=+0.00765 (ART, n=10)",
     "TALE homeodomain TF (PBX/knotted-like)",
     "None",
     "None",
     "LOW",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["ATF6",
     "Consensus TF (MR)",
     "MR beta=-0.00374 (ART, n=210)",
     "ER-stress / unfolded-protein-response (UPR) transcription factor",
     "None",
     "Discovery-stage UPR modulators (neurodegeneration focus)",
     "LOW-MODERATE",
     "Genetic: study CSVs. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["RELA",
     "Bulk consensus regulator (no direct MR in this study)",
     "Bulk TF-activity + motif enrichment (Fig S9/S12) as OA-upregulated program",
     "NF-kB p65 subunit; pathway extensively targeted",
     "Many drugs modulate NF-kB pathway indirectly",
     "Numerous IKK/NF-kB inhibitors in oncology & inflammation pipelines",
     "MODERATE (pathway druggable; direct RELA inhibitors scarce)",
     "Genetic: study figures. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["RXRB",
     "Bulk consensus regulator (nuclear receptor)",
     "Bulk TF-activity + motif enrichment (Fig S9/S12)",
     "Retinoid X receptor beta (nuclear receptor)",
     "Bexarotene (RXR agonist) approved (cutaneous T-cell lymphoma)",
     "Multiple RXR modulators in development",
     "MODERATE (validated; approved-drug precedent)",
     "Genetic: study figures. Druggability: LITERATURE KNOWLEDGE - VERIFY"],

    ["NR1H2",
     "Bulk consensus regulator (nuclear receptor)",
     "Bulk TF-activity + motif enrichment (Fig S9/S12)",
     "Liver X receptor beta (nuclear receptor)",
     "None approved (LXR)",
     "LXR agonists in cardiometabolic development",
     "MODERATE",
     "Genetic: study figures. Druggability: LITERATURE KNOWLEDGE - VERIFY"],
]

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "study_role", "genetic_evidence", "druggability_class",
                "approved_drugs", "clinical_development", "repurposing_priority", "evidence_source"])
    w.writerows(ROWS)

print("Wrote", OUT, "with", len(ROWS), "genes")
print("\nRepurposing-priority ranking (this study's evidence + literature):")
for r in sorted(ROWS, key=lambda x: {"HIGH":0,"MODERATE-HIGH":1,"MODERATE":2,"LOW-MODERATE":3,"LOW":4}[x[6].split()[0]]):
    print("  %-9s %-58s %s" % (r[0], r[3][:56], r[6]))
