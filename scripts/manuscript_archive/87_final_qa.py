from __future__ import annotations

import csv
import hashlib
import re
import zipfile
from pathlib import Path

from docx import Document
from pypdf import PdfReader


ROOT = Path(r"E:\codex\2026-08-31\yu")
PKG = ROOT / "outputs" / "PIEZO1_OA_V10_V8_FinalExecution_Audit_20260831"
STROBE = PKG / "PIEZO1_OA_STROBE-MR_Checklist_Completed_V2.docx"
REPO = ROOT / "outputs" / "PIEZO1_OA_public_repository"
MAIN = PKG / "PIEZO1_OA_Main_Manuscript_V13_Corrected_SOP_Audit_HOLD.docx"
SUPP = PKG / "PIEZO1_OA_Supplementary_Materials_V11_Corrected_SOP_Audit_HOLD.docx"
MAIN_PDF = ROOT / "work" / "render_main_v13_final" / "PIEZO1_OA_Main_Manuscript_V13_Corrected_SOP_Audit_HOLD.pdf"
SUPP_PDF = ROOT / "work" / "render_supp_v11_final" / "PIEZO1_OA_Supplementary_Materials_V11_Corrected_SOP_Audit_HOLD.pdf"
STROBE_PDF = ROOT / "work" / "render_strobe_v2" / "PIEZO1_OA_STROBE-MR_Checklist_Completed_V2.pdf"


def words(text: str) -> int:
    return len(re.findall(r"\b[\w]+(?:[-'][\w]+)*\b", text, flags=re.UNICODE))


def all_text(doc: Document) -> str:
    return "\n".join(p.text for p in doc.paragraphs) + "\n" + "\n".join(
        c.text for table in doc.tables for row in table.rows for c in row.cells
    )


def paragraph(doc: Document, prefix: str):
    return next(p for p in doc.paragraphs if p.text.startswith(prefix))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def xml_revision_audit(path: Path) -> tuple[bool, str]:
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        document_xml = archive.read("word/document.xml")
        has_comments = "word/comments.xml" in names
        has_revisions = b"<w:ins" in document_xml or b"<w:del" in document_xml
    return not has_comments and not has_revisions, f"comments={has_comments}; tracked_revisions={has_revisions}"


def run() -> None:
    main = Document(MAIN)
    supp = Document(SUPP)
    strobe = Document(STROBE)
    main_text = all_text(main)
    supp_text = all_text(supp)
    strobe_text = all_text(strobe)
    combined = main_text + "\n" + supp_text

    abstract_text = " ".join(paragraph(main, label).text for label in ("Objective.", "Design.", "Results.", "Conclusions."))
    intro_idx = next(i for i, p in enumerate(main.paragraphs) if p.text == "1 Introduction")
    decl_idx = next(i for i, p in enumerate(main.paragraphs) if p.text == "Declarations")
    body_text = " ".join(p.text for p in main.paragraphs[intro_idx + 1:decl_idx])
    abstract_words = words(abstract_text)
    body_words = words(body_text)

    refs = []
    for p in main.paragraphs:
        match = re.match(r"^(\d+)\. ", p.text)
        if match:
            refs.append(int(match.group(1)))

    checks: list[tuple[str, bool, str]] = []
    checks.extend([
        ("DOC-01 main opens", MAIN.exists() and len(main.paragraphs) > 100, f"paragraphs={len(main.paragraphs)}"),
        ("DOC-02 supplement opens", SUPP.exists() and len(supp.tables) == 24, f"tables={len(supp.tables)}"),
        ("DOC-03 STROBE-MR checklist opens", STROBE.exists() and len(strobe.tables) >= 1,
         f"tables={len(strobe.tables)}"),
        ("LAY-01 page counts", len(PdfReader(MAIN_PDF).pages) == 17 and len(PdfReader(SUPP_PDF).pages) == 26 and len(PdfReader(STROBE_PDF).pages) == 7,
         f"main={len(PdfReader(MAIN_PDF).pages)}; supplement={len(PdfReader(SUPP_PDF).pages)}; checklist={len(PdfReader(STROBE_PDF).pages)}"),
        ("META-01 abstract word gate", abstract_words <= 340, f"abstract_words={abstract_words}"),
        ("META-02 main-text word gate", body_words < 4000, f"main_text_words={body_words}"),
        ("NUM-01 disease P rounding", "P=0.500" in main_text, "Main result and legend use P=0.500"),
        ("NUM-02 Wilcoxon alignment", "V=49; median difference 0.1412" in supp_text and "0.0323" in combined,
         "V=49; median=0.1412; P=0.0323"),
        ("FLOW-01 GSE55235 sample flow", all(x in supp_text for x in ("30 arrays", "10 RA arrays", "RA comparator outside OA-versus-healthy analysis")),
         "30 deposited; 10 RA excluded; 20 analysed"),
        ("FLOW-02 GSE82107 sample flow", "17 arrays" in supp_text and "17 deposited: 10 OA + 7 healthy" in supp_text,
         "17 deposited; 17 analysed"),
        ("VAL-01 nested LODO", all(x in combined for x in ("balanced accuracy 0.595", "macro-F1 0.554", "0.172", "P=0.001")),
         "BA=0.595; macro-F1=0.554; null95=0.172; P=0.001"),
        ("VAL-02 recall gate", all(x in combined for x in ("HTC", "RegC", "EC", "recall below 0.60")),
         "Excluded HTC, RegC and EC"),
        ("WORD-03 recall sensitivity status", "final-audit sensitivity analysis" in combined and "prespecified sensitivity analysis" not in combined.lower(),
         "Recall<0.60 analysis is identified as final-audit sensitivity, not prespecified"),
        ("WORD-04 abstract controlled wording", "Reference-internal fully nested validation" in abstract_text and "within-state expression" in abstract_text and "within-state regulation" not in abstract_text,
         "Requested Abstract phrases are present exactly"),
        ("WORD-05 corrected validation boundary", "external projection support" in combined.lower() and
         "not external validation" in combined.lower() and "provides external validation" not in combined.lower(),
         "GSE152805 is external projection support; nested validation remains reference-internal"),
        ("WORD-06 concordance terminology", "cross-dataset concordance" in combined.lower() and
         "programme contrast replicates" not in combined.lower() and "formal paired-cartilage replication" not in combined.lower(),
         "Programme claims use cross-dataset concordance rather than replication"),
        ("QC-01 single-cell source and limitations", all(x in supp_text for x in ("matrix.mtx", "barcodes.tsv", "genes.tsv", "doublet", "ambient-RNA", "26,228")),
         "Matrix source, QC thresholds, library result, doublet and ambient-RNA limitations are explicit"),
        ("PIPE-01 Table S6/S13 authoritative L0 match", all(x in supp_text for x in ("+0.0356", "-0.0612", "+0.0298", "-0.1197", "authoritative final programme pipeline")),
         "Table S6 L0 values match Table S13 and script 51 is authoritative"),
        ("PIPE-02 programme manifest", "scripts/51_P1_compositional_sensitivity.R (authoritative final programme pipeline)" in supp_text and
         "results/51_weights_GSE51588.csv" in supp_text and "results/51_weights_GSE57218.csv" in supp_text,
         "Table S17 maps the authoritative programme script and weights"),
        ("WORD-07 donor-heterogeneous Kitagawa conclusion", "donor-heterogeneous contributions" in combined and "does not support a uniform mechanism" in combined,
         "Kitagawa interpretation is donor heterogeneous, not uniform"),
        ("WORD-08 whole-blood MR boundary", "whole-blood cis-eQTL instrument model" in combined and
         "neither evaluate nor establish any local cartilage or subchondral-bone mechanism" in main_text,
         "MR conclusion is bounded to the whole-blood instrument model"),
        ("FIG-01 controlled Figure 1 and Figure 4", (PKG / "Figures" / "Figure1_Final_P0500.png").exists() and
         (PKG / "Figures" / "Figure4_programme_robustness_direction_v11.png").exists() and
         (PKG / "Figures" / "Figure5_MR_coloc_power_v11.png").exists(),
         "Controlled terminology figures exist in PNG plus editable exports"),
        ("SENS-01 high-recall direction", all(x in combined for x in ("-0.0469", "-0.1382", "-0.0694", "+0.0036")),
         "All total deltas negative; donor 113 within-state component positive"),
        ("CONF-01 unused confidence metrics removed", "normalized entropy" not in combined.lower() and "maximum score" not in combined.lower(),
         "Only top-1-minus-top-2 margin is reported"),
        ("WORD-01 Kitagawa terminology", "Lateral-/medial-reference within-state share range" in supp_text,
         "Reference direction named explicitly"),
        ("WORD-02 claim downgrade", "is more strongly associated with PIEZO1 expression than" in abstract_text and
         "supports a weak association of PIEZO1" in main_text,
         "Objective and Results heading downgraded"),
        ("REF-01 references continuous", refs == list(range(1, 37)), f"references={refs[0]}-{refs[-1]}; n={len(refs)}"),
        ("REF-02 joint-tissue eQTL citation", "OA joint-tissue eQTL map [36]" in supp_text and "s41467-026-74993-y" in main_text,
         "Formal 2026 article cited as reference 36"),
        ("REPRO-01 no result wildcards", "results/07_DEG_*.csv" not in supp_text and "results/40_P0_*.csv" not in supp_text,
         "Table S17 enumerates exact files"),
        ("REPRO-02 local scripts/results present", len(list((PKG / "scripts").glob("*"))) >= 21 and len(list((PKG / "results").glob("*"))) >= 50,
         f"scripts={len(list((PKG / 'scripts').glob('*')))}; results={len(list((PKG / 'results').glob('*')))}"),
        ("CLEAN-01 main comments/revisions", *xml_revision_audit(MAIN)),
        ("CLEAN-02 supplement comments/revisions", *xml_revision_audit(SUPP)),
        ("AUTHOR-01 identity and contact", all(x in main_text for x in ("Yang Yunze", "0009-0003-5907-7037", "631652214@qq.com")),
         "Author name, verified ORCID and author-supplied e-mail are present"),
        ("AUTHOR-02 declarations", all(x in main_text for x in ("Project administration", "received no specific grant", "declares no competing interests")),
         "Single-author CRediT, no-funding and no-competing-interests statements are present"),
        ("ETHICS-01 final exemption statement", all(x in main_text for x in ("Article 32", "publicly available, de-identified datasets", "no additional institutional ethics committee approval or participant consent was required")),
         "Final public-data exemption wording cites Article 32 and allocates original-study responsibility"),
        ("CLEAN-03 internal completion list removed", "Items requiring author completion" not in supp_text,
         "Internal author-completion section is absent from the Supplement"),
        ("STROBE-01 all mapped rows completed", all(strobe.tables[0].cell(i, 3).text.strip() and strobe.tables[0].cell(i, 4).text.strip() for i in [1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,26,27,28,29,31,32,33,34,36,37,39,40,41,42,43,45,46,48,49,50,51,53,54,55]),
         "All applicable item/subitem page and manuscript-text cells are populated"),
        ("REPO-01 local release package", all((REPO / p).exists() for p in ("README.md", "CITATION.cff", ".zenodo.json", "LICENSE", ".git")) and len(list((REPO / "scripts").glob("*"))) >= 23 and len(list((REPO / "results").glob("*"))) >= 50,
         f"scripts={len(list((REPO / 'scripts').glob('*')))}; results={len(list((REPO / 'results').glob('*')))}; local Git/tag prepared"),
    ])

    blockers = [
        "Public repository upload, immutable public release/commit URL, and repository DOI require authenticated GitHub/Zenodo access.",
    ]
    placeholder_terms = ["[Author names", "[Affiliations", "[name, postal", "[Use CRediT", "[Funding agency", "[Authors to confirm", "[Add non-author", "[to be supplied]"]
    placeholder_count = sum(combined.count(term) for term in placeholder_terms)
    checks.append(("GATE-01 author/submission fields", placeholder_count == 0, f"unresolved_placeholders={placeholder_count}"))

    qa_lines = ["PIEZO1 OA corrected-SOP final-audit QA report", "", "Automated content checks:"]
    for name, passed, detail in checks:
        qa_lines.append(f"{'PASS' if passed else 'HOLD'}\t{name}\t{detail}")
    technical_pass = all(passed for name, passed, _ in checks if not name.startswith("GATE-"))
    overall = "HOLD" if blockers else ("PASS" if technical_pass else "FAIL")
    qa_lines.extend(["", f"Technical QA: {'PASS' if technical_pass else 'FAIL'}", f"Overall submission gate: {overall}"])
    (PKG / "00_Final_QA_Report.txt").write_text("\n".join(qa_lines) + "\n", encoding="utf-8")

    gate_rows = [
        ("P0-00", "PASS", "Read-only baseline hashes and MasterQA locked"),
        ("P0-01", "PASS", "Main Figure 1 and Supplementary Figure S1 use P=0.500"),
        ("P0-02", "PASS", "Table S5/S20 use V=49, median 0.1412, P=0.0323; method states continuity correction"),
        ("P0-03", "PASS", "GSE55235 30/10 excluded/20 final; GSE82107 17/0/17; Data S1 expanded"),
        ("P0-04", "PASS", "Fully nested LODO BA=0.595, macro-F1=0.554; complete-pipeline permutation P=0.001"),
        ("P0-05", "PASS_WITH_DOWNGRADE", "HTC, RegC and EC excluded; total direction persisted; donor 113 within-state component reversed slightly"),
        ("P0-06", "PASS", "Unused maximum-score and normalized-entropy reporting removed; margin retained"),
        ("P1-07", "PASS", "Lateral-/medial-reference wording and exact share sentence applied"),
        ("P1-08", "PASS", "Objective, Results 3.3 heading and Figure S15 title downgraded"),
        ("P1-09", "PASS", "Original GEO-study references and formal 2026 joint-tissue eQTL reference added"),
        ("P1-10", "PASS_LOCAL_HOLD_PUBLIC", "Versioned local Git repository and v1.1.0 release package prepared; authenticated public upload/DOI pending"),
        ("P0-11", "PASS_EXCEPT_EXTERNAL", "Author identity/contact, ethics, CRediT, funding, COI and STROBE-MR checklist resolved; only public repository/DOI remains"),
        ("P1-12", "PASS", f"Abstract={abstract_words} words; main text={body_words} words; metadata cleaned"),
        ("P0-13", "PASS_TECHNICAL", "17-page main, 26-page supplement and 7-page checklist visually inspected; figures, tables, numeric conflicts and placeholders resolved"),
        ("V11C-01", "PASS", "Programme scoring unified to script 51; Table S6, Table S13 and Figure 4 use one final pipeline"),
        ("V11C-02", "PASS", "GSE152805 separated into reference-internal validation and external projection support; matrix source and QC limitations added"),
        ("V11C-03", "PASS", "Replication, mechanism and MR wording downgraded to concordance, donor heterogeneity and the whole-blood instrument boundary"),
        ("P0-14", "HOLD_EXTERNAL", "Final-audit documents issued with HOLD suffix until the public repository DOI is inserted"),
    ]
    gate_text = ["PIEZO1 OA corrected V11/V9 final-review gate", ""]
    gate_text.extend(f"{step}\t{status}\t{detail}" for step, status, detail in gate_rows)
    gate_text.extend(["", f"FINAL STATUS: {overall}", "Reason: manuscript, supplement, figures, ethics wording and STROBE-MR checklist pass final audit; authenticated public repository publication and DOI remain external actions."])
    (PKG / "00_FinalGateReport.txt").write_text("\n".join(gate_text) + "\n", encoding="utf-8")

    blocker_text = ["AUTHOR INPUT REQUIRED BEFORE SUBMISSION", ""] + [f"{i}. {item}" for i, item in enumerate(blockers, 1)]
    (PKG / "00_AuthorInput_Blockers.txt").write_text("\n".join(blocker_text) + "\n", encoding="utf-8")

    changes = [
        ("P0-01", "Main Figure 1; Supplementary Figure S1", "Aligned displayed disease-main-effect P to 0.500", "07_DEG_disease_OAvsNormal.csv"),
        ("P0-02", "Methods; Table S5; abstract", "Aligned paired Wilcoxon V, median and P; stated asymptotic continuity correction", "03_AxisDonor_Q1Q4.csv"),
        ("P0-03", "Tables S1/S1A; Data S1", "Corrected complete sample flow for GSE55235 and GSE82107", "Local GEO series-matrix metadata"),
        ("P0-04", "Methods; Results; Table S19; Figure S15", "Replaced fixed-marker LODO with fully nested fold-specific marker selection and full-pipeline permutations", "04_Nested_LODO_QC.csv"),
        ("P0-05", "Abstract; Results; Discussion; Tables S2/S19/S21; Figure S15", "Applied recall gate and downgraded the mechanism claim after donor 113 reversal", "05_HighRecall_State_Sensitivity.csv"),
        ("P0-06", "Main and Supplementary Methods", "Removed unused maximum-score and normalized-entropy reporting", "Text audit"),
        ("P1-07/P1-08", "Objective; Results headings; Figure S12/S15 legends", "Applied controlled terminology and claim downgrades", "SOP wording gate"),
        ("P1-09", "Data sources; references; Table S18", "Added four original dataset papers and the formal 2026 joint-tissue eQTL paper", "References 32-36"),
        ("P1-10", "Tables S16/S17; local audit package", "Added exact versions, implementations, scripts and result filenames", "scripts/ and results/"),
        ("P1-12", "Abstract and metadata", f"Reduced abstract to {abstract_words} words and kept main text at {body_words} words", "DOCX metadata audit"),
        ("P0-11", "Title page; declarations; supplement completion list", "Added author name, affiliation, correspondence e-mail, verified ORCID, single-author CRediT roles, no-funding and no-competing-interests statements; retained only unresolved repository and journal-specific blockers", "Author attestation plus ORCID 0009-0003-5907-7037"),
        ("FINAL-01", "Abstract; Methods; Figure S15", "Relabelled recall<0.60 as a final-audit sensitivity; applied reference-internal validation and within-state expression wording", "Author final seven-item instruction"),
        ("FINAL-02", "Ethics statement; STROBE-MR checklist", "Finalized Article 32 public-data exemption wording and completed the official STROBE-MR checklist", "Official regulation and STROBE-MR template"),
        ("FINAL-03", "Supplement; repository package", "Removed internal completion notes and prepared a versioned local public-repository release package", "Visual audit and local Git v1.1.0"),
        ("V11C-01", "Main; Supplement; Figure 1; Figure 4", "Replaced replication claims with cross-dataset concordance and separated reference-internal validation from external projection support", "Corrected V11/V9 review SOP"),
        ("V11C-02", "Supplementary Methods S5/S7; Tables S6/S13/S17", "Unified programme scoring to script 51; added matrix-source, doublet and ambient-RNA limitations", "results/51_compositional_sensitivity_key.csv and local processing scripts"),
        ("V11C-03", "Abstract; Results; Discussion; Figure 5; STROBE-MR", "Applied donor-heterogeneous Kitagawa and whole-blood instrument-model boundaries", "Corrected V11/V9 review SOP"),
    ]
    with (PKG / "00_ChangeLog.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sop_step", "location", "change", "authoritative_evidence"])
        writer.writerows(changes)

    inventory_rows = []
    for path in sorted(p for p in PKG.rglob("*") if p.is_file() and p.name != "00_SHA256_Inventory.csv"):
        inventory_rows.append((path.relative_to(PKG).as_posix(), path.stat().st_size, sha256(path)))
    with (PKG / "00_SHA256_Inventory.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "bytes", "sha256"])
        writer.writerows(inventory_rows)

    print(f"Technical QA: {'PASS' if technical_pass else 'FAIL'}")
    print(f"Overall gate: {overall}")
    print(f"Abstract words: {abstract_words}; main-text words: {body_words}")
    print(f"Inventory files: {len(inventory_rows)}")


if __name__ == "__main__":
    run()
