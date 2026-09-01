import fs from "node:fs/promises";
import crypto from "node:crypto";
import { Workbook } from "@oai/artifact-tool";

const root = "E:/codex/2026-08-31/yu";
const pkg = `${root}/outputs/PIEZO1_OA_V10_V8_FinalExecution_Audit_20260831`;
const previewDir = `${root}/work/v10_v8_final_sequence/sheet_previews`;
await fs.mkdir(pkg, { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const sha256 = async (path) => crypto.createHash("sha256").update(await fs.readFile(path)).digest("hex");
const csvEscape = (v) => {
  const s = v == null ? "" : String(v);
  return /[",\r\n]/.test(s) ? `"${s.replaceAll('"', '""')}"` : s;
};

async function writeCsvWithQa(filename, sheetName, rows, widths = []) {
  const wb = Workbook.create();
  const sh = wb.worksheets.add(sheetName);
  const cols = Math.max(...rows.map((r) => r.length));
  const range = sh.getRangeByIndexes(0, 0, rows.length, cols);
  range.values = rows.map((r) => [...r, ...Array(cols - r.length).fill("")]);
  sh.showGridLines = false;
  sh.freezePanes.freezeRows(1);
  sh.getRangeByIndexes(0, 0, 1, cols).format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
    borders: { preset: "outside", style: "thin", color: "#9EADBA" },
  };
  sh.getRangeByIndexes(1, 0, Math.max(1, rows.length - 1), cols).format = {
    wrapText: true,
    verticalAlignment: "top",
    borders: { preset: "inside", style: "thin", color: "#D9E2F3" },
  };
  widths.forEach((w, i) => { sh.getRangeByIndexes(0, i, rows.length, 1).format.columnWidth = w; });
  sh.getRangeByIndexes(0, 0, rows.length, cols).format.autofitRows();

  const check = await wb.inspect({ kind: "table", sheetId: sheetName, range: `A1:${String.fromCharCode(64 + cols)}${rows.length}`, include: "values", tableMaxRows: Math.min(30, rows.length), tableMaxCols: cols, maxChars: 12000 });
  if (!check.ndjson.includes(String(rows[0][0]))) throw new Error(`Inspection failed for ${filename}`);
  const preview = await wb.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(`${previewDir}/${filename.replace(/\.csv$/i, "")}.png`, new Uint8Array(await preview.arrayBuffer()));
  const values = range.values;
  const csv = values.map((r) => r.map(csvEscape).join(",")).join("\r\n") + "\r\n";
  await fs.writeFile(`${pkg}/${filename}`, "\ufeff" + csv, "utf8");
}

const main = `${root}/outputs/PIEZO1_OA_Main_Manuscript_V10_Final_QA.docx`;
const supp = `${root}/outputs/PIEZO1_OA_Supplementary_Materials_V8_Final_QA.docx`;
const rows = [
  ["qa_id", "location_or_claim", "baseline_value", "authoritative_source", "source_field_or_rule", "status"],
  ["BASE-01", "Main manuscript pages", 17, "Word-rendered baseline", "page count", "LOCKED"],
  ["BASE-02", "Supplement pages", 26, "Word-rendered baseline", "page count", "LOCKED"],
  ["BASE-03", "Main figures", 5, main, "embedded figure count", "LOCKED"],
  ["BASE-04", "Supplementary figures", 15, supp, "Figure S1-S15", "LOCKED"],
  ["NUM-01", "GSE51588 disease main effect", "effect=-0.1087938612; P=0.4994569813; FDR=0.6163728728", "D:/R/projects/piezo1_oa/results/07_DEG_disease_OAvsNormal.csv", "PIEZO1 row", "SOURCE_LOCKED"],
  ["NUM-02", "GSE51588 disease-by-compartment interaction", "P=0.147", "D:/R/projects/piezo1_oa/results/07_DEG_interaction.csv", "PIEZO1 row", "SOURCE_LOCKED"],
  ["NUM-03", "Covariate-adjusted interaction", "P=0.81", "D:/R/projects/piezo1_oa/results/54_covariate_sensitivity.csv", "interaction adjusted row", "SOURCE_LOCKED"],
  ["NUM-04", "Donor-within Q1-Q4", "mean=0.123167; median=0.141223; paired-t P=0.0133911; Wilcoxon V=49, P=0.0323129", `${root}/outputs/PIEZO1_OA_V10_Final_QA_package/03_AxisDonor_Q1Q4.csv`, "scripts/71_v10_controlled_validations.R; asymptotic Wilcoxon with continuity correction", "SOURCE_LOCKED"],
  ["NUM-05", "Fully nested state-transfer LODO", "BA=0.594681; macro-F1=0.553806; null95=0.171527; permutation P=0.000999", `${pkg}/04_Nested_LODO_QC.csv`, "top-100 markers reselected inside every fold; 1,000 complete-pipeline donor-preserving permutations", "PASS"],
  ["NUM-06", "High-recall-state symmetric decomposition", "total deltas=-0.0469,-0.1382,-0.0694; within=+0.0036,-0.1013,-0.0321", `${pkg}/05_HighRecall_State_Sensitivity.csv`, "retained FC, HomC, ProC and preHTC; no reassignment; proportions renormalized", "PASS_WITH_DOWNGRADE"],
  ["NUM-07", "Overall-OA IVW weight concentration", "84.7%", "D:/R/projects/piezo1_oa/results/50_mr_ivw_weights.csv", "rs56158123", "SOURCE_LOCKED"],
  ["NUM-08", "Colocalization PP.H4", "outcome-specific", "D:/R/projects/piezo1_oa/results/41_P0_coloc_PPH3_report.csv", "PP.H4 column", "SOURCE_LOCKED"],
  ["NUM-09", "MR minimum detectable effect", "outcome-specific", "D:/R/projects/piezo1_oa/results/50_mr_power.csv; D:/R/projects/piezo1_oa/results/50_mr_power_scale_sensitivity.csv", "prespecified variance model", "SOURCE_LOCKED"],
];
await writeCsvWithQa("00_MasterQA.csv", "MasterQA", rows, [12, 30, 46, 64, 45, 18]);

const baseline = [
  "PIEZO1 OA V10/V8 BASELINE LOCK",
  "Date: 2026-08-31 (Asia/Shanghai)",
  "Main: PIEZO1_OA_Main_Manuscript_V10_Final_QA.docx",
  `Main SHA256: ${await sha256(main)}`,
  "Supplement: PIEZO1_OA_Supplementary_Materials_V8_Final_QA.docx",
  `Supplement SHA256: ${await sha256(supp)}`,
  "Rendered baseline: main 17 pages; supplement 26 pages; 5 main figures; 15 supplementary figures.",
  "Read-only copies are stored under baseline_readonly/. The V10/V8 originals are not overwritten.",
];
await fs.writeFile(`${pkg}/00_Baseline.txt`, baseline.join("\r\n") + "\r\n", "utf8");
console.log(`wrote ${pkg}/00_MasterQA.csv`);
console.log(`wrote ${pkg}/00_Baseline.txt`);
