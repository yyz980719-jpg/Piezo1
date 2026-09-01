# PIEZO1 expression across human osteoarthritis contexts

Version **2.0.0** is a reproducible research compendium for the manuscript
“Human osteoarthritis reveals context-dependent PIEZO1 expression across
chondrocyte-state programmes.” It links public inputs, authoritative analysis
scripts, frozen reference outputs, figure source data, published figure images,
software metadata, and validation checks.

## What can be reproduced

- GSE51588 paired, disease-main-effect, interaction and covariate-sensitivity analyses.
- GSE104782 preprocessing, donor-aware continuous-axis, categorical and pseudobulk analyses.
- GSE57218 transferred-signature and programme-score analyses.
- GSE152805 donor-direction, within-state and Kitagawa decomposition analyses.
- Reference-internal fully nested leave-one-donor-out validation and the final recall-based sensitivity analysis.
- Whole-blood PIEZO1 cis-eQTL MR, power, influence, scale and colocalization audits.
- Main and supplementary figure rendering from audited source-data tables.
- Exploratory immune/regulatory/genetic scripts and their frozen outputs, kept separate from the authoritative pipeline.

The frozen values used in the manuscript are under `reference/`. A fresh run
writes to `results/`, `figures/`, and `metadata/` and does not overwrite the
reference layer.

## Repository structure

- `scripts/`: portable authoritative scripts, downloader, orchestrator and validators.
- `scripts/exploratory/`: historical hypothesis-generating analyses; not part of the main inferential chain.
- `reference/exploratory_figures/`: frozen historical exploratory plots.
- `data/`: input manifest, downloader instructions and one bundled derived probe map.
- `reference/results/`: frozen result tables used by the manuscript.
- `reference/source_data/`: audited figure and table source data.
- `reference/figures/published/`: exact PNG images embedded in the final manuscript and supplement.
- `reference/metadata/`: frozen QA and validation artifacts.
- `manuscript/`: final manuscript, supplement and STROBE-MR checklist used for traceability.
- `environment/`: exact observed software versions and dependency installers.
- `docs/`: provenance map, data availability text, and DOI upload instructions.

## Quick start

Requirements: R 4.6.1, Python 3.12 or later, and the R packages listed in
`environment/R-packages.tsv`.

```bash
Rscript environment/install_R_dependencies.R
python scripts/00_download_geo.py --include-large
python scripts/00_preflight.py --stage transcriptomic
python scripts/run_pipeline.py --stage transcriptomic --rscript Rscript
```

Genetic inputs are intentionally not redistributed. Obtain eQTLGen and FinnGen
files under their original access terms, place the exact regional/subset files
at the paths in `data/input_manifest.tsv`, then run:

```bash
python scripts/00_preflight.py --stage genetic
python scripts/run_pipeline.py --stage genetic --rscript Rscript
python scripts/run_pipeline.py --stage validation --rscript Rscript
python scripts/run_pipeline.py --stage figures --rscript Rscript
python scripts/99_validate_outputs.py
python scripts/98_repository_audit.py
```

Use `python scripts/run_pipeline.py --stage all --dry-run` to inspect the full
execution order without running it. Set `PIEZO1_PROJECT_ROOT`,
`PIEZO1_DATA_DIR`, `PIEZO1_RESULTS_DIR`, `PIEZO1_FIGURES_DIR`, or
`PIEZO1_METADATA_DIR` only when using non-default locations.

## Reproducibility boundary

Public GEO files can be downloaded automatically. FinnGen R12 access currently
uses the source resource's access route, and eQTLGen files remain governed by
the source resource. This archive therefore stores checksums and exact expected
filenames, not repackaged third-party raw data. The GTEx comparison script uses
the official API and caches the response locally.

Figure scripts use audited source-data CSVs to preserve exact submission layout.
The authoritative statistical tables can be regenerated from public inputs and
are checked numerically against the frozen reference layer.

The validation compares tables semantically: row/column order and gene-list
order are ignored, while numerical values are tolerance checked. Monte Carlo
weighted-median uncertainty columns are seed-controlled for future runs but
excluded from identity checking against the historical unseeded reference.
The availability-only covariate inventory is reported separately because it
changes when optional GEO series-matrix files are present.

## Citation and DOI

Citation metadata are provided in `CITATION.cff`, `.zenodo.json`, and
`codemeta.json`. No DOI is written into this version because none has yet been
assigned. After deposition, add the assigned DOI to `CITATION.cff` and the
repository landing page without retagging the already archived release.

## Licence

Repository code is licensed under MIT. Manuscript text and third-party source
data are not relicensed by this repository; see `THIRD_PARTY_NOTICES.md`.
