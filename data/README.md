# Input data

The repository does not redistribute third-party raw data. Required paths,
source identifiers, access routes, byte sizes and SHA-256 checksums are listed
in `input_manifest.tsv`.

Run `python scripts/00_download_geo.py --include-large` to retrieve the public
GEO inputs and unpack GSE152805. The bundled
`GPL13497_probe2symbol.tsv` is a two-column derivative of the public GPL13497
annotation and is included to make probe aggregation deterministic.

Genetic inputs must be obtained from eQTLGen and FinnGen using their official
access routes. Place them under `data/mr/` with the filenames in the manifest.
The regional FinnGen tables are GRCh38 extracts for chromosome 16 matching the
analysis scripts; extraction boundaries and checksums are recorded in the
manifest and provenance documentation.

Generated intermediates:

- `data/scrna/GSE104782_seurat.rds` is created by `scripts/10_scrna_piezo1.R`.
- `data/mr/gtex_v8_wholeblood_PIEZO1.json` is created or reused by
  `scripts/50b_P1_gtex_scale_validation.py`.

