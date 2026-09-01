# PIEZO1 OA reproducible repository — Zenodo multipart upload

Version: 2.0.0  
Reserved DOI: 10.5281/zenodo.22231510  
Creator: Yang, Yunze  

These four ZIP files are a transport split of the single reproducibility archive. They are independent ZIP archives, not binary split volumes. Upload all four files to the same Zenodo record.

## Files

1. `PIEZO1_OA_v2.0.0_part1_code_and_results.zip`
   - Repository-level metadata, licensing, citation files, code, environment lockfiles, public-input manifest, frozen results and exploratory result tables.
2. `PIEZO1_OA_v2.0.0_part2_metadata_and_source_data.zip`
   - Figure metadata, source-data tables and validation metadata.
3. `PIEZO1_OA_v2.0.0_part3_figures.zip`
   - Publication-ready main and supplementary figure files.
4. `PIEZO1_OA_v2.0.0_part4_manuscript_and_exploratory_figures.zip`
   - Final manuscript materials and exploratory figure outputs.

## Reconstruction

Extract all four ZIP files into the same directory. Each archive preserves the common top-level directory `PIEZO1_OA_reproducible_repository_v2.0.0/`, so the repository is reconstructed by merging that directory across the four extractions.

The original single-file archive is retained locally as `PIEZO1_OA_reproducible_repository_v2.0.0.zip` with SHA-256:

`cd1a22eb4a8b1c59e7d2c698451f4453d9495fc419ea90cd90538b5b6d6a67d7`

Verify the four upload files against `SHA256SUMS.txt` before upload. Zenodo should show a completed upload, a file size and a calculated checksum for every file before the record is published.
