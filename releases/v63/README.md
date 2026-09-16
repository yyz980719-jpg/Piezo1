# PIEZO1–TRPV4 cartilage association: V63 reproducibility materials

This versioned research-code release supports the current selected-cell analyses and five figure source datasets. It contains code, selected-cell arrays, donor/sample mappings, archived numerical results, figure source data and replay instructions. **It does not contain the unpublished manuscript or supplementary manuscripts.** Earlier repository releases and the older Zenodo record are not replaced.

Download the release asset `PIEZO1_V63_public_reproducibility.zip`, not GitHub's automatically generated source-code archive. The latter does not contain these large selected-cell inputs. Verify `SHA256SUMS.txt` before extraction. `FILE_MANIFEST.json` lists the individual file hashes.

## Scope

The release starts from frozen selected-cell count/expression arrays and inherited labels from GSE255460 and GSE220243. It can recompute the channel correlations, donor summaries, fixed statistical families, measurement diagnostics, member/state comparisons and nested bootstrap described in the replay guide. It is **not** a FASTQ-to-results pipeline: it does not independently re-run alignment, initial cell calling, upstream QC, annotation or doublet-model fitting. Fixed upstream masks are supplied where used. The two sources retain their distinct sampling and selection systems.

The bootstrap holds the six donors, regions and labels fixed. Its percentile ranges and sign fractions are conditional cell-resampling diagnostics, not population confidence intervals or probabilities that an effect is true. Failed draws remain explicitly reported. Historical identifiers such as `external` are preserved in frozen files; they do not imply confirmatory validation. No new hypothesis families were introduced in this release.

## Start here

1. Python 3.12 and `pip install -r requirements-v63.txt`.
2. Follow `REPLAY.md` from the extracted package root.
3. Use `REPRODUCTION_MAP.csv` to find each figure's source and computational entry.
4. Read `THIRD_PARTY_NOTICES.md` and `VERIFICATION_SCOPE.md`.

Original code: MIT (see `LICENSE`). Source data, external software and prior-source annotations retain their original rights and attribution; the code license is not a blanket license over those materials. Research use only; no clinical prediction tool is provided.
