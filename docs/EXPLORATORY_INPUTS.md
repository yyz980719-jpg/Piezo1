# Exploratory-analysis input boundary

The scripts under `scripts/exploratory/` are retained to reproduce the
hypothesis-generating context shown in Supplementary Figure S10 and described
in the supplement. They are not invoked by the authoritative orchestrator.

Depending on the script, reruns may require additional R packages (`GSVA`,
`msigdbr`, `dorothea`, `decoupleR`, `gtexr`, `pROC`, platform annotation
packages), Python packages (`numpy`, `matplotlib`), PLINK plus a matching EUR LD
reference, full FinnGen R12 summary statistics, or external API access. These
resources remain subject to their own licences and access conditions. Frozen
outputs are available under `reference/exploratory_results/` so that the exact
reported exploratory layer remains inspectable even when an external service or
licensed input cannot be replayed.

The exploratory scripts have portable repository-root paths, but external
service version drift means that a fresh API-backed result is not guaranteed to
be byte-identical. Such differences do not alter the authoritative final
pipeline or its claims.

