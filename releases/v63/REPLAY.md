# Replay instructions

Run from the extracted release root in a separate working copy. Tested on Windows, Python 3.12.14, with the direct dependency pins in `requirements-v63.txt`. These are observed versions, not a complete cross-platform environment lock. Arial was used for figure rendering; font substitution can change layout, not the tabulated estimates.

## Current nested stability diagnostic

```text
python analysis/run.py --output qa/replay
python verify_bootstrap.py
```

The first command requires a fresh output directory, validates the frozen plan and input hashes, then uses four worker processes. The second checks 19 replay CSVs against archived results and checks saved resampling indices with an alternative numerical implementation. `analysis/inputs` retains 16 libraries for provenance; the fixed 12 libraries from six donors are selected by the frozen receipt. No donor is resampled.

## Earlier cell-level calculations used by the current figures

```text
python replay_modules/V44_selected_cell_bridge/code/replay_bridge.py --out qa/bridge_replay
python replay_modules/V46_members/verify_present.py
python replay_modules/V35_measurement/Review_bundle/code/replay_measurement.py qa/measurement_replay
python replay_modules/V35_measurement/Review_bundle/code/replay_review_outputs.py qa/measurement_replay qa/measurement_diagnostics_replay
python replay_modules/V50_member_measurement/run.py --output replay_modules/V50_member_measurement/qa_replay
python replay_modules/V50_member_measurement/verify.py
python replay_modules/V51_member_state/run.py --output replay_modules/V51_member_state/qa_replay
python replay_modules/V51_member_state/verify.py
python replay_modules/V61_RepC_depth/analysis/run.py --output qa/depth_replay
python numerical_verification/core_verify.py --archived
```

V44 recomputes M5/M14 correlations and fixed-family summaries and compares five tables. V46 independently recalculates 748 channel/member correlations and checks 238 donor estimates; its historical `new_test_count` describes the original module, not a new test selected for this release. V35 recomputes measurement sensitivity, including the fixed-depth thinning draws, and the second command adds absolute-correlation and detection diagnostics. V50 and V51 re-run the member measurement and inherited-state analyses and verify their archived tables. V61 re-runs the fixed-six PRG4 depth diagnostic. The final command checks 154 archived-summary identities; by itself it is not a cell-level replay.

`verify.py` for V50/V51 expects the precise `qa_replay` location above. Do not overwrite a completed analysis directory; start a new extraction for a clean repeat. Historical receipts retain original software/path provenance. Any nonportable source paths in guarded staging sections are not needed by the packaged replay entry points.

## Rebuild figures from saved estimates

```text
python prepare_and_plot.py
python plot_figure5_v63.py
```

Run these in this order: the first generates the shared figure set from saved tables; the second replaces Figure 5 with the current four-panel, stability-integrated version. Figure generation is a presentation replay, not a substitute for recalculating correlations. Manuscript builders and manuscript documents are intentionally not distributed.

## Portability changes

The V44 analysis script no longer imports a local skill module or calls its unused normality helper; the reported Shapiro/IQR diagnostics still come from the original SciPy/NumPy calculations. Five saved output tables were reproduced after this change. V46 now creates its figure output directory if missing. No eligibility threshold, input value, scoring definition, test family or random seed was changed.

Floating-point and CSV byte equality may differ with other numerical-library versions. Inspect and document discrepancies instead of silently relaxing tolerances. Retain all failed-draw denominators when interpreting the stability outputs.
