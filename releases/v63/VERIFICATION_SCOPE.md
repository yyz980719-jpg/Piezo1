# Internal verification scope

Before publication, the sanitized release staging directory was used to run the commands in `REPLAY.md`.

- Nested bootstrap: all 19 CSVs byte-identical to the archived calculation; 236 estimates checked by an alternative implementation, four failed-gate cases checked, maximum numerical difference 6.67e-16.
- Core numerical summaries: 154/154 checks passed.
- M5/M14 bridge: five preserved output tables reproduced; M5 replay maximum difference below 1e-15.
- Member map: 748 correlations and 238 donor estimates verified; maximum correlation discrepancy below 1e-15.
- Full-member measurement checks: 14 replay CSVs matched; numerical discrepancies below 1e-14.
- Inherited-state analysis: nine replay tables matched; numerical discrepancies below 1e-14.
- Fixed-six PRG4 depth analysis: 432 historical checks passed; 48 sample estimates estimable.
- Measurement and thinning diagnostic scripts completed from the packaged selected-cell inputs.

These are internal computational checks, not independent human review, external biological replication, or upstream raw-read validation. Original and replay receipts provide the detailed scope. Failed bootstrap draws are preserved, not discarded from the reported denominator.
