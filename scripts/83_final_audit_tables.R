suppressPackageStartupMessages(library(data.table))

source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
PRIOR <- PIEZO1_METADATA
OUT <- PIEZO1_METADATA

# Recompute the high-recall-state sensitivity with both Kitagawa references
# exposed, while preserving the frozen full-reference transferred labels.
qc <- fread(file.path(OUT, "04_Nested_LODO_QC.csv"))
recall <- qc[metric == "recall", .(truth = state, estimate)]
states <- sort(recall$truth)
low_recall <- recall[estimate < 0.60, truth]
retained <- setdiff(states, low_recall)
assignment <- fread(file.path(PRIOR, "01_GSE152805_StateAssignmentConfidence.csv"))
sens <- assignment[state %in% retained]

decomp_rows <- list()
direction_rows <- list()
for (dd in sort(unique(sens$donor))) {
  la <- sens[donor == dd & compartment == "lateral_nonload"]
  me <- sens[donor == dd & compartment == "medial_load"]
  pL <- as.numeric(table(factor(la$state, levels = retained))) / nrow(la)
  pM <- as.numeric(table(factor(me$state, levels = retained))) / nrow(me)
  muL <- vapply(retained, function(s) mean(la[state == s, PIEZO1]), numeric(1))
  muM <- vapply(retained, function(s) mean(me[state == s, PIEZO1]), numeric(1))
  delta <- sum(pL * muL) - sum(pM * muM)
  within_L <- sum(pL * (muL - muM)); composition_L <- sum(muM * (pL - pM))
  within_M <- sum(pM * (muL - muM)); composition_M <- sum(muL * (pL - pM))
  within_S <- (within_L + within_M) / 2; composition_S <- (composition_L + composition_M) / 2
  decomp_rows[[as.character(dd)]] <- data.table(
    donor = dd, excluded_states = paste(low_recall, collapse = ";"),
    retained_states = paste(retained, collapse = ";"),
    n_lateral = nrow(la), n_medial = nrow(me), total_delta = delta,
    within_lateral_reference = within_L, composition_lateral_reference = composition_L,
    within_medial_reference = within_M, composition_medial_reference = composition_M,
    within_symmetric = within_S, composition_symmetric = composition_S,
    within_share_symmetric = within_S / delta,
    composition_share_symmetric = composition_S / delta,
    within_share_lateral_reference = within_L / delta,
    within_share_medial_reference = within_M / delta,
    decomposition_error = abs(within_S + composition_S - delta),
    medial_higher_within = within_S < 0,
    both_components_contribute = sign(within_S) == sign(composition_S) && sign(within_S) == sign(delta))
  direction_rows[[as.character(dd)]] <- data.table(
    donor = dd, state = retained,
    n_lateral = as.numeric(table(factor(la$state, levels = retained))),
    n_medial = as.numeric(table(factor(me$state, levels = retained))),
    mean_lateral = muL, mean_medial = muM,
    delta_lateral_minus_medial = muL - muM)
}
decomp <- rbindlist(decomp_rows)
directions <- rbindlist(direction_rows)
stopifnot(max(decomp$decomposition_error) < 1e-10)
fwrite(decomp, file.path(OUT, "05_HighRecall_State_Sensitivity.csv"))
fwrite(directions, file.path(OUT, "05_HighRecall_State_Directions.csv"))
fwrite(recall[estimate < 0.60], file.path(OUT, "05_Excluded_LowRecall_States.csv"))

# Extend Data S1 with the complete deposited sample flow for GSE55235 and
# GSE82107, taken from their local GEO series-matrix metadata.
read_geo_field <- function(path, field) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  hit <- lines[startsWith(lines, field)][1]
  vals <- strsplit(hit, "\t", fixed = TRUE)[[1]][-1]
  gsub('^"|"$', "", vals)
}

old <- fread(reference_path("metadata", "Data_S1_GSM_Donor_Region_Mapping.csv"))
old <- old[!dataset %in% c("GSE55235", "GSE82107")]

f55235 <- file.path(ROOT, "data/GSE55235_series_matrix.txt.gz")
gsm55235 <- read_geo_field(f55235, "!Sample_geo_accession")
title55235 <- read_geo_field(f55235, "!Sample_title")
group55235 <- fifelse(grepl("healthy", title55235, ignore.case = TRUE), "healthy",
                      fifelse(grepl("osteoarthritic", title55235, ignore.case = TRUE), "OA", "RA"))
donor55235 <- sub(".*\\(([^)]+)\\).*", "\\1", title55235)
map55235 <- data.table(
  dataset = "GSE55235", GSM_or_sample_id = gsm55235,
  deposited_title = title55235, donor = donor55235, region = group55235,
  analysis_status = fifelse(group55235 == "RA", "Excluded", "Included"),
  reason = fifelse(group55235 == "RA",
                   "Rheumatoid-arthritis comparator outside OA-versus-healthy analysis",
                   "OA-versus-healthy synovial disease contrast"))

f82107 <- file.path(ROOT, "data/GSE82107_series_matrix.txt.gz")
gsm82107 <- read_geo_field(f82107, "!Sample_geo_accession")
title82107 <- read_geo_field(f82107, "!Sample_title")
group82107 <- fifelse(grepl("^HC", title82107), "healthy", "OA")
map82107 <- data.table(
  dataset = "GSE82107", GSM_or_sample_id = gsm82107,
  deposited_title = title82107, donor = title82107, region = group82107,
  analysis_status = "Included", reason = "OA-versus-healthy synovial disease contrast")

full_map <- rbindlist(list(old, map55235, map82107), use.names = TRUE, fill = TRUE)
setorder(full_map, dataset, GSM_or_sample_id)
fwrite(full_map, file.path(OUT, "Data_S1_GSM_Donor_Region_Mapping.csv"))

flow <- full_map[, .(deposited = .N, excluded = sum(analysis_status == "Excluded"),
                     included = sum(analysis_status == "Included")), by = dataset]
fwrite(flow, file.path(OUT, "03_GEO_SampleFlow_Audit.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "00_R_SessionInfo.txt"))
