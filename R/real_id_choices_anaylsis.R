# =============================================================================
# REAL IDENTIFICATION CHOICES – PATIENT DATA ANALYSIS
# =============================================================================
#
# Analysis of answer choices (identification responses) from clinical
# olfactory test records.  The dataset contains TDI scores and per-item
# answer choices for 352 patients.
#
# =============================================================================

# ---- External functions and parameters --------------------------------------

resolve_project_dir <- function() {
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }

  if (file.exists("globals.R")) {
    return(getwd())
  }

  if (file.exists(file.path("R", "globals.R"))) {
    return(file.path(getwd(), "R"))
  }

  stop("Could not locate globals.R relative to the script or repository root.")
}

source(file.path(resolve_project_dir(), "globals.R"))


# =============================================================================
# 0. Packages and global options
# =============================================================================

if (!requireNamespace("psych", quietly = TRUE)) {
  stop(
    "Package 'psych' is required. Install it with:\n",
    "install.packages('psych')"
  )
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  stop(
    "Package 'ComplexHeatmap' is required. Install it with:\n",
    "install.packages('ComplexHeatmap')"
  )
}

if (!requireNamespace("cABCanalysis", quietly = TRUE)) {
  stop(
    "Package 'cABCanalysis' is required. Install it with:\n",
    "install.packages('cABCanalysis')"
  )
}

options(stringsAsFactors = FALSE)

# =============================================================================
# 1. Data loading and cleaning
# =============================================================================

Patient_data <- read.csv(
  "~/.Datenplatte/Joerns Dateien/Aktuell/RiechenRandom/09Originale/Patient data from ORBIS-352 patients+anosmiaiden.csv"
)

# Drop Excel intermediate-formula artefact columns (named "X.").
Patient_data <- Patient_data[, -grep("X.", names(Patient_data))]

# Retain only the 352 data rows; trailing rows are NA artefacts from the export.
Patient_data <- Patient_data[1:352, ]

Patient_data$ID      <- 1:nrow(Patient_data)
Patient_data$Anosmia <- ifelse(Patient_data$TDI.score < 16.75, 1, 0)


# =============================================================================
# 2. Initial exploration
# =============================================================================

psych::describe(Patient_data)

par(mfrow = c(2, 1))
plot(density(Patient_data$Age,       na.rm = TRUE))
plot(density(Patient_data$TDI.score, na.rm = TRUE))
par(mfrow = c(1, 1))


# =============================================================================
# 3. Heatmap of correct-answer matrix
# =============================================================================

correct_matrix <- Patient_data[, c(
  "TDI.score",
  names(Patient_data)[grep("Pen.", names(Patient_data))]
)]
correct_matrix <- correct_matrix[, -grep("Answer", names(correct_matrix))]
correct_matrix <- na.omit(correct_matrix)

ha_rows_tdi <- ComplexHeatmap::rowAnnotation(
  "TDI" = anno_barplot(correct_matrix$TDI.score),
  show_legend = FALSE,
  width = unit(4, "cm")
)

Heatmap(
  correct_matrix[, 2:ncol(correct_matrix)],
  cluster_rows              = TRUE,
  cluster_columns           = TRUE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  col                       = ggthemes::colorblind_pal()(8),
  right_annotation          = ha_rows_tdi
)


# =============================================================================
# 4. ABC analysis of per-item correct-answer counts
# =============================================================================

sum_per_pen <- apply(correct_matrix[, -1], 2, sum)

ABC_sum_per_pen <- cABC_analysis(sum_per_pen, PlotIt = TRUE)
names(sum_per_pen[ABC_sum_per_pen$Aind])


# =============================================================================
# 5. Single-choice data
# =============================================================================

single_choices_data <- Patient_data[, c(
  "ID",
  "Anosmia",
  names(Patient_data)[grep("Answer", names(Patient_data))]
)]

single_choices_complete_data <- single_choices_data[
  !apply(single_choices_data == "", 1, any), , drop = FALSE
]

head(single_choices_complete_data)


# =============================================================================
# 6. Encode answer strings as choice positions (1–4)
# =============================================================================

# For each trial, the four response alternatives in presentation order (1–4).
trial_choices <- list(
  `1`  = c("Ananas",       "Brombeere",    "Erdbeere",    "Orange"),
  `2`  = c("Rauch",        "Klebstoff",    "Schuhleder",  "Gras"),
  `3`  = c("Honig",        "Zimt",         "Schokolade",  "Vanille"),
  `4`  = c("Pfefferminz",  "Schnittlauch", "Fichte",      "Zwiebel"),
  `5`  = c("Banane",       "Kokos",        "Walnuss",     "Kirsche"),
  `6`  = c("Pfirsich",     "Zitrone",      "Apfel",       "Grapefruit"),
  `7`  = c("Kaugummi",     "Gummib.",      "Lakritz",     "Kekse"),
  `8`  = c("Senf",         "Gummi",        "Menthol",     "Terpentin"),
  `9`  = c("Zwiebel",      "Sauerkraut",   "Möhren",      "Knoblauch"),
  `10` = c("Zigarette",    "Wein",         "Kaffee", 	  "Kerzenrauch"),
  `11` = c("Melone",  	   "Apfel",        "Orange",      "Pfirsich"),
  `12` = c("Gewürznelke",  "Pfeffer",      "Zimt",        "Senf"),
  `13` = c("Ananas",       "Pflaume",      "Pfirsich",    "Birne"),
  `14` = c("Kamille",      "Rose",         "Himbeere",    "Kirsche"),
  `15` = c("Honig",        "Rum",          "Anis",        "Fichte"),
  `16` = c("Brot",         "Schinken",     "Kase",        "Fisch")
)

# Lowercase, strip German umlauts, and remove trailing punctuation so that
# truncated forms (e.g. "Gummib.") and umlaut-less encodings ("Mohren") resolve.
normalize_str <- function(x) {
  x <- trimws(tolower(x))
  x <- gsub("ä", "a", x, fixed = TRUE)
  x <- gsub("ö", "o", x, fixed = TRUE)
  x <- gsub("ü", "u", x, fixed = TRUE)
  x <- gsub("ß", "ss", x, fixed = TRUE)
  x <- gsub("ffee", "ffe", x, fixed = TRUE)
  sub("[[:punct:]]+$", "", x)
}


# Return position 1–4 of answer in choices, or NA if no match found.
# Exact match is tried first; prefix match handles abbreviations.
match_answer_to_position <- function(answer, choices) {
  if (is.na(answer) || !nzchar(trimws(answer))) return(NA_integer_)

  a_norm <- normalize_str(answer)
  c_norm <- normalize_str(choices)

  exact <- which(a_norm == c_norm)
  if (length(exact) == 1L) return(as.integer(exact))

  prefix <- which(startsWith(a_norm, c_norm))
  if (length(prefix) == 1L) return(as.integer(prefix))

  NA_integer_
}


answer_cols <- grep("Answer", names(single_choices_complete_data), value = TRUE)
trial_nums  <- as.integer(regmatches(answer_cols, regexpr("\\d+", answer_cols)))
answer_cols <- answer_cols[order(trial_nums)]

encoded_matrix <- single_choices_complete_data[, c("ID", "Anosmia")]

for (i in seq_along(answer_cols)) {
  encoded_matrix[[paste0("Pen.", i)]] <- vapply(
    single_choices_complete_data[[answer_cols[i]]],
    match_answer_to_position,
    choices     = trial_choices[[i]],
    FUN.VALUE   = integer(1)
  )
}

na_counts <- colSums(is.na(encoded_matrix[, -(1:2), drop = FALSE]))
if (any(na_counts > 0)) {
  warning(
    "Some answer strings could not be matched to a choice position:\n",
    paste(names(na_counts[na_counts > 0]),
          na_counts[na_counts > 0],
          sep = ": ", collapse = "\n")
  )
}

head(encoded_matrix)

# Coerce all columns to integer (Anosmia from ifelse() is double by default).
encoded_matrix <- as.data.frame(lapply(encoded_matrix, as.integer))

cat(sprintf("Response patterns available from %s patients\n", dim(encoded_matrix)[1]))


# =============================================================================
# 7. Heatmap of encoded choice matrix
# =============================================================================

pen_cols <- grep("^Pen\\.", names(encoded_matrix), value = TRUE)

# =============================================================================
# 8. Save and re-read encoded matrix
# =============================================================================

write.csv(encoded_matrix, "encoded_choice_matrix.csv", row.names = FALSE)

encoded_matrix_loaded <- read.csv("encoded_choice_matrix.csv")


# =============================================================================
# 9. Choice-sequence strings
# =============================================================================

pen_cols_loaded <- grep("^Pen\\.", names(encoded_matrix_loaded), value = TRUE)
pen_cols_loaded <- pen_cols_loaded[
  order(as.integer(regmatches(pen_cols_loaded, regexpr("\\d+", pen_cols_loaded))))
]

sequence_strings <- apply(
  encoded_matrix_loaded[, pen_cols_loaded],
  1,
  function(row) paste0("c(", paste(row, collapse = ", "), ")")
)

names(sequence_strings) <- encoded_matrix_loaded$ID

head(sequence_strings)


# =============================================================================
# 10. Compact sequence pattern frequency analysis
# =============================================================================

# Strip all non-digit characters → "c(4, 3, 2, ...)" becomes "4324..."
compact_sequences <- gsub("[^0-9]", "", sequence_strings)

pattern_counts <- sort(table(compact_sequences), decreasing = TRUE)

n_patients   <- length(compact_sequences)
n_unique     <- length(pattern_counts)
n_duplicated <- sum(pattern_counts > 1)

cat(sprintf(
  "Patients: %d | Unique patterns: %d (%.1f%%) | Patterns shared by ≥2 patients: %d\n",
  n_patients, n_unique, 100 * n_unique / n_patients, n_duplicated
))

# Under purely random independent choice (4 alternatives, 16 trials):
# E[collision pairs] ≈ C(n,2) / 4^16  (birthday-problem approximation)
n_possible         <- 4^16
expected_collision_pairs <- choose(n_patients, 2) / n_possible
observed_collision_pairs <- sum(choose(as.integer(pattern_counts), 2))

cat(sprintf(
  "Expected collision pairs under random choice: %.5f\nObserved collision pairs: %d\n",
  expected_collision_pairs, observed_collision_pairs
))

# Barplot – patterns shared by at least 2 patients
repeated <- pattern_counts[pattern_counts > 1]

if (length(repeated) == 0) {
  message("All response patterns are unique — consistent with random choice behaviour.")
} else {
  df_rep <- data.frame(
    pattern = factor(names(repeated), levels = names(repeated)),
    count   = as.integer(repeated)
  )

  print(
    ggplot(df_rep, aes(x = pattern, y = count)) +
      geom_col(fill = "steelblue", width = 0.7) +
      labs(
        title    = "Repeated response patterns in the Sniffin' Sticks identification test",
        subtitle = sprintf(
          "%d of %d patients share a sequence pattern with at least one other patient",
          sum(df_rep$count), n_patients
        ),
        x = "Sequence pattern",
        y = "Number of patients"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))
  )
}


# =============================================================================
# 11. Foreseeable ("deterministic strategy") response patterns
# =============================================================================
#
# Rationale: a patient unable to smell the stimuli may abandon guessing and
# fall back on a simple rule for filling in the forced-choice answers.  Such
# rules produce sequences that are foreseeable a priori, e.g. always the same
# alternative (1111111111111111), block-wise runs (1111222233334444 or
# 4444333322221111), cyclic sweeps (1234123412341234), alternations
# (1212121212121212) or zigzag sweeps (1234321234321234).
#
# All templates below are enumerated a priori from the design of the test
# (L trials, K response alternatives) and not derived from the observed data,
# so the test is confirmatory.  For each template we report the number of
# exact hits, the hit rate, the number of patients within Hamming distance 1
# and 2 (near misses), the smallest observed distance, and the number of
# exact hits expected under uniform random choice (n * K^-L).
#
# =============================================================================

L <- length(pen_cols_loaded)   # trials per patient (16)
K <- 4L                        # forced-choice alternatives per trial

# Scoring key of the Sniffin' Sticks identification test, i.e. the position of
# the correct alternative in each trial (also used in Sections 12 and 13).
correct_key <- c(
  4, 3, 2, 1, 1, 2, 3, 4,
  4, 3, 2, 1, 1, 2, 3, 4
)


# -----------------------------------------------------------------------------
# 11.1 A priori catalogue of foreseeable templates
# -----------------------------------------------------------------------------

rotate  <- function(x, k) x[((seq_along(x) - 1L + k) %% length(x)) + 1L]
recycle <- function(base, n = L) rep_len(base, n)

templates <- data.frame(
  family  = character(0),
  pattern = character(0),
  stringsAsFactors = FALSE
)

add_template <- function(family, values) {
  if (length(values) != L) return(invisible(NULL))
  templates <<- rbind(
    templates,
    data.frame(
      family  = family,
      pattern = paste(values, collapse = ""),
      stringsAsFactors = FALSE
    )
  )
}

# (a) Constant – the same alternative in every trial
for (d in seq_len(K)) add_template("Constant", rep(d, L))

# (b) Cyclic sweep, ascending and descending, all K starting points
#     e.g. 1234123412341234 / 4321432143214321
for (p in 0:(K - 1)) {
  add_template("Cyclic ascending",  recycle(rotate(seq_len(K), p)))
  add_template("Cyclic descending", recycle(rotate(rev(seq_len(K)), p)))
}

# (c) Doubled cyclic sweep – each alternative held for two trials
#     e.g. 1122334411223344 / 4433221144332211
for (p in 0:(K - 1)) {
  add_template("Doubled cyclic ascending",  recycle(rep(rotate(seq_len(K), p), each = 2)))
  add_template("Doubled cyclic descending", recycle(rep(rotate(rev(seq_len(K)), p), each = 2)))
}

# (d) Block patterns – each alternative held for L/K consecutive trials
#     e.g. 1111222233334444 / 4444333322221111
if (L %% K == 0) {
  for (p in 0:(K - 1)) {
    add_template("Block ascending",  rep(rotate(seq_len(K), p),      each = L / K))
    add_template("Block descending", rep(rotate(rev(seq_len(K)), p), each = L / K))
  }
}

# (e) Alternation between two alternatives, all ordered pairs
#     e.g. 1212121212121212, 4141414141414141
for (a in seq_len(K)) {
  for (b in seq_len(K)) {
    if (a != b) add_template("Two-alternative alternation", recycle(c(a, b)))
  }
}

# (f) Two constant halves – one alternative for the first half of the test,
#     another for the second, e.g. 1111111122222222
if (L %% 2 == 0) {
  for (a in seq_len(K)) {
    for (b in seq_len(K)) {
      if (a != b) add_template("Two constant halves", c(rep(a, L / 2), rep(b, L / 2)))
    }
  }
}

# (g) Zigzag ("ping-pong") sweeps, e.g. 1234321234321234
add_template("Zigzag ascending",  recycle(c(seq_len(K), rev(seq_len(K))[-c(1, K)])))
add_template("Zigzag descending", recycle(c(rev(seq_len(K)), seq_len(K)[-c(1, K)])))

# (h) Palindromic ("tent"/"V") sweeps, e.g. 1234432112344321 / 4321123443211234
add_template("Palindromic sweep", recycle(c(seq_len(K), rev(seq_len(K)))))
add_template("Palindromic sweep", recycle(c(rev(seq_len(K)), seq_len(K))))

templates <- templates[!duplicated(templates$pattern), , drop = FALSE]

# The V-shaped template 4321123443211234 coincides with the scoring key of the
# test.  A hit there is therefore indistinguishable from genuinely correct
# responding and must be reported separately rather than as guessing strategy.
templates$equals_key <- templates$pattern == paste(correct_key, collapse = "")

n_templates <- nrow(templates)


# -----------------------------------------------------------------------------
# 11.2 Matching the observed sequences against the catalogue
# -----------------------------------------------------------------------------

# Only complete sequences can be matched; incomplete ones carry "NA" entries
# that are lost when the digits are extracted (Section 10).
observed_sequences <- compact_sequences[nchar(compact_sequences) == L]
n_obs <- length(observed_sequences)

if (n_obs < length(compact_sequences)) {
  message(sprintf(
    "%d of %d sequences were incomplete and excluded from the pattern analysis.",
    length(compact_sequences) - n_obs, length(compact_sequences)
  ))
}

obs_mat <- matrix(
  as.integer(unlist(strsplit(observed_sequences, "", fixed = TRUE))),
  nrow  = n_obs,
  byrow = TRUE
)

tpl_mat <- matrix(
  as.integer(unlist(strsplit(templates$pattern, "", fixed = TRUE))),
  nrow  = n_templates,
  byrow = TRUE
)

# D[i, j] = Hamming distance between patient i and template j
D <- matrix(NA_integer_, nrow = n_obs, ncol = n_templates)
for (j in seq_len(n_templates)) {
  D[, j] <- colSums(t(obs_mat) != tpl_mat[j, ])
}

hits_per_template <- colSums(D == 0)

pattern_results <- data.frame(
  Family               = templates$family,
  Pattern              = templates$pattern,
  Hits                 = as.integer(hits_per_template),
  Hit_rate_percent     = round(100 * hits_per_template / n_obs, 2),
  Within_distance_1    = as.integer(colSums(D <= 1)),
  Within_distance_2    = as.integer(colSums(D <= 2)),
  Min_distance         = as.integer(apply(D, 2, min)),
  Expected_hits_random = signif(n_obs / K^L, 3),
  Equals_scoring_key   = templates$equals_key,
  stringsAsFactors     = FALSE
)

pattern_results <- pattern_results[
  order(-pattern_results$Hits, pattern_results$Min_distance, pattern_results$Family),
  ,
  drop = FALSE
]

cat(sprintf(
  "\nForeseeable-pattern analysis: %d a priori templates tested against %d complete response sequences (L = %d trials, K = %d alternatives)\n\n",
  n_templates, n_obs, L, K
))

print(pattern_results, row.names = FALSE)

write.csv(pattern_results, "foreseeable_pattern_results.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 11.3 Summary statistics for reporting
# -----------------------------------------------------------------------------

# Probability that a uniformly random sequence lies within Hamming distance d
# of one given template.
prob_within <- function(d) sum(choose(L, 0:d) * (K - 1)^(0:d)) / K^L

patient_min_distance <- apply(D, 1, min)

n_hits_any    <- sum(patient_min_distance == 0)
n_within1_any <- sum(patient_min_distance <= 1)
n_within2_any <- sum(patient_min_distance <= 2)

p_any_template <- n_templates * prob_within(0)   # templates are distinct
binom_any <- binom.test(n_hits_any, n_obs, p = p_any_template, alternative = "greater")

cat(sprintf(
  paste0(
    "\nPatients matching any of the %d templates exactly: %d/%d (%.2f%%); ",
    "expected under uniform random choice: %.4g (p = %.3g, exact binomial, one-sided)\n",
    "Patients within Hamming distance 1 of a template: %d/%d (%.2f%%); expected <= %.3g\n",
    "Patients within Hamming distance 2 of a template: %d/%d (%.2f%%); expected <= %.3g\n",
    "Distance to the nearest template: median %d, minimum %d (chance level %.1f)\n"
  ),
  n_templates, n_hits_any, n_obs, 100 * n_hits_any / n_obs,
  n_obs * p_any_template, binom_any$p.value,
  n_within1_any, n_obs, 100 * n_within1_any / n_obs, n_obs * n_templates * prob_within(1),
  n_within2_any, n_obs, 100 * n_within2_any / n_obs, n_obs * n_templates * prob_within(2),
  median(patient_min_distance), min(patient_min_distance), L * (K - 1) / K
))

# Hits broken down by template family (excluding the scoring-key template,
# where a hit cannot be attributed to a response strategy).
family_input <- pattern_results[!pattern_results$Equals_scoring_key, , drop = FALSE]
family_input$Templates <- 1L

family_summary <- aggregate(
  cbind(Templates, Hits, Within_distance_1, Within_distance_2) ~ Family,
  data = family_input,
  FUN  = sum
)
family_summary$Hit_rate_percent <- round(100 * family_summary$Hits / n_obs, 2)

cat("\nHits by template family (scoring-key template excluded):\n\n")
print(family_summary, row.names = FALSE)

if (any(pattern_results$Equals_scoring_key)) {
  key_row <- pattern_results[pattern_results$Equals_scoring_key, , drop = FALSE]
  cat(sprintf(
    "\nNote: template %s equals the scoring key (%d hits = fully correct answer sheets, not a guessing strategy).\n",
    key_row$Pattern[1], key_row$Hits[1]
  ))
}


# -----------------------------------------------------------------------------
# 11.4 Template-free indicators of rule-based responding
# -----------------------------------------------------------------------------
#
# The catalogue above is necessarily incomplete.  The following indicators
# capture the same idea without fixing a template: local transition rules
# (stay / step up / step down), long constant runs, and use of only part of
# the response scale.  Each is compared with uniform random choice, either
# analytically or against simulated random sequences.

set.seed(42)
n_sim <- 20000L
sim_mat <- matrix(sample.int(K, n_sim * L, replace = TRUE), ncol = L)

rate_stay <- mean(obs_mat[, -L] == obs_mat[, -1])
rate_up   <- mean(((obs_mat[, -L] %% K) + 1L) == obs_mat[, -1])
rate_down <- mean(obs_mat[, -L] == ((obs_mat[, -1] %% K) + 1L))

n_transitions <- n_obs * (L - 1)

transition_summary <- data.frame(
  Transition_rule = c("stay (x[i+1] = x[i])",
                      "step up (x[i+1] = x[i] + 1 mod K)",
                      "step down (x[i+1] = x[i] - 1 mod K)"),
  Observed_rate   = round(c(rate_stay, rate_up, rate_down), 4),
  Expected_rate   = round(1 / K, 4),
  p_value         = signif(c(
    binom.test(round(rate_stay * n_transitions), n_transitions, 1 / K)$p.value,
    binom.test(round(rate_up   * n_transitions), n_transitions, 1 / K)$p.value,
    binom.test(round(rate_down * n_transitions), n_transitions, 1 / K)$p.value
  ), 3),
  stringsAsFactors = FALSE
)

cat("\nTransition rules across consecutive trials (two-sided exact binomial):\n\n")
print(transition_summary, row.names = FALSE)

longest_run  <- function(v) max(rle(v)$lengths)
n_distinct   <- function(v) length(unique(v))

obs_longest  <- apply(obs_mat, 1, longest_run)
obs_distinct <- apply(obs_mat, 1, n_distinct)
sim_longest  <- apply(sim_mat, 1, longest_run)
sim_distinct <- apply(sim_mat, 1, n_distinct)

p_long_run   <- mean(sim_longest  >= 4)   # >= 4 identical consecutive answers
p_few_used   <- mean(sim_distinct <= 2)   # at most two of the four alternatives

cat(sprintf(
  paste0(
    "\nLongest run of identical answers: observed mean %.2f vs. %.2f under random choice\n",
    "Patients with a run of >= 4 identical answers: %d/%d (%.2f%%), expected %.2f (p = %.3g)\n",
    "Alternatives used per patient: observed mean %.2f vs. %.2f under random choice\n",
    "Patients using <= 2 of the %d alternatives: %d/%d (%.2f%%), expected %.2f (p = %.3g)\n"
  ),
  mean(obs_longest), mean(sim_longest),
  sum(obs_longest >= 4), n_obs, 100 * mean(obs_longest >= 4), n_obs * p_long_run,
  binom.test(sum(obs_longest >= 4), n_obs, p_long_run, alternative = "greater")$p.value,
  mean(obs_distinct), mean(sim_distinct),
  K, sum(obs_distinct <= 2), n_obs, 100 * mean(obs_distinct <= 2), n_obs * p_few_used,
  binom.test(sum(obs_distinct <= 2), n_obs, p_few_used, alternative = "greater")$p.value
))


# -----------------------------------------------------------------------------
# 11.5 Figure – distance to the nearest foreseeable template
# -----------------------------------------------------------------------------

sim_min_distance <- rep(L, n_sim)
for (j in seq_len(n_templates)) {
  sim_min_distance <- pmin(sim_min_distance, colSums(t(sim_mat) != tpl_mat[j, ]))
}

# Proportions per group, tabulated over the common distance range.
dist_levels <- 0:max(patient_min_distance, sim_min_distance)

prop_of <- function(x) {
  as.numeric(table(factor(x, levels = dist_levels))) / length(x)
}

df_dist <- rbind(
  data.frame(distance   = dist_levels,
             proportion = prop_of(patient_min_distance),
             source     = "Observed patients"),
  data.frame(distance   = dist_levels,
             proportion = prop_of(sim_min_distance),
             source     = "Random choice (simulated)")
)

hamming_plot_foreseeable_patterns <- 
  ggplot2::ggplot(df_dist, ggplot2::aes(x = distance, y = proportion, fill = source)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(preserve = "single"),
      width    = 0.8
    ) +
    ggplot2::scale_x_continuous(breaks = dist_levels) +
    ggplot2::scale_fill_manual(values = c("Observed patients" = "dodgerblue",
                                          "Random choice (simulated)" = "grey60")) +
    ggplot2::labs(
      title    = "Distance to the nearest foreseeable\nresponse pattern",
      subtitle = sprintf(
        "Hamming distance to the closest of %d a priori templates; \n %d patients, %d simulated random sequences",
        n_templates, n_obs, n_sim
      ),
      x    = "Hamming distance to nearest template (0 = exact match)",
      y    = "Proportion of sequences",
      fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position.inside = TRUE, legend.position = c(.3,.8), axis.title.x = element_text(size = 10))



print(hamming_plot_foreseeable_patterns)

ggplot2::ggsave(filename = "hamming_plot_foreseeable_patterns.svg", plot = hamming_plot_foreseeable_patterns, height = 8, width = 5)


# =============================================================================
# 12. Correct-answer matrix with per-patient and per-item counts
# =============================================================================

# `correct_key` is defined in Section 11.

pen_mat <- as.matrix(encoded_matrix[, pen_cols])

# 1 = correct, 0 = wrong, NA preserved
correct_mat <- sweep(
  pen_mat,
  2,
  correct_key,
  FUN = "=="
) * 1L

# Scores
row_score <- rowSums(correct_mat, na.rm = TRUE)
col_score <- colSums(correct_mat, na.rm = TRUE)


# ---------------------------------------------------------
# Annotations
# ---------------------------------------------------------

ha_right <- rowAnnotation(
  Score = anno_barplot(
    row_score,
    baseline = 0,
    gp = gpar(fill = "lightblue", col = NA)
  ),
  width = unit(3, "cm")
)

ha_top <- HeatmapAnnotation(
  `Correct answer` = anno_text(
    correct_key, rot = 0,
    gp = gpar(fontface = "bold", fontsize = 9)
  ),
  `N correct` = anno_barplot(
    col_score,
    baseline = 0,
    gp = gpar(fill = "lightblue", col = NA)
  ),
  height = unit(2, "cm")
)


# ---------------------------------------------------------
# Heatmap
# ---------------------------------------------------------

pen_names <- c("Orange", "Shoe leather", "Cinnamon", "Peppermint", "Banana", "Lemon", "Licorice", "Turpentine",
               "Garlic", "Coffee", "Apple", "Clove", "Pineapple", "Rose", "Anise", "Fish")

correct_mat_odors <- correct_mat
colnames(correct_mat_odors) <- pen_names

create_heatmap <- function() {
  ht1 <- Heatmap(
    correct_mat_odors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(4, "cm"),
    col = c("0" = "lightgray", "1" = "dodgerblue"),
    right_annotation = ha_right,
    top_annotation = ha_top,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 6),
    name = "Response",
    heatmap_legend_param = list(
      title = "Response",
      at = c(0, 1),
      labels = c("Wrong", "Correct"),
      direction = "horizontal",
      title_position = "topcenter"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(pen_mat[i, j])) {
        grid.text(pen_mat[i, j], x = x, y = y,
                  gp = gpar(fontsize = 8, col = "black"))
      }
    }
  )

  grid.newpage()
  draw(ht1,
       heatmap_legend_side = "bottom",
       padding = unit(c(25, 2, 15, 2), "mm"))

  grid::upViewport(0)

  grid.text(
    "Response patterns of anosmic patients in the Sniffin' Sticks identification test",
    x = unit(4, "mm"),
    y = unit(1, "npc") - unit(5, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 13, fontface = "plain", fontfamily = "Arial", col = "#222222")
  )
}

create_heatmap()

svg(filename = "Correct_answer_matrix.svg", width = 12, height = 13)
create_heatmap()
dev.off()

# =============================================================================
# 13. Encoded-choice heatmap with correct-key and score annotations
# =============================================================================

choice_colors <- ggplot2::alpha(ggthemes::colorblind_pal()(8)[2:5], 0.7)
# choice_colors <- viridis::turbo(4)[1:4]
names(choice_colors) <- as.character(1:4)

pen_mat_odors <- pen_mat
colnames(pen_mat_odors) <- pen_names



ha_top_choices <- HeatmapAnnotation(
  `Correct answer 1` = anno_text(
    correct_key, rot = 0,
    gp = gpar(fontface = "plain", fontsize = 9)
  ),
  `Correct answer` = as.character(correct_key),
  `N correct` = anno_barplot(
    col_score,
    baseline = 0,
    gp = gpar(fill = "lightblue", col = NA)
  ),
  col = list(`Correct answer` = choice_colors),
  show_legend = FALSE,
  height = unit(3, "cm")
)

create_choice_heatmap <- function() {
  ht2 <- Heatmap(
    pen_mat_odors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = choice_colors, 
    right_annotation = ha_right,
    top_annotation = ha_top_choices,
    row_dend_width = unit(4, "cm"),
    column_dend_height = unit(4, "cm"),
    name = "Choice position",
    heatmap_legend_param = list(
      title = "Choice position",
      direction = "horizontal",
      title_position = "topcenter"
    ),
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 6),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(pen_mat[i, j])) {
        grid.text(pen_mat[i, j], x = x, y = y,
                  gp = gpar(fontsize = 8, col = "black"))
      }
    }
  )

  grid.newpage()
  draw(ht2,
       heatmap_legend_side = "bottom",
       padding = unit(c(25, 2, 15, 2), "mm"))
  
  gb<-grid.grabExpr(draw(ht2,
                         heatmap_legend_side = "bottom",
                         padding = unit(c(25, 2, 15, 2), "mm")))
  

  grid::upViewport(0)

  grid.text(
    "Answer choices in the Sniffin' Sticks identification test",
    x = unit(4, "mm"),
    y = unit(1, "npc") - unit(5, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 13, fontface = "plain", fontfamily = "Arial", col = "#222222")
  )
}

create_choice_heatmap()

svg(filename = "Choice_matrix.svg", width = 12, height = 13)
create_choice_heatmap()
dev.off()

# =============================================================================
# 14. ABC analysis of per-item correct-answer counts
# =============================================================================

sum_per_pen_single <- apply(correct_mat, 2, sum)

ABC_sum_per_pen_single <- cABC_analysis(sum_per_pen_single, PlotIt = TRUE)
names(sum_per_pen_single[ABC_sum_per_pen_single$Aind])


# =============================================================================
# End of file
# =============================================================================

