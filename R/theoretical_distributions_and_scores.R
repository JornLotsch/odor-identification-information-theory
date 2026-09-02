# =============================================================================
# THEORETICAL DISTRIBUTIONS FOR CLASSICAL AND INFORMATIVENESS SCORES
# =============================================================================
#
# Mathematical definitions
#
# Classical score:
#   S ~ Binomial(L, 1/k)
#
# Block entropy:
#   H_m(X) = -sum_b p_b log2(p_b)
#   where p_b is the empirical frequency of overlapping blocks of length m.
#
# Informativeness:
#   I(X) = w1 H1(X) + w2 H2(X) + w3 H3(X)
#
# Default weights:
#   w1 = w2 = w3 = 1
#
# Reference distributions:
#   1. Unconstrained: each position independently sampled from 1,...,k.
#   2. Constrained: random permutation of exactly L/k copies of each symbol.
#
# Important:
#   The alphabet size k is always passed explicitly. It is never inferred
#   from the symbols observed in an individual simulated sequence.
#
# =============================================================================

# ---- External functions and parameters ---------------------------------------------------

source("globals.R")


# =============================================================================
# 0. Packages and global options
# =============================================================================

if (!requireNamespace("pbmcapply", quietly = TRUE)) {
  stop(
    "Package 'pbmcapply' is required. Install it with:\n",
    "install.packages('pbmcapply')"
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. Install it with:\n",
    "install.packages('ggplot2')"
  )
}

options(stringsAsFactors = FALSE)


# =============================================================================
# 1. Reproducible random-number helpers
# =============================================================================

# pbmclapply uses parallel random streams. To ensure that serial and parallel
# calculations generate exactly the same sequence for every replicate, each
# replicate receives its own deterministic seed.

make_replicate_seeds <- function(n_sim, seed) {
  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed)) {
    stop("'seed' must be one non-missing numeric value.")
  }
  
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
  
  sample.int(
    n = .Machine$integer.max,
    size = n_sim,
    replace = FALSE
  )
}

with_seed <- function(seed, expr) {
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(as.integer(seed))
  eval(substitute(expr), envir = parent.frame())
}


# =============================================================================
# 2. Classical chance distribution
# =============================================================================

chance_distribution <- function(L, k, alpha = c(0.05, 0.01)) {
  stopifnot(length(L) == 1L, length(k) == 1L)
  stopifnot(L >= 1, k >= 2, L == as.integer(L), k == as.integer(k))
  
  p <- 1 / k
  
  expected_score <- L * p
  variance_score <- L * p * (1 - p)
  sd_score <- sqrt(variance_score)
  
  # s_crit = min{s in 0:L : P(S >= s) <= alpha}
  # pbinom(s - 1, lower.tail=FALSE) = P(S >= s)
  s_vals         <- 0L:L
  tail_prob_ge   <- pbinom(s_vals - 1L, size = L, prob = p, lower.tail = FALSE)
  
  critical_values <- vapply(alpha, function(a) {
    eligible_scores <- s_vals[tail_prob_ge <= a]
    
    if (length(eligible_scores) == 0L) {
      NA_integer_
    } else {
      min(eligible_scores)
    }
  }, integer(1))
  
  critical_tail_probabilities <- vapply(
    critical_values,
    function(s_crit) {
      if (is.na(s_crit)) {
        NA_real_
      } else {
        pbinom(
          s_crit - 1L,
          size = L,
          prob = p,
          lower.tail = FALSE
        )
      }
    },
    numeric(1)
  )
  
  names(critical_values) <- paste0("alpha_", alpha)
  names(critical_tail_probabilities) <- paste0("alpha_", alpha)
  
  scores <- 0:L
  probs <- dbinom(scores, size = L, prob = p)
  p_ge <- pbinom(scores - 1L, size = L, prob = p, lower.tail = FALSE)
  
  list(
    L = L,
    k = k,
    p = p,
    expected_score = expected_score,
    variance_score = variance_score,
    sd_score = sd_score,
    critical_values = critical_values,
    critical_tail_probabilities = critical_tail_probabilities,
    scores = scores,
    probs = probs,
    p_ge = p_ge
  )
}


create_chance_table <- function(L, k, alpha = c(0.05, 0.01)) {
  dist <- chance_distribution(L, k, alpha = alpha)
  
  tbl <- data.frame(
    L = L,
    k = k,
    s_correct = dist$scores,
    p_exact = dist$probs,
    p_ge = dist$p_ge
  )
  
  tbl$signif <- ""
  
  for (a in sort(alpha)) {
    label <- if (a == 0.01) {
      "**"
    } else if (a == 0.05) {
      "*"
    } else {
      paste0("p<", a)
    }
    
    tbl$signif[tbl$p_ge <= a & tbl$signif == ""] <- label
  }
  
  tbl
}


plot_chance_distribution <- function(L, k, s_observed = NULL,
                                     title_prefix = "") {
  dist <- chance_distribution(L, k)
  
  df <- data.frame(score = dist$scores, prob = dist$probs)
  
  alpha_05 <- dist$critical_values["alpha_0.05"]
  alpha_01 <- dist$critical_values["alpha_0.01"]
  
  vline_df <- data.frame(
    xint  = c(dist$expected_score, alpha_05, alpha_01),
    label = c("Expected", "p < 0.05", "p < 0.01")
  )
  
  vline_df <- vline_df[is.finite(vline_df$xint), , drop = FALSE]
  
  colour_vals  <- c("Expected" = "blue", "p < 0.05" = "orange",
                    "p < 0.01" = "red",  "Observed"  = "darkgreen")
  linetype_vals <- c("Expected" = "dashed", "p < 0.05" = "dashed",
                     "p < 0.01" = "dashed", "Observed"  = "solid")
  active_breaks <- vline_df$label
  
  if (!is.null(s_observed)) {
    if (!s_observed %in% dist$scores) {
      stop("s_observed must be an integer between 0 and L.")
    }
    vline_df <- rbind(vline_df,
                      data.frame(xint = s_observed, label = "Observed"))
    active_breaks <- c(active_breaks, "Observed")
  }
  
  ggplot(df, aes(x = score, y = prob)) +
    geom_col(fill = "lightgray", colour = "gray50", width = 0.2) +
    geom_vline(
      data      = vline_df,
      aes(xintercept = xint, colour = label, linetype = label),
      linewidth = 0.8
    ) +
    scale_colour_manual(
      name   = NULL,
      values = colour_vals,
      breaks = active_breaks
    ) +
    scale_linetype_manual(
      name   = NULL,
      values = linetype_vals,
      breaks = active_breaks
    ) +
    scale_x_continuous(breaks = dist$scores) +
    labs(
      title = paste0(title_prefix, "Chance distribution (L=", L, ", k=", k, ")"),
      x     = "Number of correct answers",
      y     = "Probability"
    ) +
    theme_bw() +
    theme(legend.position = "top")
}


# =============================================================================
# 3. Entropy and informativeness functions
# =============================================================================

entropy_from_counts <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[counts > 0]
  
  if (length(counts) == 0L) {
    return(0)
  }
  
  probabilities <- counts / sum(counts)
  
  -sum(probabilities * log2(probabilities))
}


block_entropy_general <- function(types, m = 1L, k) {
  types <- as.integer(types)
  L <- length(types)
  
  if (length(k) != 1L || k < 2L || k != as.integer(k)) {
    stop("'k' must be one integer greater than or equal to 2.")
  }
  
  if (m < 1L || m != as.integer(m)) {
    stop("'m' must be a positive integer.")
  }
  
  if (m > L) {
    stop("Block length m cannot exceed sequence length.")
  }
  
  if (any(types < 1L | types > k)) {
    stop("All sequence symbols must lie in 1,...,k.")
  }
  
  n_blocks <- L - m + 1L
  
  block_keys <- vapply(
    seq_len(n_blocks),
    function(i) {
      paste(types[i:(i + m - 1L)], collapse = "_")
    },
    character(1)
  )
  
  counts <- table(block_keys)
  
  entropy_from_counts(counts)
}


informativeness_score_general <- function(
    types,
    k,
    weights = c(1, 1, 1)) {
  
  types <- as.integer(types)
  L <- length(types)
  
  if (L < 3L) {
    stop("The current I(X) definition requires L >= 3.")
  }
  
  if (length(weights) != 3L || any(!is.finite(weights))) {
    stop("'weights' must contain three finite values.")
  }
  
  H1 <- block_entropy_general(types, m = 1L, k = k)
  H2 <- block_entropy_general(types, m = 2L, k = k)
  H3 <- block_entropy_general(types, m = 3L, k = k)
  
  I_total <- sum(weights * c(H1, H2, H3))
  I_per_trial <- I_total / L
  
  list(
    H1 = H1,
    H2 = H2,
    H3 = H3,
    I_total = I_total,
    I_per_trial = I_per_trial,
    score = I_total
  )
}


validate_equal_frequency_sequence <- function(types, L, k) {
  types <- as.integer(types)
  
  if (length(types) != L) {
    stop("Sequence length does not equal L.")
  }
  
  if (any(types < 1L | types > k)) {
    stop("All sequence symbols must lie in 1,...,k.")
  }
  
  if (L %% k != 0L) {
    stop("Equal-frequency constraint requires L %% k == 0.")
  }
  
  expected_count <- L / k
  observed_counts <- tabulate(types, nbins = k)
  
  if (!all(observed_counts == expected_count)) {
    stop(
      "The sequence does not contain exactly L/k copies of every symbol."
    )
  }
  
  TRUE
}


# =============================================================================
# 4. Deterministic sequence generation
# =============================================================================

generate_unconstrained_sequence <- function(L, k, seed) {
  with_seed(
    seed,
    sample.int(k, size = L, replace = TRUE)
  )
}


generate_constrained_sequence <- function(L, k, seed) {
  if (L %% k != 0L) {
    stop("Constrained sequences require L %% k == 0.")
  }
  
  with_seed(
    seed,
    sample(rep(seq_len(k), each = L / k), replace = FALSE)
  )
}


# =============================================================================
# 5. Constrained reference maximum
# =============================================================================

compute_I_max <- function(
    L,
    k,
    n_sim = 100000L,
    seed = 123L,
    mc.cores = 1L,
    weights = c(1, 1, 1)) {
  
  if (L %% k != 0L) {
    stop("Constrained designs require L %% k == 0.")
  }
  
  seeds <- make_replicate_seeds(n_sim, seed)
  
  score_one <- function(i) {
    X_sim <- generate_constrained_sequence(L, k, seeds[i])
    
    informativeness_score_general(
      X_sim,
      k = k,
      weights = weights
    )$score
  }
  
  I_values <- pbmclapply(
    seq_len(n_sim),
    score_one,
    mc.cores = mc.cores,
    mc.set.seed = FALSE
  )
  
  I_values <- unlist(I_values, use.names = FALSE)
  
  list(
    I_max = max(I_values),
    I_values = I_values,
    n_sim = n_sim,
    L = L,
    k = k,
    seed = seed,
    weights = weights
  )
}


# =============================================================================
# 6. Exact enumeration for small designs
# =============================================================================

enumerate_all_sequences <- function(L, k, weights = c(1, 1, 1)) {
  
  if (L < 3L) {
    stop("The current I(X) definition requires L >= 3.")
  }
  
  all_sequences <- expand.grid(
    rep(list(seq_len(k)), L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  sequence_matrix <- as.matrix(all_sequences)
  
  I_values <- apply(sequence_matrix, 1, function(X) {
    informativeness_score_general(
      types = as.integer(X),
      k = k,
      weights = weights
    )$score
  })
  
  list(
    L = L,
    k = k,
    sequences = sequence_matrix,
    I_values = I_values,
    I_max_feasible = max(I_values),
    n_sequences = nrow(sequence_matrix),
    exact = TRUE
  )
}


get_feasible_reference <- function(
    L,
    k,
    weights = c(1, 1, 1),
    max_sequences = 1e6) {
  
  n_possible <- k^L
  
  if (n_possible > max_sequences) {
    return(list(
      L = L,
      k = k,
      I_max_feasible = NA_real_,
      I_values = NULL,
      n_sequences = n_possible,
      exact = FALSE,
      reason = paste0(
        "Sequence space too large for exhaustive enumeration: ",
        k, "^", L, " = ", format(n_possible, scientific = FALSE)
      )
    ))
  }
  
  enumerate_all_sequences(
    L = L,
    k = k,
    weights = weights
  )
}


# =============================================================================
# 6. Reference distributions for I(X)
# =============================================================================

summarize_numeric_distribution <- function(x) {
  qs <- quantile(
    x,
    probs = c(0.50, 0.90, 0.95, 0.99),
    names = TRUE,
    type = 7
  )
  
  list(
    mean = mean(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    quantiles = qs
  )
}


compute_score_distributions <- function(
    L,
    k,
    n_sim = 100000L,
    seed = 123L,
    n_sim_max = 100000L,
    seed_max = seed + 1000L,
    mc.cores = 1L,
    weights = c(1, 1, 1)) {
  
  if (L < 3L) {
    stop("The current score definition requires L >= 3.")
  }
  
  if (L %% k != 0L) {
    warning(
      "L is not divisible by k; constrained distributions will be NULL."
    )
  }
  
  # Fixed normalization value for relative efficiency.
  # This is the maximum I observed in the constrained reference sample.
  if (L %% k == 0L) {
    I_max_reference <- compute_I_max(
      L = L,
      k = k,
      n_sim = n_sim_max,
      seed = seed_max,
      mc.cores = mc.cores,
      weights = weights
    )$I_max
  } else {
    I_max_reference <- NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # Unconstrained reference distribution
  # ---------------------------------------------------------------------------
  
  uncon_seeds <- make_replicate_seeds(n_sim, seed)
  
  score_unconstrained <- function(i) {
    X_sim <- generate_unconstrained_sequence(L, k, uncon_seeds[i])
    
    I_value <- informativeness_score_general(
      X_sim,
      k = k,
      weights = weights
    )$score
    
    I_per_trial <- I_value / L
    
    c(
      I = I_value,
      I_per_trial = I_per_trial
    )
  }
  
  uncon_res <- pbmclapply(
    seq_len(n_sim),
    score_unconstrained,
    mc.cores = mc.cores,
    mc.set.seed = FALSE
  )
  
  uncon_matrix <- do.call(rbind, uncon_res)
  
  I_unconstrained <- uncon_matrix[, "I"]
  
  # ---------------------------------------------------------------------------
  # Constrained reference distribution
  # ---------------------------------------------------------------------------
  
  I_constrained <- NULL
  
  if (L %% k == 0L) {
    constrained_seeds <- make_replicate_seeds(n_sim, seed + 1L)
    
    score_constrained <- function(i) {
      X_sim <- generate_constrained_sequence(
        L,
        k,
        constrained_seeds[i]
      )
      
      I_value <- informativeness_score_general(
        X_sim,
        k = k,
        weights = weights
      )$score
      
      c(
        I = I_value,
        I_per_trial = I_value / L
      )
    }
    
    con_res <- pbmclapply(
      seq_len(n_sim),
      score_constrained,
      mc.cores = mc.cores,
      mc.set.seed = FALSE
    )
    
    con_matrix <- do.call(rbind, con_res)
    
    I_constrained <- con_matrix[, "I"]
  }
  
  list(
    L = L,
    k = k,
    weights = weights,
    n_sim = n_sim,
    n_sim_max = n_sim_max,
    seed = seed,
    seed_max = seed_max,
    I_max_reference = I_max_reference,
    I_unconstrained = I_unconstrained,
    I_constrained = I_constrained,
    uncon_seeds = uncon_seeds,
    constrained_seeds = if (L %% k == 0L) constrained_seeds else NULL,
    summary_I_unconstrained = summarize_numeric_distribution(
      I_unconstrained
    ),
    summary_I_constrained = if (!is.null(I_constrained)) {
      summarize_numeric_distribution(I_constrained)
    } else {
      NULL
    }
  )
}


# =============================================================================
# 7. Relative efficiency E(X)
# =============================================================================

# E(X) = I(X) / I_max^(c)
#
# I_max^(c) is the maximum informativeness observed in the constrained
# Monte Carlo sample, stored as dist_obj$I_max_reference.

compute_E <- function(I_obs, dist_obj) {
  stopifnot(is.numeric(I_obs), length(I_obs) >= 1L)
  
  I_max <- dist_obj$I_max_reference
  
  if (!is.finite(I_max) || I_max <= 0) {
    stop("I_max_reference is not a positive finite value.")
  }
  
  I_obs / I_max
}


# =============================================================================
# Feasible-design efficiency for exhaustively enumerated small designs
# =============================================================================

# E_feasible(X) = I(X) / I_max_feasible
#
# This differs from E(X), which is defined relative to the constrained,
# equal-frequency reference maximum.

compute_E_feasible <- function(I_obs, feasible_ref) {
  stopifnot(is.numeric(I_obs), length(I_obs) >= 1L)
  
  I_max_feasible <- feasible_ref$I_max_feasible
  
  if (!is.finite(I_max_feasible) || I_max_feasible <= 0) {
    return(NA_real_)
  }
  
  I_obs / I_max_feasible
}


# =============================================================================
# 8. Empirical upper-tail tables for I(X)
# =============================================================================

create_score_tail_table <- function(
    score_vec,
    score_name = "Score",
    alpha = c(0.05, 0.01)) {
  
  score_vec <- score_vec[is.finite(score_vec)]
  
  if (length(score_vec) == 0L) {
    stop("score_vec contains no finite values.")
  }
  
  score_vec <- sort(score_vec)
  n <- length(score_vec)
  
  # Use the distinct observed simulated values as thresholds.
  thresholds <- sort(unique(score_vec))
  
  p_exact <- vapply(
    thresholds,
    function(t) mean(score_vec == t),
    numeric(1)
  )
  
  p_ge <- vapply(
    thresholds,
    function(t) mean(score_vec >= t),
    numeric(1)
  )
  
  tbl <- data.frame(
    score = thresholds,
    p_exact = p_exact,
    p_ge = p_ge
  )
  
  tbl$signif <- ""
  tbl$signif[tbl$p_ge <= 0.01] <- "**"
  tbl$signif[tbl$p_ge <= 0.05 & tbl$p_ge > 0.01] <- "*"
  
  tbl$score_name <- score_name
  tbl$alpha_05_threshold <- as.numeric(
    quantile(score_vec, 0.95, type = 7)
  )
  tbl$alpha_01_threshold <- as.numeric(
    quantile(score_vec, 0.99, type = 7)
  )
  
  tbl
}


# =============================================================================
# 9. Summary tables for I(X)
# =============================================================================

create_score_summary_table <- function(dist_obj, score_type = "I") {
  if (!score_type %in% c("I")) {
    stop("score_type must be 'I'.")
  }
  
  get_summary <- function(design) {
    if (design == "Unconstrained") {
      dist_obj$summary_I_unconstrained
    } else {
      dist_obj$summary_I_constrained
    }
  }
  
  designs <- c("Unconstrained", "Constrained")
  
  rows <- lapply(designs, function(design) {
    s <- get_summary(design)
    
    if (is.null(s)) {
      return(data.frame(
        L = dist_obj$L,
        k = dist_obj$k,
        score_type = score_type,
        design = design,
        mean = NA_real_,
        sd = NA_real_,
        min = NA_real_,
        max = NA_real_,
        median = NA_real_,
        p90 = NA_real_,
        p95 = NA_real_,
        p99 = NA_real_
      ))
    }
    
    data.frame(
      L = dist_obj$L,
      k = dist_obj$k,
      score_type = score_type,
      design = design,
      mean = s$mean,
      sd = s$sd,
      min = s$min,
      max = s$max,
      median = unname(s$quantiles["50%"]),
      p90 = unname(s$quantiles["90%"]),
      p95 = unname(s$quantiles["95%"]),
      p99 = unname(s$quantiles["99%"])
    )
  })
  
  do.call(rbind, rows)
}


# =============================================================================
# 10. Plots with exact score-axis thresholds
# =============================================================================

plot_score_distribution <- function(
    dist_obj,
    score_type = "I",
    title_prefix = "",
    observed_score = NULL,
    show_thresholds = TRUE) {
  
  if (!score_type %in% c("I")) stop("score_type must be 'I'.")
  
  x_uncon     <- dist_obj$I_unconstrained
  x_con       <- dist_obj$I_constrained
  score_label <- "Informativeness score I(X) [bit]"
  
  x_uncon <- x_uncon[is.finite(x_uncon)]
  if (!is.null(x_con)) x_con <- x_con[is.finite(x_con)]
  
  dens_uncon <- density(x_uncon)
  dens_df <- data.frame(x = dens_uncon$x, y = dens_uncon$y, grp = "Unconstrained")
  
  if (!is.null(x_con)) {
    dens_con <- density(x_con)
    dens_df  <- rbind(dens_df,
                      data.frame(x = dens_con$x, y = dens_con$y, grp = "Constrained"))
  }
  
  # All possible colour/linetype values across density lines and vlines.
  # The threshold lines reference the unconstrained distribution.
  colour_vals <- c(
    "Unconstrained"                  = "gray30",
    "Constrained"                    = "blue",
    "Unconstrained mean"             = "blue",
    "Unconstrained 95th percentile"  = "orange",
    "Unconstrained 99th percentile"  = "red",
    "Observed"                       = "darkgreen"
  )
  linetype_vals <- c(
    "Unconstrained"                  = "solid",
    "Constrained"                    = "dashed",
    "Unconstrained mean"             = "dashed",
    "Unconstrained 95th percentile"  = "dashed",
    "Unconstrained 99th percentile"  = "dashed",
    "Observed"                       = "solid"
  )
  
  active_breaks <- if (!is.null(x_con)) c("Unconstrained", "Constrained") else "Unconstrained"
  
  p <- ggplot() +
    geom_line(
      data      = dens_df,
      aes(x = x, y = y, colour = grp, linetype = grp),
      linewidth = 0.9
    )
  
  if (show_thresholds) {
    threshold_mean <- mean(x_uncon)
    threshold_95   <- as.numeric(quantile(x_uncon, 0.95, type = 7))
    threshold_99   <- as.numeric(quantile(x_uncon, 0.99, type = 7))
    
    vline_df <- data.frame(
      xint = c(threshold_mean, threshold_95, threshold_99),
      lbl  = c("Unconstrained mean",
               "Unconstrained 95th percentile",
               "Unconstrained 99th percentile")
    )
    active_breaks <- c(active_breaks, vline_df$lbl)
    
    p <- p + geom_vline(
      data      = vline_df,
      aes(xintercept = xint, colour = lbl, linetype = lbl),
      linewidth = 0.8
    )
  }
  
  if (!is.null(observed_score)) {
    obs_df <- data.frame(xint = observed_score, lbl = "Observed")
    active_breaks <- c(active_breaks, "Observed")
    p <- p + geom_vline(
      data      = obs_df,
      aes(xintercept = xint, colour = lbl, linetype = lbl),
      linewidth = 1.1
    )
  }
  
  p +
    scale_colour_manual(
      name   = NULL,
      values = colour_vals,
      breaks = active_breaks
    ) +
    scale_linetype_manual(
      name   = NULL,
      values = linetype_vals,
      breaks = active_breaks
    ) +
    labs(
      title = paste0(title_prefix, "Distribution of ", score_label,
                     " (L=", dist_obj$L, ", k=", dist_obj$k, ")"),
      x = score_label,
      y = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "top", legend.text = element_text(size = 8))
}


plot_score_reference <- function(
    score_vec,
    score_name = "Score",
    observed_score = NULL,
    title_prefix = "") {
  
  score_vec <- score_vec[is.finite(score_vec)]
  if (length(score_vec) == 0L) stop("No finite values available for plotting.")
  
  dens    <- density(score_vec)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  mean_score <- mean(score_vec)
  p95_score  <- as.numeric(quantile(score_vec, 0.95, type = 7))
  p99_score  <- as.numeric(quantile(score_vec, 0.99, type = 7))
  
  vline_df <- data.frame(
    xint  = c(mean_score, p95_score, p99_score),
    label = c("Mean", "95th percentile", "99th percentile")
  )
  
  colour_vals  <- c("Mean" = "blue", "95th percentile" = "orange",
                    "99th percentile" = "red", "Observed" = "darkgreen")
  linetype_vals <- c("Mean" = "dashed", "95th percentile" = "dashed",
                     "99th percentile" = "dashed", "Observed" = "solid")
  active_breaks <- c("Mean", "95th percentile", "99th percentile")
  
  if (!is.null(observed_score)) {
    vline_df <- rbind(vline_df,
                      data.frame(xint = observed_score, label = "Observed"))
    active_breaks <- c(active_breaks, "Observed")
  }
  
  ggplot() +
    geom_line(data = dens_df, aes(x = x, y = y),
              colour = "gray30", linewidth = 0.9) +
    geom_vline(
      data      = vline_df,
      aes(xintercept = xint, colour = label, linetype = label),
      linewidth = 0.8
    ) +
    scale_colour_manual(
      name   = NULL,
      values = colour_vals,
      breaks = active_breaks
    ) +
    scale_linetype_manual(
      name   = NULL,
      values = linetype_vals,
      breaks = active_breaks
    ) +
    labs(
      title = paste0(title_prefix, "Distribution of ", score_name),
      x     = score_name,
      y     = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "top", legend.text = element_text(size = 8))
}


# =============================================================================
# 11. Complete analysis and output for one design
# =============================================================================

analyze_score_distributions <- function(
    dist_obj,
    output_prefix,
    test_name = "",
    observed_I = NULL,
    save_plots = TRUE) {
  
  L <- dist_obj$L
  k <- dist_obj$k
  
  cat("\n=============================================================\n")
  cat("Analysis:", test_name, "(L =", L, ", k =", k, ")\n")
  cat("=============================================================\n")
  
  # Summary tables
  summary_I <- create_score_summary_table(dist_obj, score_type = "I")
  
  write.csv(summary_I, paste0(output_prefix, "_summary_I.csv"), row.names = FALSE)
  
  print(summary_I, row.names = FALSE)
  
  # Combined density plots (unconstrained + constrained overlaid)
  plots <- list()
  
  for (st in c("I")) {
    p <- plot_score_distribution(dist_obj, score_type = st,
                                 title_prefix = paste0(test_name, ": "))
    print(p)
    plots[[paste0("combined_", st)]] <- p
    
    if (save_plots) {
      ggsave(
        filename = paste0(output_prefix, "_combined_", st, "_density.svg"),
        plot     = p,
        width    = 10,
        height   = 7.8,
        dpi      = 180
      )
    }
  }
  
  # Tail tables and per-distribution density plots
  score_sets <- list(
    I_unconstrained = list(
      values = dist_obj$I_unconstrained,
      observed = observed_I
    ),
    I_constrained = list(
      values = dist_obj$I_constrained,
      observed = observed_I
    )
  )
  
  for (set_name in names(score_sets)) {
    values <- score_sets[[set_name]]$values
    
    if (is.null(values)) {
      next
    }
    
    values <- values[is.finite(values)]
    
    if (length(values) == 0L) {
      next
    }
    
    if (grepl("^I_", set_name)) {
      score_name <- paste0(
        "I(X) - ",
        ifelse(grepl("unconstrained", set_name),
               "unconstrained", "constrained"),
        " [bit]"
      )
    } else {
      stop("Unexpected score set name.")
    }
    
    tail_table <- create_score_tail_table(
      values,
      score_name = score_name
    )
    
    write.csv(
      tail_table,
      paste0(output_prefix, "_", set_name, "_tail_table.csv"),
      row.names = FALSE
    )
    
    p_ref <- plot_score_reference(
      values,
      score_name     = score_name,
      observed_score = score_sets[[set_name]]$observed,
      title_prefix   = paste0(test_name, ": ")
    )
    print(p_ref)
    plots[[set_name]] <- p_ref
    
    if (save_plots) {
      ggsave(
        filename = paste0(output_prefix, "_", set_name, "_density.svg"),
        plot     = p_ref,
        width    = 10,
        height   = 7.8,
        dpi      = 180
      )
    }
  }
  
  invisible(
    list(
      summary_I = summary_I,
      plots     = plots
    )
  )
}


# =============================================================================
# 12. Top-scoring sequence extraction
# =============================================================================

# For each design type (unconstrained / constrained), the score vector and
# seed vector are parallel: score[i] was produced by seed[i]. The function
# exports up to top_n distinct sequences, ordered by decreasing I(X).

save_top_sequences <- function(dist_obj, output_prefix, top_n = 10L) {
  L <- dist_obj$L
  k <- dist_obj$k
  
  configs <- list(
    list(
      label  = "unconstrained",
      I_vec  = dist_obj$I_unconstrained,
      seeds  = dist_obj$uncon_seeds,
      gen_fn = generate_unconstrained_sequence
    ),
    list(
      label  = "constrained",
      I_vec  = dist_obj$I_constrained,
      seeds  = dist_obj$constrained_seeds,
      gen_fn = generate_constrained_sequence
    )
  )
  
  all_lines <- character(0)
  
  for (cfg in configs) {
    if (is.null(cfg$I_vec) || is.null(cfg$seeds)) {
      next
    }
    
    ordered_idx <- order(cfg$I_vec, decreasing = TRUE)
    
    unique_idx <- integer(0)
    seen_sequences <- character(0)
    
    for (i in ordered_idx) {
      seq_i <- cfg$gen_fn(L, k, cfg$seeds[i])
      sequence_key <- paste(seq_i, collapse = "")
      
      if (!sequence_key %in% seen_sequences) {
        unique_idx <- c(unique_idx, i)
        seen_sequences <- c(seen_sequences, sequence_key)
      }
      
      if (length(unique_idx) >= top_n) {
        break
      }
    }
    
    n_top <- length(unique_idx)
    
    header <- sprintf(
      "--- Top %d unique %s sequences by I score (L=%d, k=%d) ---",
      n_top, cfg$label, L, k
    )
    cat(header, "\n")
    all_lines <- c(all_lines, header)
    
    colheader <- sprintf(
      "%-6s %-10s %-12s %s",
      "Rank", "I", "seed", "sequence"
    )
    cat(colheader, "\n")
    all_lines <- c(all_lines, colheader)
    
    for (rank in seq_along(unique_idx)) {
      i <- unique_idx[rank]
      seq_i <- cfg$gen_fn(L, k, cfg$seeds[i])
      I_i <- cfg$I_vec[i]
      
      line <- sprintf(
        "%-6d %-10.4f %-12d %s",
        rank,
        I_i,
        cfg$seeds[i],
        paste(seq_i, collapse = "")
      )
      
      cat(line, "\n")
      all_lines <- c(all_lines, line)
    }
    
    all_lines <- c(all_lines, "")
    cat("\n")
  }
  
  out_file <- paste0(output_prefix, "_top_sequences.txt")
  writeLines(all_lines, out_file)
  
  cat("Top sequences saved to:", out_file, "\n")
  invisible(all_lines)
}


# =============================================================================
# 13. Run all designs
# =============================================================================

# Set this to 1 for exact reproducibility across platforms and core counts.
# Increase it after validation if desired.
MC_CORES <- 1L

N_SIM <- 100000L
N_SIM_MAX <- 100000L

WEIGHTS <- c(1, 1, 1)

# ----------------------------------------------------------------------------- 
# Analyzed distributions
# -----------------------------------------------------------------------------

designs <- list(
  L16_k4 = list(L = 16L, k = 4L),
  L18_k6 = list(L = 18L, k = 6L),
  L12_k6 = list(L = 12L, k = 6L),
  L16_k6 = list(L = 16L, k = 6L),
  L24_k6 = list(L = 24L, k = 6L),
  L40_k4 = list(L = 40L, k = 4L),
  L3_k4  = list(L = 3L,  k = 4L),
  L4_k4  = list(L = 4L,  k = 4L)
)

chance_plots <- list()

for (design_name in names(designs)) {
  design <- designs[[design_name]]
  
  chance_obj <- chance_distribution(
    L = design$L,
    k = design$k
  )
  
  chance_table <- create_chance_table(
    L = design$L,
    k = design$k
  )
  
  write.csv(
    chance_table,
    paste0("chance_table_", design_name, ".csv"),
    row.names = FALSE
  )
  
  p_chance <- plot_chance_distribution(L = design$L, k = design$k)
  print(p_chance)
  chance_plots[[design_name]] <- p_chance
  
  ggsave(
    filename = paste0("chance_distribution_", design_name, ".svg"),
    plot     = p_chance,
    width    = 10,
    height   = 7.8,
    dpi      = 180
  )
}

# ----------------------------------------------------------------------------- 
# Information-theoretic distributions
# -----------------------------------------------------------------------------

distribution_objects  <- list()
analysis_results      <- list()

for (design_name in names(designs)) {
  design <- designs[[design_name]]
  
  cat(
    "\n\nRunning score distributions for ",
    design_name,
    " (L=",
    design$L,
    ", k=",
    design$k,
    ")\n",
    sep = ""
  )
  
  dist_obj <- compute_score_distributions(
    L = design$L,
    k = design$k,
    n_sim = N_SIM,
    seed = 123L,
    n_sim_max = N_SIM_MAX,
    seed_max = 1123L,
    mc.cores = MC_CORES,
    weights = WEIGHTS
  )
  
  distribution_objects[[design_name]] <- dist_obj
  
  analysis_results[[design_name]] <- analyze_score_distributions(
    dist_obj = dist_obj,
    output_prefix = paste0("scores_", design_name),
    test_name = design_name,
    observed_I = NULL,
    save_plots = TRUE
  )
  
  save_top_sequences(
    dist_obj      = dist_obj,
    output_prefix = paste0("scores_", design_name),
    top_n         = 10L
  )
}


# =============================================================================
# 14. Summary table: mean and 95% reference interval of I(X)
# =============================================================================

# The 95% interval is a reference interval for simulated sequences:
# 2.5th to 97.5th percentile of the corresponding I(X) distribution.
# It is not a confidence interval around the Monte Carlo mean.

summarize_I_reference_distribution <- function(
    design_name,
    dist_obj,
    distribution_type = c("Unconstrained", "Constrained")) {
  
  distribution_type <- match.arg(distribution_type)
  
  I_values <- if (distribution_type == "Unconstrained") {
    dist_obj$I_unconstrained
  } else {
    dist_obj$I_constrained
  }
  
  if (is.null(I_values)) {
    return(data.frame(
      Design = design_name,
      L = dist_obj$L,
      k = dist_obj$k,
      Reference_set = distribution_type,
      N_sequences = NA_integer_,
      Mean_I_bits = NA_real_,
      SD_I_bits = NA_real_,
      CI_95_lower_bits = NA_real_,
      CI_95_upper_bits = NA_real_,
      Max_I_bits = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  I_values <- I_values[is.finite(I_values)]
  
  if (length(I_values) == 0L) {
    return(data.frame(
      Design = design_name,
      L = dist_obj$L,
      k = dist_obj$k,
      Reference_set = distribution_type,
      N_sequences = NA_integer_,
      Mean_I_bits = NA_real_,
      SD_I_bits = NA_real_,
      CI_95_lower_bits = NA_real_,
      CI_95_upper_bits = NA_real_,
      Max_I_bits = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  I_interval <- quantile(
    I_values,
    probs = c(0.025, 0.975),
    type = 7,
    names = FALSE
  )
  
  data.frame(
    Design = design_name,
    L = dist_obj$L,
    k = dist_obj$k,
    Reference_set = distribution_type,
    N_sequences = length(I_values),
    Mean_I_bits = mean(I_values),
    SD_I_bits = sd(I_values),
    CI_95_lower_bits = I_interval[1],
    CI_95_upper_bits = I_interval[2],
    Max_I_bits = max(I_values),
    stringsAsFactors = FALSE
  )
}


# One row for unconstrained and, where available, constrained distributions.
I_reference_summary_table <- do.call(
  rbind,
  unlist(
    lapply(names(distribution_objects), function(design_name) {
      
      dist_obj <- distribution_objects[[design_name]]
      
      list(
        summarize_I_reference_distribution(
          design_name = design_name,
          dist_obj = dist_obj,
          distribution_type = "Unconstrained"
        ),
        summarize_I_reference_distribution(
          design_name = design_name,
          dist_obj = dist_obj,
          distribution_type = "Constrained"
        )
      )
    }),
    recursive = FALSE
  )
)

row.names(I_reference_summary_table) <- NULL


# Add a Word-friendly descriptive column.
I_reference_summary_table$`I(X), mean [95% reference interval]` <- ifelse(
  is.finite(I_reference_summary_table$Mean_I_bits),
  sprintf(
    "%.2f [%.2f–%.2f]",
    I_reference_summary_table$Mean_I_bits,
    I_reference_summary_table$CI_95_lower_bits,
    I_reference_summary_table$CI_95_upper_bits
  ),
  NA_character_
)

I_reference_summary_table$`Maximum I(X) (bits)` <- ifelse(
  is.finite(I_reference_summary_table$Max_I_bits),
  sprintf(
    "%.4f",
    I_reference_summary_table$Max_I_bits
  ),
  NA_character_
)


# Round numeric values in the long output table.
numeric_columns_I_reference <- vapply(
  I_reference_summary_table,
  is.numeric,
  logical(1)
)

I_reference_summary_table[numeric_columns_I_reference] <- lapply(
  I_reference_summary_table[numeric_columns_I_reference],
  round,
  digits = 4
)


cat("\n=============================================================\n")
cat("Reference-distribution summary: I(X) mean and 95% interval\n")
cat("=============================================================\n")

print(I_reference_summary_table, row.names = FALSE)


# Long CSV: includes mean, SD, lower and upper reference limits.
write.csv(
  I_reference_summary_table,
  file = "I_reference_distribution_summary_full.csv",
  row.names = FALSE,
  na = ""
)


# Compact publication-oriented table.
I_reference_summary_compact <- I_reference_summary_table[, c(
  "Design",
  "L",
  "k",
  "Reference_set",
  "I(X), mean [95% reference interval]",
  "Maximum I(X) (bits)"
)]

print(I_reference_summary_compact, row.names = FALSE)

write.csv(
  I_reference_summary_compact,
  file = "I_reference_distribution_summary_compact.csv",
  row.names = FALSE,
  na = ""
)



# =============================================================================
# 14. Create publication plot
# =============================================================================

combined_distributions_plot <-
  (
    (
      chance_plots$L16_k4 +
        theme_plot() +
        theme(
          legend.position.inside = TRUE,
          legend.position = c(.8, .95)
        )
    ) |
    (
      (
        analysis_results$L16_k4$plots$I_constrained +
          theme_plot() +
          theme(
            legend.position.inside = TRUE,
            legend.position = c(.2, .9)
          ) +
          xlim(7, 9.7)
      ) /
      (
        analysis_results$L16_k4$plots$I_unconstrained +
          theme_plot() +
          theme(legend.position = "none") +
          xlim(7, 9.7)
      )
    )
  ) +
  plot_annotation(
    title = "Chance and Reference Distributions of Sequence Informativeness",
    subtitle = "16-trial, 4-alternative forced-choice odor identification test",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0, face = "plain"),
      plot.subtitle = element_text(size = 11, hjust = 0),
      plot.tag = element_text(face = "bold", size = 16)
    )
  )




print(combined_distributions_plot)

ggsave(filename = "combined_distributions_plot.svg", plot = combined_distributions_plot, height = 10, width = 12)


# =============================================================================
# 15. Analysis of published, modified, and candidate sequences
# =============================================================================

# -----------------------------------------------------------------------------
# Enter sequences here
#
# Each list element must contain:
#   test_name : label displayed in Table 1
#   reference : left empty for manual insertion in Word
#   L         : sequence length
#   k         : number of response alternatives
#   sequence  : vector of correct-answer positions, coded 1,...,k
# -----------------------------------------------------------------------------

sequence_list <- list(
  
  Sniffin_Sticks_published = list(
    test_name = "Sniffin' Sticks published",
    reference = "",
    L = 16L,
    k = 4L,
    sequence = c(
      4, 3, 2, 1, 1, 2, 3, 4, 4, 3, 2, 1, 1, 2, 3, 4
    )
  ),
  
  Sniffin_Sticks_constrained_optimum = list(
    test_name = "Sniffin' Sticks constrained optimum",
    reference = "",
    L = 16L,
    k = 4L,
    sequence = c(
      2, 3, 3, 4, 1, 1, 4, 2, 1, 3, 1, 2, 4, 4, 3, 2
    )
  ),
  
  UPSIT = list(
    test_name = "UPSIT",
    reference = "",
    L = 40L,
    k = 4L,
    sequence = c(2, 2, 4, 4, 3, 2, 1, 2, 3, 2, 
                 3, 2, 1, 4, 2, 4, 1, 1, 2, 3, 
                 1, 1, 2, 1, 2, 3, 4, 2, 2, 4, 
                 4, 3, 3, 1, 4, 4, 1, 4, 2, 1
    )
  ),
  
  UPSIT_constrained_optimum  = list(
    test_name = "UPSIT constrained optimum",
    reference = "",
    L = 40L,
    k = 4L,
    sequence = c(3, 2, 4, 3, 3, 1, 3, 2, 2, 3, 4, 4, 2, 2, 2, 1, 3, 1, 2, 3, 
                 3, 3, 4, 1, 4, 2, 1, 4, 3, 1, 1, 2, 4, 4, 4, 1, 2, 1, 1, 4)
  ),
  
  Short_olfactory_test = list(
    test_name = "Short olfactory test",
    reference = "",
    L = 3L,
    k = 4L,
    sequence = c(
      2,3,4
    )
  ),
  
  New_design_L16_k6 = list(
    test_name = "K6-Modified Sniffin' Sticks, L=16",
    reference = "",
    L = 16L,
    k = 6L,
    sequence = c(
      6, 5, 3, 3, 2, 6, 2, 2, 4, 1, 3, 6, 4, 5, 1, 5 
    )
  ),
  
  New_design_L12_k6 = list(
    test_name = "K6-Modified Sniffin' Sticks, L=12",
    reference = "",
    L = 12L,
    k = 6L,
    sequence = c(
      3, 4, 1, 1, 5, 3, 6, 6, 5, 4, 2, 2
    )
  ),
  
  New_design_L18_k6 = list(
    test_name = "K6-Modified Sniffin' Sticks, L=18",
    reference = "",
    L = 18L,
    k = 6L,
    sequence = c(
      3, 1, 5, 5, 1, 3, 2, 2, 4, 3, 4, 1, 6, 6, 2, 5, 4, 6  
    )
  )
  
  
  
  
  # ---------------------------------------------------------------------------
)


# =============================================================================
# 16. Helper function: analyse one sequence
# =============================================================================

analyse_sequence <- function(sequence_info,
                             distribution_objects,
                             weights = WEIGHTS) {
  
  test_name <- sequence_info$test_name
  reference <- sequence_info$reference
  L <- as.integer(sequence_info$L)
  k <- as.integer(sequence_info$k)
  X <- as.integer(sequence_info$sequence)
  
  # Validate the supplied sequence.
  if (length(X) != L) {
    stop(
      "Sequence length does not match L for: ",
      test_name,
      ". Expected L = ",
      L,
      ", but sequence has length ",
      length(X),
      "."
    )
  }
  
  if (any(X < 1L | X > k)) {
    stop(
      "Sequence values must be integers from 1 to k for: ",
      test_name
    )
  }
  
  # The distribution object must have been generated earlier in the main script.
  distribution_key <- paste0("L", L, "_k", k)
  
  if (!distribution_key %in% names(distribution_objects)) {
    stop(
      "No reference distribution found for ",
      test_name,
      ". Expected distribution_objects[['",
      distribution_key,
      "']]. Add this L/k combination to 'designs' and rerun the simulations."
    )
  }
  
  dist_obj <- distribution_objects[[distribution_key]]
  
  # Entropy and total informativeness score.
  score_obj <- informativeness_score_general(
    X,
    k = k,
    weights = weights
  )
  
  observed_H1 <- score_obj$H1
  observed_H2 <- score_obj$H2
  observed_H3 <- score_obj$H3
  observed_I <- score_obj$score
  observed_I_per_trial <- score_obj$I_per_trial
  
  # Position counts are reported in fixed order: positions 1,...,k.
  position_counts <- tabulate(X, nbins = k)
  
  # Equal-frequency criterion.
  is_constrained_obs <- (L %% k == 0L) &&
    all(position_counts == L %/% k)
  
  reference_type <- if (is_constrained_obs) {
    "Constrained"
  } else {
    "Unconstrained"
  }
  
  # Select the design-appropriate empirical I(X) reference distribution.
  I_reference <- if (is_constrained_obs) {
    dist_obj$I_constrained
  } else {
    dist_obj$I_unconstrained
  }
  
  if (is.null(I_reference)) {
    stop(
      "The required reference distribution is not available for: ",
      test_name
    )
  }
  
  # Percentile rank in the appropriate reference distribution.
  percentile_I <- mean(
    I_reference <= observed_I,
    na.rm = TRUE
  ) * 100
  
  # Relative efficiency is available only when a constrained reference set exists.
  # For designs where L is not divisible by k, I_max_reference is NA and E(X)
  # is therefore not defined.
  if (is.finite(dist_obj$I_max_reference) &&
      dist_obj$I_max_reference > 0) {
    
    observed_E <- compute_E(
      observed_I,
      dist_obj
    )
    
  } else {
    
    observed_E <- NA_real_
    
  }
  
  # Exact feasible-design benchmark for small sequence spaces.
  # For L=3 and k=4, all 4^3 = 64 sequences are enumerated.
  feasible_ref <- get_feasible_reference(
    L = L,
    k = k,
    weights = weights,
    max_sequences = 1e6
  )
  
  observed_E_feasible <- compute_E_feasible(
    observed_I,
    feasible_ref
  )
  
  feasible_reference_type <- if (isTRUE(feasible_ref$exact)) {
    "Exact enumeration of all k^L sequences"
  } else {
    "Not calculated: k^L too large"
  }
  
  
  # Exact empirical upper-tail probability in the selected reference set.
  p_ge_I <- mean(
    I_reference >= observed_I,
    na.rm = TRUE
  )
  
  # Upper-tail reference benchmarks.
  I_p95 <- as.numeric(quantile(I_reference, 0.95, type = 7))
  I_p99 <- as.numeric(quantile(I_reference, 0.99, type = 7))
  
  # Text values for Word-ready Table 1.
  position_count_text <- paste(position_counts, collapse = ",")
  balanced_text <- ifelse(is_constrained_obs, "Yes", "No")
  
  I_percentile_text <- sprintf(
    "%.4f [%.1f%% percentile]",
    observed_I,
    percentile_I
  )
  
  list(
    test_name = test_name,
    reference = reference,
    L = L,
    k = k,
    sequence = X,
    position_counts = position_counts,
    is_constrained = is_constrained_obs,
    reference_type = reference_type,
    H1 = observed_H1,
    H2 = observed_H2,
    H3 = observed_H3,
    I = observed_I,
    I_per_trial = observed_I_per_trial,
    E = observed_E,
    E_feasible = observed_E_feasible,
    I_max_feasible = feasible_ref$I_max_feasible,
    n_feasible_sequences = feasible_ref$n_sequences,
    feasible_reference_exact = feasible_ref$exact,
    feasible_reference_type = feasible_reference_type,
    percentile_I = percentile_I,
    p_ge_I = p_ge_I,
    I_p95 = I_p95,
    I_p99 = I_p99,
    table_row = data.frame(
      `Odor identification test` = test_name,
      Reference = reference,
      L = L,
      k = k,
      `Position counts` = position_count_text,
      `Equal-frequency criterion fulfilled` = balanced_text,
      `H1 (bits)` = observed_H1,
      `H2 (bits)` = observed_H2,
      `H3 (bits)` = observed_H3,
      `I(X) (bits)` = observed_I,
      `I(X) percentile` = percentile_I,
      `I(X) [percentile]` = I_percentile_text,
      `I(X)/L (bits per trial)` = observed_I_per_trial,
      `Relative efficiency E(X)` = observed_E,
      `Feasible-design efficiency E_feasible(X)` = observed_E_feasible,
      `Feasible reference` = feasible_reference_type,
      `Reference distribution` = reference_type,
      `Empirical upper-tail P[I >= observed]` = p_ge_I,
      `I(X) 95th percentile` = I_p95,
      `I(X) 99th percentile` = I_p99,
      check.names = FALSE
    )
  )
}


# =============================================================================
# 17. Analyse all supplied sequences using lapply()
# =============================================================================

sequence_results <- lapply(
  sequence_list,
  analyse_sequence,
  distribution_objects = distribution_objects,
  weights = WEIGHTS
)

names(sequence_results) <- names(sequence_list)


# =============================================================================
# 18. Print individual sequence results
# =============================================================================

for (result_name in names(sequence_results)) {
  result <- sequence_results[[result_name]]
  
  cat("\n=============================================================\n")
  cat("Sequence:", result$test_name, "\n")
  cat("=============================================================\n")
  
  cat("L =", result$L, "; k =", result$k, "\n")
  cat("Position counts:", paste(result$position_counts, collapse = ", "), "\n")
  cat("Equal-frequency criterion:", ifelse(result$is_constrained, "Yes", "No"), "\n")
  cat("Reference distribution:", result$reference_type, "\n")
  
  cat(sprintf("H1(X)       = %.4f bits\n", result$H1))
  cat(sprintf("H2(X)       = %.4f bits\n", result$H2))
  cat(sprintf("H3(X)       = %.4f bits\n", result$H3))
  cat(sprintf(
    "I(X)        = %.4f bits [%.1f%% percentile]\n",
    result$I,
    result$percentile_I
  ))
  cat(sprintf(
    "I(X)/L      = %.4f bits per trial\n",
    result$I_per_trial
  ))
  if (is.finite(result$E)) {
    cat(sprintf(
      "E(X)        = %.4f\n",
      result$E
    ))
  } else {
    cat("E(X)        = NA (no constrained reference distribution: L is not divisible by k)\n")
  }
  if (is.finite(result$E_feasible)) {
    cat(sprintf(
      "E_feasible  = %.4f (exact maximum across %d sequences)\n",
      result$E_feasible,
      result$n_feasible_sequences
    ))
  } else {
    cat("E_feasible  = NA (sequence space not exhaustively enumerated)\n")
  }
  
  cat(sprintf(
    "P(I >= observed) = %.6f\n",
    result$p_ge_I
  ))
}


# =============================================================================
# 19. Assemble Table 1 for publication
# =============================================================================

table_1_publication <- do.call(
  rbind,
  lapply(sequence_results, function(x) x$table_row)
)

row.names(table_1_publication) <- NULL

# Round numeric columns for a clean export.
numeric_columns <- vapply(table_1_publication, is.numeric, logical(1))

table_1_publication[numeric_columns] <- lapply(
  table_1_publication[numeric_columns],
  round,
  digits = 4
)

cat("\n=============================================================\n")
cat("Table 1: publication-ready results\n")
cat("=============================================================\n")

print(table_1_publication, row.names = FALSE)

write.csv(
  table_1_publication,
  file = "Table_1_odor_identification_tests.csv",
  row.names = FALSE,
  na = ""
)


# =============================================================================
# 20. Write a compact Word-oriented table
# =============================================================================

# This compact version corresponds most closely to your current Table 1 layout.

table_1_compact <- table_1_publication[, c(
  "Odor identification test",
  "Reference",
  "L",
  "k",
  "Position counts",
  "Equal-frequency criterion fulfilled",
  "I(X) [percentile]",
  "I(X)/L (bits per trial)",
  "Relative efficiency E(X)",
  "Feasible-design efficiency E_feasible(X)"
)]

print(table_1_compact, row.names = FALSE)

write.csv(
  table_1_compact,
  file = "Table_1_odor_identification_tests_compact.csv",
  row.names = FALSE,
  na = ""
)

# =============================================================================
# 21. Maximum information across designs (L × k overview plot)
# =============================================================================

# Dense grids: constrained points (L divisible by k) plus non-divisible endpoints
overview_grid <- rbind(
  data.frame(L = c(3L, seq(4L, 40L, by = 4L)),          k = 4L),
  data.frame(L = c(3L, 4L, seq(6L, 36L, by = 6L), 40L), k = 6L)
)

# Reuse already-computed values where available; fall back to unconstrained max
get_cached_I_max <- function(L_val, k_val) {
  for (d in distribution_objects) {
    if (d$L == L_val && d$k == k_val) {
      if (is.finite(d$I_max_reference)) return(d$I_max_reference)
      if (!is.null(d$I_unconstrained))  return(max(d$I_unconstrained, na.rm = TRUE))
    }
  }
  NA_real_
}

overview_grid$I_max <- mapply(get_cached_I_max, overview_grid$L, overview_grid$k)

missing_idx <- which(!is.finite(overview_grid$I_max))

# Helper: constrained max when L%%k==0, otherwise unconstrained max
compute_overview_I_max <- function(L_val, k_val, n_sim = 20000L, seed = 123L) {
  if (L_val %% k_val == 0L) {
    compute_I_max(
      L = L_val, k = k_val, n_sim = n_sim, seed = seed,
      mc.cores = 1L, weights = WEIGHTS
    )$I_max
  } else {
    seeds <- make_replicate_seeds(n_sim, seed)
    I_vals <- vapply(seq_len(n_sim), function(i) {
      informativeness_score_general(
        generate_unconstrained_sequence(L_val, k_val, seeds[i]),
        k = k_val, weights = WEIGHTS
      )$score
    }, numeric(1L))
    max(I_vals, na.rm = TRUE)
  }
}

if (length(missing_idx) > 0L) {
  cat(sprintf(
    "Computing I_max for %d new (L, k) pairs...\n", length(missing_idx)
  ))
  new_vals <- pbmcapply::pbmclapply(
    missing_idx,
    function(i) compute_overview_I_max(overview_grid$L[i], overview_grid$k[i]),
    mc.cores    = max(1L, parallel::detectCores() - 1L),
    mc.set.seed = FALSE
  )
  overview_grid$I_max[missing_idx] <- unlist(new_vals)
}

overview_grid$k_label <- factor(paste0("k = ", overview_grid$k))

design_overview_plot <- ggplot(
    overview_grid,
    aes(x = L, y = I_max, colour = k_label, group = k_label)
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0) +
  scale_x_continuous(breaks = sort(unique(overview_grid$L))) +
  labs(
    title   = "Maximum Constrained Information Score by Design",
    subtitle = "Maximum over 20,000 Monte Carlo sequences; constrained (L divisible by k) or unconstrained otherwise",
    x       = "Number of test items (L)",
    y       = "Maximum I(X) [bit]",
    colour  = "Alternatives"
  ) +
  theme_plot() +
  theme(legend.position.inside = TRUE, legend.position = c(.8,.2)) +
  scale_color_manual(values = ggthemes::colorblind_pal()(8)[2:3])

print(design_overview_plot)

ggsave(
  filename = "design_overview_I_max.svg",
  plot     = design_overview_plot,
  width    = 9,
  height   = 5,
  dpi      = 180
)


# =============================================================================
# End of file
# =============================================================================