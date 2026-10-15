# =============================================================================
# OBC 2.0 — Inferential Alignment: IJCL Robustness Extension
# Version: v3 targeted robustness
#
# PURPOSE
# -------
# This script DOES NOT rerun the full OBC parsing pipeline.
# It reads the outputs/caches already produced by the original analysis and
# adds three analyses requested for a stronger IJCL submission:
#
# A. ESTIMAND-SPECIFIC SENSITIVITY
#    Same corpus design, different target estimands:
#      - Role main effect
#      - Period main effect
#      - Role x Period interaction
#    Uses existing simulation replications. NO resimulation required.
#
# B. CLUSTER-AWARE ROBUSTNESS FOR EXPERIMENT 2
#    Repeats only the hierarchical-source-breadth experiment and compares:
#      - conventional quasi-Poisson inference
#      - CR2 cluster-robust inference clustered by Proceeding
#    This tests whether the low-source-breadth result is merely an artefact of
#    ignoring Proceeding-level dependence.
#
# C. ALTERNATIVE-ROLE / FULL-TIME-SPAN SENSITIVITY
#    Repeats the fixed-margin joint-support experiment with:
#      - witness vs defendant
#      - full OBC span 1720–1913
#      - two equal calendar-length periods:
#          1720–1816 (97 years)
#          1817–1913 (97 years)
#    This tests whether the main result depends on sparse late-period judge data.
#
# IMPORTANT
# ---------
# Run this AFTER the original pipeline. It expects the folder
# OBC_IJCL_InferentialAlignment to contain the original CSV files and cache.
# =============================================================================


# =============================================================================
# 0. USER SETTINGS
# =============================================================================

USER_HOME <- Sys.getenv("USERPROFILE")
if (!nzchar(USER_HOME)) USER_HOME <- path.expand("~")

DOWNLOADS_DIR <- file.path(USER_HOME, "Downloads")

# Original output folder from the first analysis.
RES_DIR <- file.path(DOWNLOADS_DIR, "OBC_IJCL_InferentialAlignment")

# New results go into a separate folder so nothing in the original analysis
# is overwritten.
ROBUST_DIR <- file.path(
  DOWNLOADS_DIR,
  "OBC_IJCL_InferentialAlignment_ROBUSTNESS"
)

# Main settings
RNG_SEED <- 20260924L

FEATURES <- c("modal_aux_n", "first_person_n")

# Existing design
TOKEN_BUDGET_FOR_ESTIMAND_FIGURE <- 250000L

# Cluster-aware Experiment 2
HIER_TOKEN_BUDGET <- 250000L
N_UNITS <- 200L
PROCEEDINGS_PER_CELL <- c(3L, 5L, 10L, 20L, 40L)
N_REPS_CLUSTER <- 500L

# Alternative-role sensitivity
ALT_ROLES <- c("witness", "defendant")
ALT_EARLY <- 1720:1816
ALT_LATE  <- 1817:1913
ALT_TOKEN_BUDGET <- 250000L
ALIGNMENT_DELTAS <- c(0.00, 0.25, 0.50, 0.75, 0.90)
N_REPS_ALT <- 500L

# If you only want a quick diagnostic run first, change both to 100.
# For manuscript results, 500 is recommended.


# =============================================================================
# 1. PACKAGES
# =============================================================================

needed <- c(
  "data.table",
  "ggplot2",
  "clubSandwich"
)

to_install <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
}

library(data.table)
library(ggplot2)
library(clubSandwich)

dir.create(ROBUST_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. LOCATE ORIGINAL RESULT FOLDER
# =============================================================================

required_files <- c(
  "03_speaker_trial_analysis_units.csv",
  "06_simulation_truths.csv",
  "08_primary_simulation_all_replications.csv",
  "10_hierarchical_simulation_all_replications.csv",
  file.path("cache", "OBC2_utterance_features.rds")
)

check_res_dir <- function(path) {
  dir.exists(path) &&
    all(file.exists(file.path(path, required_files)))
}

if (!check_res_dir(RES_DIR)) {

  if (.Platform$OS.type == "windows") {
    message(
      "Original result folder was not found automatically.\n",
      "Please select OBC_IJCL_InferentialAlignment."
    )

    chosen <- utils::choose.dir(
      default = DOWNLOADS_DIR,
      caption = "Select OBC_IJCL_InferentialAlignment"
    )

    if (!is.na(chosen) && nzchar(chosen)) {
      RES_DIR <- chosen
    }
  }
}

if (!check_res_dir(RES_DIR)) {
  stop(
    "\nCould not find all original analysis files under:\n  ",
    RES_DIR,
    "\n\nRequired files:\n  ",
    paste(required_files, collapse = "\n  "),
    "\n"
  )
}

message("Original analysis folder: ", RES_DIR)
message("New robustness folder: ", ROBUST_DIR)


# =============================================================================
# 3. GENERIC HELPERS
# =============================================================================

cell_counts_for_delta <- function(n_units, delta) {

  if (n_units %% 4L != 0L) {
    stop("n_units must be divisible by 4.")
  }

  if (delta < 0 || delta >= 1) {
    stop("delta must satisfy 0 <= delta < 1.")
  }

  half <- n_units / 2L
  n_diag <- as.integer(round((1 + delta) * n_units / 4))
  n_off  <- as.integer(half - n_diag)

  if (n_off < 1L) {
    stop("delta is too high for the chosen n_units.")
  }

  data.table(
    role = c("role0", "role0", "role1", "role1"),
    period = c("early", "late", "early", "late"),
    n = c(n_diag, n_off, n_off, n_diag)
  )
}


rate_truth_2x2 <- function(
    dt,
    role_col = "role",
    period_col = "period",
    rate_col,
    role0,
    role1
) {

  z <- dt[
    ,
    .(rate = mean(get(rate_col), na.rm = TRUE)),
    by = c(role_col, period_col)
  ]

  setnames(z, c(role_col, period_col), c("role_tmp", "period_tmp"))

  key <- paste(z$role_tmp, z$period_tmp, sep = "__")
  rr <- setNames(z$rate, key)

  required <- c(
    paste(role0, "early", sep = "__"),
    paste(role0, "late", sep = "__"),
    paste(role1, "early", sep = "__"),
    paste(role1, "late", sep = "__")
  )

  if (!all(required %in% names(rr))) {
    stop(
      "Missing at least one role x period cell when calculating truth."
    )
  }

  if (any(rr[required] <= 0)) {
    stop("At least one empirical cell rate is non-positive.")
  }

  r00 <- rr[paste(role0, "early", sep = "__")]
  r01 <- rr[paste(role0, "late",  sep = "__")]
  r10 <- rr[paste(role1, "early", sep = "__")]
  r11 <- rr[paste(role1, "late",  sep = "__")]

  data.table(
    term = c(
      "role_main",
      "period_main",
      "role_period_interaction"
    ),
    truth = as.numeric(c(
      log(r10) - log(r00),
      log(r01) - log(r00),
      log(r11) - log(r10) - log(r01) + log(r00)
    ))
  )
}


safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}


# =============================================================================
# PART A. ESTIMAND-SPECIFIC SENSITIVITY
# =============================================================================
#
# This part uses EXISTING replications, so it is fast.
# It directly tests the theoretical claim:
#   Alignment is a relation between corpus design and a target estimand,
#   not a single global property of the corpus.
# =============================================================================

message("\n[A] Estimand-specific sensitivity ...")

primary_rep <- fread(
  file.path(
    RES_DIR,
    "08_primary_simulation_all_replications.csv"
  )
)

term_labels <- c(
  "rolejudge" = "Role main effect",
  "period2late" = "Period main effect",
  "rolejudge:period2late" = "Role × Period interaction"
)

primary_rep[
  ,
  estimand := unname(term_labels[term])
]

if (any(is.na(primary_rep$estimand))) {
  stop("Unexpected term names in the primary simulation file.")
}

estimand_summary <- primary_rep[
  ,
  .(
    n_rep = .N,
    bias = safe_mean(estimate - truth),
    rmse = sqrt(safe_mean((estimate - truth)^2)),
    mean_se = safe_mean(se),
    mean_ci_width = safe_mean(ci_high - ci_low),
    coverage = safe_mean(covered),
    sign_recovery = safe_mean(sign_recovered)
  ),
  by = .(
    feature,
    token_budget,
    delta,
    term,
    estimand
  )
]

# RMSE relative to the perfectly crossed delta=0 design.
estimand_summary[
  ,
  rmse_at_delta0 := rmse[delta == 0][1],
  by = .(feature, token_budget, term)
]

estimand_summary[
  ,
  rmse_index := rmse / rmse_at_delta0
]

fwrite(
  estimand_summary,
  file.path(
    ROBUST_DIR,
    "14_estimand_specific_summary.csv"
  )
)

# Endpoint table: delta=.90 compared with delta=0.
estimand_endpoint <- estimand_summary[
  delta %in% c(0, 0.90),
  .(
    rmse_delta0 = rmse[delta == 0][1],
    rmse_delta90 = rmse[delta == 0.90][1],
    percent_increase = 100 * (
      rmse[delta == 0.90][1] /
        rmse[delta == 0][1] - 1
    ),
    ci_width_delta0 = mean_ci_width[delta == 0][1],
    ci_width_delta90 = mean_ci_width[delta == 0.90][1],
    coverage_delta0 = coverage[delta == 0][1],
    coverage_delta90 = coverage[delta == 0.90][1]
  ),
  by = .(
    feature,
    token_budget,
    term,
    estimand
  )
]

fwrite(
  estimand_endpoint,
  file.path(
    ROBUST_DIR,
    "15_estimand_specific_endpoint_change.csv"
  )
)

plot_est <- estimand_summary[
  token_budget == TOKEN_BUDGET_FOR_ESTIMAND_FIGURE
]

p_est <- ggplot(
  plot_est,
  aes(
    x = delta,
    y = rmse_index,
    linetype = estimand,
    shape = estimand,
    group = estimand
  )
) +
  geom_hline(
    yintercept = 1,
    linewidth = 0.45,
    linetype = "dotted"
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_wrap(~ feature) +
  labs(
    title = "The same corpus design has different costs for different estimands",
    subtitle = paste0(
      "RMSE indexed to the fully crossed design (delta = 0); ",
      TOKEN_BUDGET_FOR_ESTIMAND_FIGURE,
      " tokens, fixed N and fixed marginal Role/Period balance"
    ),
    x = "Joint-distribution misalignment (delta)",
    y = "RMSE relative to delta = 0",
    linetype = "Target estimand",
    shape = "Target estimand"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    ROBUST_DIR,
    "Figure_5_Estimand_Specific_RMSE_Index.png"
  ),
  p_est,
  width = 9.2,
  height = 5.5,
  dpi = 400
)


# =============================================================================
# PART B. CLUSTER-AWARE ROBUSTNESS FOR EXPERIMENT 2
# =============================================================================
#
# Question:
#   Is the source-breadth result merely a consequence of using conventional
#   model-based SEs while observations are clustered within Proceedings?
#
# Strategy:
#   Repeat Experiment 2 with the SAME data-generating design.
#   Fit the same quasi-Poisson model, but calculate:
#     1. conventional model-based CI
#     2. CR2 small-sample cluster-robust CI, clustered by Proceeding
#
# Important:
#   CR2 changes uncertainty estimates, not the quasi-Poisson point estimate.
#   Therefore:
#     - if RMSE still falls with broader source coverage, that pattern cannot
#       be explained away by SE correction;
#     - if CR2 repairs low-breadth CI coverage, the manuscript should say that
#       conventional intervals were anti-conservative under source clustering.
# =============================================================================

message("\n[B] Cluster-aware robustness for Experiment 2 ...")

units <- fread(
  file.path(
    RES_DIR,
    "03_speaker_trial_analysis_units.csv"
  )
)

units[
  ,
  role := tolower(role)
]

units[
  ,
  period2 := factor(
    period2,
    levels = c("early", "late")
  )
]

primary <- units[
  role %in% c("witness", "judge") &
    period2 %in% c("early", "late") &
    tokens > 0
]

primary[
  ,
  role := factor(
    role,
    levels = c("witness", "judge")
  )
]

for (feat in FEATURES) {
  rate_col <- paste0(feat, "_rate")
  primary[
    ,
    (rate_col) := (get(feat) + 0.5) / (tokens + 1)
  ]
}

truth_primary <- fread(
  file.path(
    RES_DIR,
    "06_simulation_truths.csv"
  )
)

truth_primary[
  ,
  term_generic := fcase(
    term == "rolejudge", "role_main",
    term == "period2late", "period_main",
    term == "rolejudge:period2late",
      "role_period_interaction",
    default = NA_character_
  )
]

draw_hier_cell <- function(
    source,
    role_value,
    period_value,
    n,
    feature,
    tokens_per_unit,
    proceedings_per_cell
) {

  pool <- source[
    role == role_value &
      period2 == period_value &
      !is.na(proceeding) &
      nzchar(proceeding)
  ]

  avail_sources <- unique(pool$proceeding)

  if (!length(avail_sources)) {
    stop(
      "No Proceedings in cell: ",
      role_value,
      " / ",
      period_value
    )
  }

  k <- min(
    proceedings_per_cell,
    length(avail_sources)
  )

  selected_sources <- sample(
    avail_sources,
    size = k,
    replace = FALSE
  )

  pool2 <- pool[
    proceeding %in% selected_sources
  ]

  idx <- sample.int(
    nrow(pool2),
    size = n,
    replace = TRUE
  )

  d <- pool2[idx]

  rate_col <- paste0(feature, "_rate")
  lambda <- pmax(
    d[[rate_col]],
    0
  ) * tokens_per_unit

  data.table(
    role = factor(
      rep(role_value, n),
      levels = c("witness", "judge")
    ),
    period2 = factor(
      rep(as.character(period_value), n),
      levels = c("early", "late")
    ),
    tokens = rep(tokens_per_unit, n),
    count = rpois(n, lambda = lambda),
    source_proceeding = d$proceeding,
    source_trial = d$trial_key,
    source_unit = d$speaker_trial_id
  )
}


simulate_hierarchical <- function(
    source,
    feature,
    token_budget,
    n_units,
    proceedings_per_cell
) {

  if (token_budget %% n_units != 0L) {
    stop("token_budget must be divisible by n_units.")
  }

  tokens_per_unit <- token_budget / n_units
  n_cell <- n_units / 4L

  specs <- data.table(
    role = c(
      "witness", "witness",
      "judge", "judge"
    ),
    period2 = c(
      "early", "late",
      "early", "late"
    )
  )

  pieces <- vector(
    "list",
    nrow(specs)
  )

  for (j in seq_len(nrow(specs))) {
    pieces[[j]] <- draw_hier_cell(
      source = source,
      role_value = specs$role[j],
      period_value = specs$period2[j],
      n = n_cell,
      feature = feature,
      tokens_per_unit = tokens_per_unit,
      proceedings_per_cell = proceedings_per_cell
    )
  }

  sim <- rbindlist(pieces)

  stopifnot(
    nrow(sim) == n_units,
    sum(sim$tokens) == token_budget,
    all(table(sim$role) == n_units / 2),
    all(table(sim$period2) == n_units / 2)
  )

  sim
}


extract_naive_and_cr2 <- function(sim) {

  fit <- tryCatch(
    glm(
      count ~ role * period2 +
        offset(log(tokens)),
      family = quasipoisson(link = "log"),
      data = sim
    ),
    error = function(e) NULL
  )

  wanted <- c(
    "rolejudge",
    "period2late",
    "rolejudge:period2late"
  )

  generic <- c(
    "role_main",
    "period_main",
    "role_period_interaction"
  )

  if (is.null(fit)) {
    return(
      rbindlist(
        lapply(
          c("naive", "CR2"),
          function(m) {
            data.table(
              model = m,
              term = wanted,
              term_generic = generic,
              estimate = NA_real_,
              se = NA_real_,
              df = NA_real_,
              ci_low = NA_real_,
              ci_high = NA_real_,
              converged = FALSE
            )
          }
        )
      )
    )
  }

  # ----- conventional model-based inference -----
  sm <- summary(fit)$coefficients

  naive <- data.table(
    model = "naive",
    term = wanted,
    term_generic = generic,
    estimate = NA_real_,
    se = NA_real_,
    df = Inf,
    ci_low = NA_real_,
    ci_high = NA_real_,
    converged = isTRUE(fit$converged)
  )

  present <- intersect(
    wanted,
    rownames(sm)
  )

  if (length(present)) {

    idx_out <- match(
      present,
      naive$term
    )

    est <- sm[
      present,
      "Estimate"
    ]

    se0 <- sm[
      present,
      "Std. Error"
    ]

    naive[
      idx_out,
      `:=`(
        estimate = est,
        se = se0,
        ci_low = est - 1.96 * se0,
        ci_high = est + 1.96 * se0
      )
    ]
  }

  # ----- CR2 Proceeding-clustered inference -----
  cr2 <- data.table(
    model = "CR2",
    term = wanted,
    term_generic = generic,
    estimate = NA_real_,
    se = NA_real_,
    df = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    converged = isTRUE(fit$converged)
  )

  ct <- tryCatch(
    clubSandwich::coef_test(
      fit,
      vcov = "CR2",
      cluster = sim$source_proceeding,
      test = "Satterthwaite"
    ),
    error = function(e) NULL
  )

  if (!is.null(ct)) {

    ct_df <- as.data.frame(ct)
    ct_df$term <- rownames(ct_df)

    for (tt in wanted) {

      rr <- ct_df[
        ct_df$term == tt,
        ,
        drop = FALSE
      ]

      if (!nrow(rr)) next

      est <- as.numeric(rr$beta[1])
      se1 <- as.numeric(rr$SE[1])

      # clubSandwich normally returns df_Satt for Satterthwaite.
      df1 <- if ("df_Satt" %in% names(rr)) {
        as.numeric(rr$df_Satt[1])
      } else if ("df" %in% names(rr)) {
        as.numeric(rr$df[1])
      } else {
        NA_real_
      }

      crit <- if (
        is.finite(df1) &&
          df1 > 0
      ) {
        qt(0.975, df = df1)
      } else {
        1.96
      }

      ii <- match(tt, cr2$term)

      cr2[
        ii,
        `:=`(
          estimate = est,
          se = se1,
          df = df1,
          ci_low = est - crit * se1,
          ci_high = est + crit * se1
        )
      ]
    }
  }

  rbind(
    naive,
    cr2,
    fill = TRUE
  )
}


set.seed(RNG_SEED + 100L)

cluster_grid <- CJ(
  feature = FEATURES,
  proceedings_per_cell = PROCEEDINGS_PER_CELL,
  rep = seq_len(N_REPS_CLUSTER),
  sorted = TRUE
)

message(
  "Cluster-aware fits to run: ",
  format(
    nrow(cluster_grid),
    big.mark = ","
  )
)

cluster_list <- vector(
  "list",
  nrow(cluster_grid)
)

for (i in seq_len(nrow(cluster_grid))) {

  if (
    i == 1L ||
      i %% 250L == 0L ||
      i == nrow(cluster_grid)
  ) {
    message(
      "  Cluster robustness ",
      i,
      " / ",
      nrow(cluster_grid)
    )
  }

  g <- cluster_grid[i]

  sim <- simulate_hierarchical(
    source = primary,
    feature = g$feature,
    token_budget = HIER_TOKEN_BUDGET,
    n_units = N_UNITS,
    proceedings_per_cell =
      g$proceedings_per_cell
  )

  z <- extract_naive_and_cr2(sim)

  z[
    ,
    `:=`(
      feature = g$feature,
      token_budget = HIER_TOKEN_BUDGET,
      proceedings_per_cell =
        g$proceedings_per_cell,
      rep = g$rep,
      realized_proceedings =
        uniqueN(sim$source_proceeding),
      realized_trials =
        uniqueN(sim$source_trial)
    )
  ]

  cluster_list[[i]] <- z
}

cluster_res <- rbindlist(
  cluster_list,
  fill = TRUE
)

rm(cluster_list)
gc()

cluster_res <- merge(
  cluster_res,
  truth_primary[
    ,
    .(
      feature,
      term_generic,
      truth
    )
  ],
  by = c(
    "feature",
    "term_generic"
  ),
  all.x = TRUE
)

cluster_res[
  ,
  error := estimate - truth
]

cluster_res[
  ,
  covered :=
    is.finite(ci_low) &
    is.finite(ci_high) &
    ci_low <= truth &
    ci_high >= truth
]

fwrite(
  cluster_res,
  file.path(
    ROBUST_DIR,
    "16_cluster_robust_all_replications.csv"
  )
)

cluster_summary <- cluster_res[
  ,
  .(
    n_rep = .N,
    valid_estimate_rate =
      mean(is.finite(estimate)),
    valid_ci_rate =
      mean(
        is.finite(ci_low) &
          is.finite(ci_high)
      ),
    bias =
      safe_mean(error),
    rmse =
      sqrt(safe_mean(error^2)),
    mean_se =
      safe_mean(se),
    mean_ci_width =
      safe_mean(ci_high - ci_low),
    coverage =
      safe_mean(covered),
    mean_realized_proceedings =
      safe_mean(realized_proceedings),
    mean_realized_trials =
      safe_mean(realized_trials),
    median_df =
      safe_median(df)
  ),
  by = .(
    model,
    feature,
    proceedings_per_cell,
    term_generic
  )
][
  order(
    feature,
    term_generic,
    proceedings_per_cell,
    model
  )
]

fwrite(
  cluster_summary,
  file.path(
    ROBUST_DIR,
    "17_cluster_robust_summary.csv"
  )
)

p_cluster_cov <- ggplot(
  cluster_summary[
    term_generic ==
      "role_period_interaction"
  ],
  aes(
    x = proceedings_per_cell,
    y = coverage,
    linetype = model,
    shape = model,
    group = model
  )
) +
  geom_hline(
    yintercept = 0.95,
    linewidth = 0.45,
    linetype = "dotted"
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_wrap(~ feature) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    title = "Does cluster-robust inference repair low-source coverage?",
    subtitle = "Same simulated corpora; conventional versus Proceeding-clustered CR2 confidence intervals",
    x = "Requested Proceedings per Role × Period cell",
    y = "Empirical 95% CI coverage",
    linetype = "Inference",
    shape = "Inference"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    ROBUST_DIR,
    "Figure_6_Cluster_Robust_Coverage.png"
  ),
  p_cluster_cov,
  width = 9.2,
  height = 5.5,
  dpi = 400
)

# RMSE does not depend on the SE estimator for naive vs CR2, but we save a
# separate source-breadth plot to make that point explicit.
p_cluster_rmse <- ggplot(
  unique(
    cluster_summary[
      model == "naive" &
        term_generic ==
          "role_period_interaction",
      .(
        feature,
        proceedings_per_cell,
        rmse
      )
    ]
  ),
  aes(
    x = proceedings_per_cell,
    y = rmse
  )
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_wrap(
    ~ feature,
    scales = "free_y"
  ) +
  labs(
    title = "Point-estimate instability persists independently of SE correction",
    subtitle = "RMSE of the Role × Period interaction across source-breadth conditions",
    x = "Requested Proceedings per Role × Period cell",
    y = "RMSE"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(
    ROBUST_DIR,
    "Figure_7_Source_Breadth_RMSE_Robustness.png"
  ),
  p_cluster_rmse,
  width = 8.8,
  height = 5.3,
  dpi = 400
)


# =============================================================================
# PART C. ALTERNATIVE ROLE CONTRAST ACROSS THE FULL 1720–1913 SPAN
# =============================================================================
#
# Why witness vs defendant?
#   Both roles have substantial support throughout all five OBC periods,
#   including 1880–1913. This avoids the sharp late-period sparsity of judges.
#
# Why split at 1816/1817?
#   1720–1913 contains 194 calendar years. The split creates two equal-length
#   97-year windows. It is used only as a robustness partition, not proposed
#   as a substantive historical periodization.
# =============================================================================

message(
  "\n[C] Alternative-role/full-span sensitivity ..."
)

utt_path <- file.path(
  RES_DIR,
  "cache",
  "OBC2_utterance_features.rds"
)

utt <- readRDS(utt_path)
setDT(utt)

# Reconstruct safe hierarchical IDs exactly as in the original pipeline.
utt[
  ,
  trial_key := fifelse(
    !is.na(trial) &
      nzchar(trial),
    paste(
      proceeding,
      trial,
      sep = "::"
    ),
    paste(
      proceeding,
      event,
      sep = "::"
    )
  )
]

utt[
  ,
  speaker_local := fifelse(
    !is.na(speaker) &
      nzchar(speaker),
    speaker,
    paste0(
      "anonymous_",
      utterance_n
    )
  )
]

utt[
  ,
  speaker_trial_id := paste(
    trial_key,
    speaker_local,
    sep = "::"
  )
]

utt[
  ,
  alt_period := fcase(
    year %in% ALT_EARLY, "early",
    year %in% ALT_LATE, "late",
    default = NA_character_
  )
]

alt_units <- utt[
  role %in% ALT_ROLES &
    !is.na(alt_period) &
    wc > 0,
  .(
    tokens = sum(
      wc,
      na.rm = TRUE
    ),
    utterances = .N,
    modal_aux_n = sum(
      modal_aux_n,
      na.rm = TRUE
    ),
    first_person_n = sum(
      first_person_n,
      na.rm = TRUE
    )
  ),
  by = .(
    proceeding,
    trial_key,
    speaker_trial_id,
    role,
    alt_period
  )
]

alt_units[
  ,
  role := factor(
    role,
    levels = ALT_ROLES
  )
]

alt_units[
  ,
  alt_period := factor(
    alt_period,
    levels = c(
      "early",
      "late"
    )
  )
]

alt_support <- alt_units[
  ,
  .(
    units = .N,
    tokens = sum(tokens),
    proceedings =
      uniqueN(proceeding),
    trials =
      uniqueN(trial_key)
  ),
  by = .(
    role,
    alt_period
  )
]

fwrite(
  alt_support,
  file.path(
    ROBUST_DIR,
    "18_alternative_role_fullspan_support.csv"
  )
)

print(alt_support)

if (
  nrow(alt_support) < 4L ||
    any(alt_support$units < 50L)
) {
  warning(
    "Alternative-role support is unexpectedly sparse. ",
    "Inspect 18_alternative_role_fullspan_support.csv."
  )
}

for (feat in FEATURES) {
  rate_col <- paste0(
    feat,
    "_rate"
  )

  alt_units[
    ,
    (rate_col) :=
      (get(feat) + 0.5) /
      (tokens + 1)
  ]
}

alt_truth <- rbindlist(
  lapply(
    FEATURES,
    function(feat) {

      rate_col <- paste0(
        feat,
        "_rate"
      )

      z <- copy(alt_units)
      z[
        ,
        period_generic :=
          as.character(alt_period)
      ]

      tr <- rate_truth_2x2(
        dt = z,
        role_col = "role",
        period_col = "period_generic",
        rate_col = rate_col,
        role0 = ALT_ROLES[1],
        role1 = ALT_ROLES[2]
      )

      tr[
        ,
        feature := feat
      ]

      tr
    }
  )
)
# 修正列名：alt_truth 中原来叫 term，
# 后续 merge 使用的是 term_generic
if ("term" %in% names(alt_truth) &&
    !"term_generic" %in% names(alt_truth)) {
  data.table::setnames(
    alt_truth,
    old = "term",
    new = "term_generic"
  )
}

# 检查
print(names(alt_truth))
print(alt_truth)
fwrite(
  alt_truth,
  file.path(
    ROBUST_DIR,
    "19_alternative_role_truths.csv"
  )
)


draw_alt_cell <- function(
    source,
    role_value,
    period_value,
    n,
    feature,
    tokens_per_unit
) {

  pool <- source[
    role == role_value &
      alt_period == period_value
  ]

  if (!nrow(pool)) {
    stop(
      "No source units in alternative cell: ",
      role_value,
      " / ",
      period_value
    )
  }

  idx <- sample.int(
    nrow(pool),
    size = n,
    replace = TRUE
  )

  d <- pool[idx]

  rate_col <- paste0(
    feature,
    "_rate"
  )

  lambda <- pmax(
    d[[rate_col]],
    0
  ) * tokens_per_unit

  data.table(
    role = factor(
      rep(role_value, n),
      levels = ALT_ROLES
    ),
    period = factor(
      rep(period_value, n),
      levels = c(
        "early",
        "late"
      )
    ),
    tokens =
      rep(tokens_per_unit, n),
    count =
      rpois(n, lambda = lambda),
    source_proceeding =
      d$proceeding,
    source_trial =
      d$trial_key
  )
}


simulate_alt <- function(
    source,
    feature,
    token_budget,
    n_units,
    delta
) {

  if (
    token_budget %% n_units != 0L
  ) {
    stop(
      "token_budget must be divisible by n_units."
    )
  }

  tokens_per_unit <-
    token_budget / n_units

  cc <- cell_counts_for_delta(
    n_units,
    delta
  )

  cc[
    role == "role0",
    role_value := ALT_ROLES[1]
  ]

  cc[
    role == "role1",
    role_value := ALT_ROLES[2]
  ]

  pieces <- vector(
    "list",
    nrow(cc)
  )

  for (j in seq_len(nrow(cc))) {

    pieces[[j]] <- draw_alt_cell(
      source = source,
      role_value =
        cc$role_value[j],
      period_value =
        cc$period[j],
      n = cc$n[j],
      feature = feature,
      tokens_per_unit =
        tokens_per_unit
    )
  }

  sim <- rbindlist(pieces)

  stopifnot(
    nrow(sim) == n_units,
    sum(sim$tokens) ==
      token_budget,
    all(
      table(sim$role) ==
        n_units / 2
    ),
    all(
      table(sim$period) ==
        n_units / 2
    )
  )

  sim
}


fit_alt <- function(sim) {

  fit <- tryCatch(
    glm(
      count ~ role * period +
        offset(log(tokens)),
      family =
        quasipoisson(link = "log"),
      data = sim
    ),
    error = function(e) NULL
  )

  # With factor levels witness, defendant and early, late:
  wanted <- c(
    paste0(
      "role",
      ALT_ROLES[2]
    ),
    "periodlate",
    paste0(
      "role",
      ALT_ROLES[2],
      ":periodlate"
    )
  )

  generic <- c(
    "role_main",
    "period_main",
    "role_period_interaction"
  )

  out <- data.table(
    term = wanted,
    term_generic = generic,
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    converged = FALSE
  )

  if (is.null(fit)) {
    return(out)
  }

  sm <- summary(fit)$coefficients

  present <- intersect(
    wanted,
    rownames(sm)
  )

  if (length(present)) {

    ii <- match(
      present,
      out$term
    )

    est <- sm[
      present,
      "Estimate"
    ]

    se0 <- sm[
      present,
      "Std. Error"
    ]

    out[
      ii,
      `:=`(
        estimate = est,
        se = se0,
        ci_low = est -
          1.96 * se0,
        ci_high = est +
          1.96 * se0
      )
    ]
  }

  out[
    ,
    converged :=
      isTRUE(fit$converged)
  ]

  out
}


set.seed(RNG_SEED + 200L)

alt_grid <- CJ(
  feature = FEATURES,
  delta = ALIGNMENT_DELTAS,
  rep = seq_len(N_REPS_ALT),
  sorted = TRUE
)

message(
  "Alternative-role fits to run: ",
  format(
    nrow(alt_grid),
    big.mark = ","
  )
)

alt_list <- vector(
  "list",
  nrow(alt_grid)
)

for (i in seq_len(nrow(alt_grid))) {

  if (
    i == 1L ||
      i %% 500L == 0L ||
      i == nrow(alt_grid)
  ) {
    message(
      "  Alternative-role simulation ",
      i,
      " / ",
      nrow(alt_grid)
    )
  }

  g <- alt_grid[i]

  sim <- simulate_alt(
    source = alt_units,
    feature = g$feature,
    token_budget =
      ALT_TOKEN_BUDGET,
    n_units = N_UNITS,
    delta = g$delta
  )

  z <- fit_alt(sim)

  z[
    ,
    `:=`(
      feature = g$feature,
      token_budget =
        ALT_TOKEN_BUDGET,
      delta = g$delta,
      rep = g$rep,
      n_unique_proceedings =
        uniqueN(
          sim$source_proceeding
        ),
      n_unique_trials =
        uniqueN(
          sim$source_trial
        )
    )
  ]

  alt_list[[i]] <- z
}

alt_res <- rbindlist(
  alt_list,
  fill = TRUE
)

rm(alt_list)
gc()

alt_res <- merge(
  alt_res,
  alt_truth[
    ,
    .(
      feature,
      term_generic,
      truth
    )
  ],
  by = c(
    "feature",
    "term_generic"
  ),
  all.x = TRUE
)

alt_res[
  ,
  `:=`(
    error = estimate - truth,
    covered =
      is.finite(ci_low) &
      is.finite(ci_high) &
      ci_low <= truth &
      ci_high >= truth
  )
]

fwrite(
  alt_res,
  file.path(
    ROBUST_DIR,
    "20_alternative_role_all_replications.csv"
  )
)

alt_summary <- alt_res[
  ,
  .(
    n_rep = .N,
    estimable_rate =
      mean(is.finite(estimate)),
    convergence_rate =
      mean(converged, na.rm = TRUE),
    bias =
      safe_mean(error),
    rmse =
      sqrt(safe_mean(error^2)),
    mean_se =
      safe_mean(se),
    mean_ci_width =
      safe_mean(ci_high - ci_low),
    coverage =
      safe_mean(covered),
    mean_unique_proceedings =
      safe_mean(n_unique_proceedings),
    mean_unique_trials =
      safe_mean(n_unique_trials)
  ),
  by = .(
    feature,
    token_budget,
    delta,
    term_generic
  )
]

alt_summary[
  ,
  rmse_at_delta0 :=
    rmse[delta == 0][1],
  by = .(
    feature,
    term_generic
  )
]

alt_summary[
  ,
  rmse_index :=
    rmse / rmse_at_delta0
]

fwrite(
  alt_summary,
  file.path(
    ROBUST_DIR,
    "21_alternative_role_summary.csv"
  )
)

p_alt <- ggplot(
  alt_summary[
    term_generic ==
      "role_period_interaction"
  ],
  aes(
    x = delta,
    y = rmse_index
  )
) +
  geom_hline(
    yintercept = 1,
    linewidth = 0.45,
    linetype = "dotted"
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.4) +
  facet_wrap(~ feature) +
  labs(
    title = "Alternative-role sensitivity across the full OBC time span",
    subtitle = paste0(
      "Witness vs defendant; 1720–1816 versus 1817–1913; ",
      "fixed ",
      ALT_TOKEN_BUDGET,
      " tokens, N = ",
      N_UNITS,
      " and fixed marginal balance"
    ),
    x = "Joint-distribution misalignment (delta)",
    y = "Interaction RMSE relative to delta = 0"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(
    ROBUST_DIR,
    "Figure_8_Alternative_Role_FullSpan_RMSE.png"
  ),
  p_alt,
  width = 8.8,
  height = 5.3,
  dpi = 400
)


# =============================================================================
# 4. AUTO-GENERATED MANUSCRIPT SUMMARY
# =============================================================================

message(
  "\nGenerating manuscript-ready numerical summary ..."
)

summary_file <- file.path(
  ROBUST_DIR,
  "22_MANUSCRIPT_NUMERICAL_SUMMARY.txt"
)

con <- file(
  summary_file,
  open = "wt",
  encoding = "UTF-8"
)

writeLines(
  c(
    "OBC INFERENTIAL ALIGNMENT — ROBUSTNESS SUMMARY",
    "================================================",
    "",
    "A. ESTIMAND-SPECIFIC SENSITIVITY",
    paste0(
      "Token budget highlighted in Figure 5: ",
      TOKEN_BUDGET_FOR_ESTIMAND_FIGURE
    ),
    ""
  ),
  con
)

for (feat in FEATURES) {

  writeLines(
    paste0("Feature: ", feat),
    con
  )

  zz <- estimand_endpoint[
    feature == feat &
      token_budget ==
        TOKEN_BUDGET_FOR_ESTIMAND_FIGURE
  ]

  for (j in seq_len(nrow(zz))) {
    writeLines(
      sprintf(
        "  %s: RMSE %.4f -> %.4f (%+.1f%%); CI width %.4f -> %.4f; coverage %.3f -> %.3f",
        zz$estimand[j],
        zz$rmse_delta0[j],
        zz$rmse_delta90[j],
        zz$percent_increase[j],
        zz$ci_width_delta0[j],
        zz$ci_width_delta90[j],
        zz$coverage_delta0[j],
        zz$coverage_delta90[j]
      ),
      con
    )
  }

  writeLines("", con)
}

writeLines(
  c(
    "",
    "B. CLUSTER-AWARE EXPERIMENT 2",
    ""
  ),
  con
)

for (feat in FEATURES) {

  writeLines(
    paste0("Feature: ", feat),
    con
  )

  zz <- cluster_summary[
    feature == feat &
      term_generic ==
        "role_period_interaction" &
      proceedings_per_cell %in%
        c(
          min(PROCEEDINGS_PER_CELL),
          max(PROCEEDINGS_PER_CELL)
        )
  ]

  for (k in c(
    min(PROCEEDINGS_PER_CELL),
    max(PROCEEDINGS_PER_CELL)
  )) {

    zzk <- zz[
      proceedings_per_cell == k
    ]

    writeLines(
      paste0(
        "  Proceedings/cell = ",
        k
      ),
      con
    )

    for (j in seq_len(nrow(zzk))) {

      writeLines(
        sprintf(
          "    %s: RMSE %.4f; CI width %.4f; coverage %.3f; valid CI %.3f",
          zzk$model[j],
          zzk$rmse[j],
          zzk$mean_ci_width[j],
          zzk$coverage[j],
          zzk$valid_ci_rate[j]
        ),
        con
      )
    }
  }

  writeLines("", con)
}

writeLines(
  c(
    "",
    "C. ALTERNATIVE ROLE / FULL-SPAN SENSITIVITY",
    paste0(
      "Roles: ",
      ALT_ROLES[1],
      " vs ",
      ALT_ROLES[2]
    ),
    "Periods: 1720–1816 vs 1817–1913",
    ""
  ),
  con
)

for (feat in FEATURES) {

  zz <- alt_summary[
    feature == feat &
      term_generic ==
        "role_period_interaction" &
      delta %in% c(
        0,
        0.90
      )
  ]

  if (nrow(zz) == 2L) {

    z0 <- zz[
      delta == 0
    ]

    z9 <- zz[
      delta == 0.90
    ]

    writeLines(
      sprintf(
        "  %s: interaction RMSE %.4f -> %.4f (index %.2f); CI width %.4f -> %.4f; coverage %.3f -> %.3f",
        feat,
        z0$rmse,
        z9$rmse,
        z9$rmse_index,
        z0$mean_ci_width,
        z9$mean_ci_width,
        z0$coverage,
        z9$coverage
      ),
      con
    )
  }
}

close(con)


# =============================================================================
# 5. SESSION INFO + COMPLETION
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    ROBUST_DIR,
    "99_sessionInfo_robustness.txt"
  )
)

cat(
  "\n============================================================\n",
  "TARGETED IJCL ROBUSTNESS ANALYSES COMPLETED\n",
  "============================================================\n",
  "Output folder:\n",
  normalizePath(
    ROBUST_DIR,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\nKey files to send back to ChatGPT:\n",
  "  14_estimand_specific_summary.csv\n",
  "  15_estimand_specific_endpoint_change.csv\n",
  "  17_cluster_robust_summary.csv\n",
  "  18_alternative_role_fullspan_support.csv\n",
  "  21_alternative_role_summary.csv\n",
  "  22_MANUSCRIPT_NUMERICAL_SUMMARY.txt\n",
  "  Figure_5_Estimand_Specific_RMSE_Index.png\n",
  "  Figure_6_Cluster_Robust_Coverage.png\n",
  "  Figure_7_Source_Breadth_RMSE_Robustness.png\n",
  "  Figure_8_Alternative_Role_FullSpan_RMSE.png\n",
  "\nYou do NOT need to rerun the original XML parsing or original simulations.\n",
  sep = ""
)
