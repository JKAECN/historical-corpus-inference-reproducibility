# Historical Corpus Inference Reproducibility

Reproducibility materials for the manuscript:

**Balanced for what? Estimand-specific evidence structure in corpus design**

The study uses **Old Bailey Corpus 2.0 (OBC 2.0)** to examine how joint cell support, token depth, analysis-unit breadth, and higher-level source breadth affect the precision available for explicitly defined corpus-linguistic estimands.

## Important: download and extract the ZIP archives first

Several materials in this repository are distributed as `.zip` archives. They should be **downloaded and extracted locally before use**. The R scripts expect ordinary folders and files; they do not read directly from the ZIP archives.

The repository does **not** redistribute the raw OBC 2.0 corpus. OBC 2.0 should be obtained separately from its official source when reconstruction from the original corpus files is required.

## Which script should be run?

There are three stages in the analysis history:

1. **Original/base pipeline**
2. **Robustness extension (v3)**: `OBC2_Inferential_Alignment_IJCL_ROBUSTNESS_v3.R`
3. **Major-revision analysis (v4)**: `OBC2_Inferential_Alignment_IJCL_MAJOR_REVISION_v4.R`

The important distinction is that **v4 does not take the v3 output as its input**. Both v3 and v4 depend on files produced by the original/base pipeline.

Therefore:

- To reproduce the **main numerical results reported in the current manuscript**, extract the base archive and run **v4**.
- To reproduce the earlier **robustness analyses and audit trail**, run **v3** after the base pipeline.
- If reproducing the full development sequence, the chronological order is **base pipeline → v3 robustness → v4 major revision**, but **v3 is not a technical prerequisite for v4**.

The current manuscript's principal simulation results, analytic benchmarks, fixed-budget analysis, source-breadth analysis, smoothing sensitivity checks, and Monte Carlo standard errors are based on **v4**.

## Repository contents

### `OBC_IJCL_InferentialAlignment.zip`

Base analysis/output archive. After extraction, the folder should be named:

```text
OBC_IJCL_InferentialAlignment/
```

The most important downstream input is:

```text
03_speaker_trial_analysis_units.csv
```

The v3 robustness script also expects:

```text
03_speaker_trial_analysis_units.csv
06_simulation_truths.csv
08_primary_simulation_all_replications.csv
10_hierarchical_simulation_all_replications.csv
cache/OBC2_utterance_features.rds
```

### `OBC2_Inferential_Alignment_IJCL_ROBUSTNESS_v3.R`

Targeted robustness extension. It is designed to be run **after the original/base pipeline**.

It adds:

- estimand-specific sensitivity analyses;
- Proceeding-clustered robustness checks for the earlier hierarchical-source experiment;
- an alternative-role/full-time-span sensitivity analysis.

The script uses:

```r
20260924
```

and requires:

```r
data.table
ggplot2
clubSandwich
```

Its outputs are written to:

```text
OBC_IJCL_InferentialAlignment_ROBUSTNESS/
```

### `OBC_IJCL_InferentialAlignment_ROBUSTNESS.zip`

Archived output from the v3 robustness analysis.

These materials are retained for transparency and to document the analysis development. Some earlier robustness checks, including the CR2 comparison, are not treated as core evidence in the final manuscript.

### `OBC2_Inferential_Alignment_IJCL_MAJOR_REVISION_v4.R`

This is the **main analysis script for the current manuscript**.

It reads:

```text
OBC_IJCL_InferentialAlignment/03_speaker_trial_analysis_units.csv
```

and implements:

1. estimands defined outside model coding;
2. Gamma-Poisson empirical-Bayes rate calibration;
3. analytic variance benchmarks for log-rate contrasts;
4. fixed-total-token breadth-depth allocation;
5. source-matched Proceeding breadth under a fixed source-weighted target;
6. Monte Carlo standard errors;
7. raw-rate, add-0.5, and empirical-Bayes sensitivity analyses.

The script uses:

```r
set.seed(20260925)
```

and requires:

```r
data.table
ggplot2
```

Its main output files include:

```text
01_EB_parameters.csv
02_rate_diagnostics.csv
03_truth_sensitivity.csv
04_primary_EB_5000.csv
05_smoothing_sensitivity.csv
06_fixed_budget_allocation.csv
07_source_matched_breadth.csv
Figure_1_Benchmark.png
Figure_2_Allocation.png
Figure_3_SourceBreadth.png
99_sessionInfo.txt
```

The main and fixed-budget simulation conditions use **5,000 replications**; the source-breadth analysis uses **10,000 replications**; and the rate-calibration sensitivity conditions use **2,000 replications**.

### `OBC_IJCL_MAJOR_REVISION_v4.zip`

Archived output produced by the v4 major-revision analysis. This is the output archive most directly associated with the numerical results reported in the current manuscript.

### `IJCL_reproducibility_patch.zip`

Additional reproducibility documentation and reporting materials prepared during the final methodological audit. This archive is supplementary to the v4 results and is not needed to reproduce the numerical results already reported from v4 unless the added diagnostic/reporting files are specifically required.

## Recommended reproduction workflow

### A. Reproduce the current manuscript results

1. Download `OBC_IJCL_InferentialAlignment.zip`.
2. Extract it locally.
3. Keep the extracted folder name as `OBC_IJCL_InferentialAlignment`.
4. Place that folder in the location expected by the script (by default, the user's `Downloads` directory), or select the derived CSV manually if prompted.
5. Run `OBC2_Inferential_Alignment_IJCL_MAJOR_REVISION_v4.R`.
6. The script creates `OBC_IJCL_MAJOR_REVISION_v4/`.
7. Compare the generated files with `OBC_IJCL_MAJOR_REVISION_v4.zip`.

### B. Reproduce the robustness extension

1. Extract `OBC_IJCL_InferentialAlignment.zip`.
2. Confirm that the required original-pipeline files listed above are present.
3. Run `OBC2_Inferential_Alignment_IJCL_ROBUSTNESS_v3.R`.
4. The script writes results to `OBC_IJCL_InferentialAlignment_ROBUSTNESS/`.
5. Compare the generated files with `OBC_IJCL_InferentialAlignment_ROBUSTNESS.zip`.

## Relationship between v3 and v4

The version numbers reflect the **development history of the manuscript**, not a chained software pipeline.

```text
Original/base pipeline
        |
        +----> v3 robustness extension
        |
        +----> v4 major-revision analysis
```

Thus, **v3 does not need to finish before v4 can run**. What both scripts require is the derived output from the original/base pipeline. For the final manuscript, **v4 is the primary analysis**, while v3 is retained as a robustness and provenance record.

## Fixed-budget allocation note

In the fixed-budget experiment, integer cell allocation means that nominal and realised imbalance are not identical for every value of `N`. The current manuscript therefore distinguishes descriptive rounded conditions from exact matched comparisons and does not attribute all changes in the rounded extreme-imbalance curve solely to greater breadth.

## Software environment

The manuscript analyses were run in:

```text
R 4.6.1
Windows 11
```

The exact package and platform information for the v4 analysis is recorded in:

```text
99_sessionInfo.txt
```

## Data and code availability

OBC 2.0 is a publicly available historical corpus. This repository provides the analysis code, derived analysis materials, simulation outputs, figures, and reproducibility documentation used for the manuscript.

Repository:

https://github.com/JKAECN/historical-corpus-inference-reproducibility

## Citation

If this repository is used before the article has a formal bibliographic record, please cite the manuscript and this repository together. Full publication details can be added here after acceptance.
