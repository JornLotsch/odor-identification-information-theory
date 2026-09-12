
# Information theoretic analysis of randomness and structure in clinical olfactory test design

This repository is a computational companion to the accompanying manuscript:

Lötsch J, Hummel T, Himmelspach A. Information theoretic analysis of randomness and structure in clinical olfactory test design.

It contains the code used to implement the theoretical reference calculations and the empirical analysis of real identification-choice data described in the paper. The repository is intentionally kept concise and does not duplicate the manuscript; it provides the reproducible computational workflow that supports the paper.

## Example figure

The figure below shows the constrained reference distribution of the informativeness score $I(X)$ for a 18-trial, 6-alternative design. It is the example output directly generated from the theoretical analysis workflow and is included here to illustrate the kind of distribution used in the paper.

![Constrained informativeness distribution for a 18-trial, 6-alternative design](./scores_L18_k6_I_constrained_density.svg)

## Repository scope

The active workflow for this project is limited to the following files:

- `R/globals.R` — shared package setup, helper functions, and common plotting utilities
- `R/theoretical_distributions_and_scores.R` — theoretical chance and informativeness distributions for test sequences
- `R/real_id_choices_anaylsis.R` — empirical analysis of patient response patterns and foreseeable rule-based answering strategies

This README intentionally excludes legacy analyses and unrelated Python or historical project files.

## Scientific role of the code

### 1) Theoretical score and reference-distribution analysis

The script `R/theoretical_distributions_and_scores.R` implements the methodological part of the paper concerned with representing an odor-identification test as a sequence of correct-answer positions and quantifying its structural complexity with information-theoretic measures.

It includes:

- classical chance distributions for forced-choice tests with $L$ trials and $k$ alternatives
- unconstrained and constrained random-sequence reference models
- Monte Carlo generation of the distribution of informativeness scores $I(X)$
- reference intervals and summary tables used to interpret empirical and theoretical values
- plots comparing observed or expected score distributions with random-choice baselines

### Sequence files for independent evaluation

The [Sequences](./Sequences) folder provides the first 100 unique sequences ranked by informativeness score $I$ for each available test design, or all available unique sequences when fewer than 100 were retained. These files are provided for readers who want to inspect the calculated sequences or independently evaluate possible response keys. They are named
`scores_L< L >_k< k >_top_sequences.txt`, where $L$ is the number of trials and $k$ is the number of response alternatives.

Each file contains the unconstrained sequences first. When a constrained equal-frequency analysis is possible, the constrained sequences follow after a heading such as:

```text
--- Top 100 unique constrained sequences by I score (L=18, k=6) ---
```

The constrained section is the relevant section when selecting a balanced test response key: every response position occurs equally often in these sequences. To use a file, select the file matching the desired $L$ and $k$, then read the sequence entries under the constrained heading. The digits in each sequence identify the response position selected on each trial.

An equal-frequency constraint can be applied only when $L$ is a multiple of $k$; a constrained section is included only for designs for which constrained sequences were retained. If no constrained section is present, the file contains only a heading such as:

```text
--- Top 100 unique unconstrained sequences by I score (L=18, k=6) ---
```

In that case, the listed sequences are still ranked candidates for independent evaluation, but they do not guarantee equal use of all response positions. The unconstrained sequences are also useful as a reference for comparing balanced and unrestricted sequence structures.

This component supports the paper's methodological framing: the sequence of correct-answer positions can be characterized in terms of balance and unpredictability, not only by overall score.

### Interpretation and design implications

The analysis treats the ordered sequence of correct response positions as a measurable property of a forced-choice test design, in addition to the usual test length $L$ and number of alternatives $k$. The block-entropy measures $H1$, $H2$, and $H3$ capture marginal balance and local sequence predictability; their weighted combination gives the informativeness score $I(X)$ and its per-trial form. The resulting candidates can therefore be ranked to identify balanced, locally unpredictable, and maximally or near-maximally informative response keys within a fixed test format.

This sequence analysis does not replace conventional clinical scoring. Participant performance remains the number of correctly identified odors, with its usual binomial chance interpretation. Informativeness describes the response key and the extent to which its structure avoids predictable positional patterns that could support non-olfactory response strategies.

For a fixed format, response keys can be generated and evaluated before administration, then selected from a predefined catalogue or deployed randomly from that catalogue. The established odor items, response format, and interpretation of participant scores remain unchanged. Reporting the implemented sequence, or its catalogue identifier, makes the administered key reproducible without requiring every future administration to use the same fixed sequence. The workflow also supports retrospective comparison of published test versions and comparison of candidate formats using a common information-based scale.

For reproducible reporting, describe the conventional test parameters together with response-position balance, $H1$, $H2$, $H3$, total or per-trial informativeness, relative efficiency when available, and the selected sequence or catalogue identifier. Shannon-based entropy measures, binomial chance probabilities, and Hamming distance serve complementary purposes: they characterize design structure, participant-score chance performance, and similarity to foreseeable response patterns, respectively.

### 2) Empirical analysis of real response choices

The script `R/real_id_choices_anaylsis.R` implements the empirical evaluation part of the paper. It converts observed patient answer strings into choice positions, summarizes response patterns, and compares them with random-choice expectations.

It includes:

- encoding answer strings into trial-wise choice positions
- compact summary of repeated pattern structures
- assessment of foreseeable deterministic or rule-based response behavior
- Hamming-distance analysis against a priori pattern templates
- reporting of pattern statistics relative to uniform random choice

This component is intended to evaluate whether observed patient response sequences show structure beyond what would be expected under random responding.

## How to run

From the repository root in R:

```r
source("R/theoretical_distributions_and_scores.R")
source("R/real_id_choices_anaylsis.R")
```

The scripts generate several CSV tables and plotting outputs in the working directory.

## Dependencies

Package dependencies are centralized in `R/globals.R`. The current analysis workflow relies on the following packages:

- `ggplot2`
- `ggthemes`
- `pbmcapply`
- `patchwork`
- `ComplexHeatmap`
- `cABCanalysis`
- `grid`

If needed, install missing packages before running the scripts:

```r
install.packages(c(
  "ggplot2", "ggthemes", "pbmcapply", "patchwork",
  "ComplexHeatmap", "cABCanalysis", "grid"
))
```

## A priori response-key template families

The empirical analysis in `R/real_id_choices_anaylsis.R` generates a catalogue of sequence families that are evaluated as foreseeable response-key patterns. These are the families used to formalize the idea that a participant may answer by following a positional rule rather than by odor perception. The following examples are representative of the generated catalogue in the code, using the 16-trial, 4-alternative format of the clinical task.

| Family | Example sequence | Interpretation |
| --- | --- | --- |
| Constant | 1111111111111111 | Same response position across all trials |
| Cyclic ascending | 1234123412341234 | Repeating cycle across the response positions |
| Cyclic descending | 4321432143214321 | Reverse repeating cycle |
| Doubled cyclic ascending | 1122334411223344 | Two-trial blocks of each position |
| Block ascending | 1111222233334444 | Large positional blocks |
| Two-alternative alternation | 1212121212121212 | Repeated alternation between two positions |
| Two constant halves | 1111111122222222 | One positional rule for the first half, another for the second |
| Zigzag ascending | 1234321234321234 | Alternating sweep with local reversals |
| Palindromic sweep | 1234432112344321 | Symmetric sweep-like pattern |

These are generated algorithmically in the code and can be cited in the paper as the a priori families used to separate conventional participant performance from the informativeness of the response-key sequence itself.

## Main outputs

Running the scripts produces reproducible summaries and result files such as:

- `chance_table_*.csv`
- `scores_*.csv`
- `scores_*_summary_I.csv`
- `scores_*_top_sequences.txt`
- `I_reference_distribution_summary_full.csv`
- `I_reference_distribution_summary_compact.csv`
- `foreseeable_pattern_results.csv`

## Reader and reviewer note

This repository is intended as a technical companion to the paper. It complements the manuscript by making the analytical workflow transparent, inspectable, and reproducible. It is not a substitute for the paper itself, which remains the authoritative source for scientific interpretation, discussion, and contextualization.

## Citation

Please cite the manuscript as:

Lötsch J, Hummel T, Himmelspach A. Information theoretic analysis of randomness and structure in clinical olfactory test design. (in preparation)

Use this repository as the corresponding computational companion to the paper and the implementation layer used to generate the reported analyses.