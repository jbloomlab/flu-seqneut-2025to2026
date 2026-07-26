# Does collapsing a strain's barcodes give better titers?

**NOTE: This analysis run on 2026-07-26. The decisive precision metric is unavailable on
this dataset (only one serum was measured on more than one plate), so this comparison is
weaker than it looks and should be re-run once more repeat sera exist.**

**This directory addresses one question: should the several barcodes carrying the same
viral strain be analyzed as separate measurements of that strain, or should their counts
be summed into a single measurement?** That is the
[`collapse_strain_barcodes`](../../seqneut-pipeline/README.md#collapse_strain_barcodes)
option in `seqneut-pipeline`. Almost every strain in this library carries 2 barcodes (a
few carry 3), so with the option off (the default) a plate gives 2-3 neutralization curves
per strain and their titers are averaged; with it on, the barcodes' counts are pooled
before any QC or curve fitting and a plate gives one curve per strain.

The two rationales pull in opposite directions, and the pipeline's own documentation
declines to make a universal recommendation:

- **Against collapsing:** several barcodes per strain are internal replicates. Averaging
  them reduces noise, and disagreement between them flags measurements that should not be
  trusted.
- **For collapsing:** the noise is largely set by how many virions of each barcoded strain
  land in each well, so pooling counts raises the effective input per strain and reduces
  counting noise. Pooling also weights each barcode by its depth, which is the
  statistically efficient way to combine them.

Which rationale wins is an empirical question about a particular library and dataset, so
it has to be measured rather than argued. The comparison here does that by running the
pipeline both ways on the same reads and comparing the two sets of titers.

## What it found

**The two settings produce all but indistinguishable titers, and on this dataset the
metric that would decide the question is unavailable. The choice is a trade-off, and it
is a low-stakes one.** In detail:

- The two settings give **nearly the same titers**: r = 0.989 on log2 over 26,470 matched
  uncensored pairs, typical disagreement 1.06-fold, 0.9% of titers differing by more than
  2-fold, 0.04% by more than 4-fold, none by more than 10-fold. Agreement is *closer* here
  than in the DRIVE comparison, on 4x as many titers.
- **No practically meaningful bias.** The median log2 fold-change is exactly zero; the
  mean is -0.022 (collapsing reads 0.985x, i.e. 1.5% lower). The Wilcoxon test is
  significant (p = 1.5e-4) purely because n = 26,470 — the magnitude is nil. Do not read
  the p-value as evidence of a scale difference.
- **The decisive precision metric could not be computed.** Plate-to-plate reproducibility
  needs sera measured on more than one plate, and only **1 of 307 sera** qualifies
  (`seqneut25to26_sera_4018`, on `plate16-2` and `plate25`), giving a single comparison
  against the 8 the script requires. The two runs came out at 0.788 vs 0.790 log2 spread —
  a 0.002 difference from one serum, which is noise, not a result. See the warning about
  the report's own wording in "Scope of these results" below.
- **Collapsing does fit visibly better curves** (median R² 0.961 vs 0.947; 23.0% vs 31.4%
  of curves below R² 0.9; median RMSD 0.061 vs 0.072), confirming that pooling counts
  genuinely reduces noise. But this is mechanism, not outcome: it does not demonstrably
  propagate into better titers, because per-barcode counts are already deep (median ~4,000
  per barcode-well, so Poisson noise is small for the typical measurement). Collapsing's
  *worst* fit is worse (min R² -1.22 vs -0.51).
- **The headroom is small either way.** Within-plate scatter between a strain's barcodes
  is 0.284 log2, against 0.43 log2 between plates. Because the plate-level term dominates
  and averaging does not reduce it, averaging the 2 barcodes per plate can cut total
  scatter by at most **7.8%** even in principle. A few-percent difference in precision
  would be fully explained by this alone.
- **Titer-matrix coherence is a dead heat**: on an identical 225 x 55 submatrix the first
  5 principal components explain 90.0% of variance without collapsing and 89.8% with it.
  This is the one quality metric here that uses all sera rather than only repeat ones, and
  it finds nothing.
- **What collapsing gives up** is per-titer uncertainty (99.7% of its titers end up with a
  single replicate, leaving `titer_sem` empty, against 4.8% without collapsing), the
  ability to detect bad measurements from replicate disagreement, and redundancy against a
  failed curve fit. This loss is *more* severe here than in DRIVE (99.7% vs 82%), because
  this library carries 2 barcodes per strain rather than 2-3.
- **Neither setting changes the downstream dataset's scope**: `results/final_titer_data/`
  retains 302 sera x 91 viruses either way.

Read `report.html` for the full picture, and the "Interpreting the results" section below
before drawing a conclusion — several tempting arguments in both directions do not hold up.
One caveat applies throughout: the collapsed run kept count thresholds that apply to
*summed* counts, making them roughly 2x looser, so part of its apparently cleaner QC is
that confound rather than the collapsing.

**Bottom line for this repository:** nothing here justifies switching away from the
default on quality grounds, and nothing here condemns it either. The one concrete,
non-speculative cost of collapsing on this dataset is the loss of `titer_sem` on
essentially every titer, which matters if any downstream analysis wants per-titer
uncertainty. The one concrete benefit is better-conditioned curve fits. Both are small.

## The tool

Although this directory holds one comparison, the script is general: it compares **any two
commits** of this repository, so the same machinery can be pointed at other analysis
choices — changed QC thresholds, a new pipeline version, added plates — without
modification. See "Reusing this on another comparison" below.

`compare_pipeline_runs.py` is a **byte-identical copy** of the script in the
[flu-seqneut-DRIVEstudy](https://github.com/jbloomlab/flu-seqneut-DRIVEstudy) repository's
directory of the same name (verified by md5). It was deliberately not modified, so that
the two comparisons are provably the same tool and their `summary.csv` files can be
stacked without worrying about drift. Two consequences worth knowing are recorded under
"Known artifacts" below.

## Usage

[`compare_pipeline_runs.py`](compare_pipeline_runs.py) compares the results of two
pipeline runs held at two git commits, and writes an HTML report plus machine-readable
tables. It reads the committed result CSVs out of each commit with `git show`, so it
re-runs nothing and never touches the working tree. This works because `.gitignore`
re-includes the key result CSVs, so past runs remain fully recoverable from history.

The command that produced the results in this directory:

```bash
D=non-pipeline_analyses/collapse_vs_noncollapse_barcode_comparison
python $D/compare_pipeline_runs.py \
    --run-a f584e886 --label-a no_collapse \
    --run-b 69656005 --label-b collapse \
    --comparison-id collapse_vs_noncollapse_barcode_comparison \
    --out $D --count-depth
```

Any commit-ish works (SHA, branch, tag), but **prefer explicit SHAs for a comparison you
intend to commit.** A branch name or `HEAD` silently re-points as work continues, so
re-running the "same" command later can compare something else — which is why the command
above is pinned rather than using the branch names the two runs were made on.

Outputs:

| file | what it is |
| --- | --- |
| `report.html` | the report; interactive charts, open in a browser |
| `summary.csv` | headline metrics in a fixed schema, safe to concatenate across comparisons |
| `tables/*.csv` | every table behind every number, including ones too long for the HTML |

`--count-depth` reads `results/barcode_counts/` to quantify how many barcode-wells are
count-limited, sampling 400 of the 2,681 committed count files. It is off by default
because those are thousands of small files; it was included here because with the
precision metric unavailable, section 7's bound on how much counting noise there is to
remove is one of the few well-powered mechanism checks left.

`--full-tables` would also write **bulk per-observation tables** (>= 10,000 rows). It was
not used, because those are near-verbatim copies of pipeline output already committed at
both commits. For this comparison the omitted tables are `concordance_per_titer.csv`
(26,470 rows) and `curvefit_per_curve.csv` (84,335 rows); both are named in the console
output and at the foot of `report.html`, never dropped silently.

Requires `pandas`, `numpy`, `scipy`, `altair`, `jinja2`, `pyyaml` — all present in the
`seqneut-pipeline` conda environment. The report's charts load the vega JavaScript from
a CDN, so viewing it needs a network connection; the data itself is embedded in the file.

## Scope of these results, and when to re-run

**The findings above are specific to two commits and to the data available when they were
run**, and are a snapshot rather than a standing conclusion:

| | |
| --- | --- |
| run A | `f584e886` ("update to `seqneut-pipeline` 8.0.0 and re-run"), `collapse_strain_barcodes: false` |
| run B | `69656005` ("re-run pipeline with *collapse_strain_barcodes: true*"), `collapse_strain_barcodes: true` |
| data | 28 plates, 307 sera, 91 strains at 2-3 barcodes each (mostly 2) |
| generated | 2026-07-26 |

The two commits are parent and child, and the only difference between them outside
`results/` and the rebuilt `auspice/` JSONs is the single added config line. The
`seqneut-pipeline` submodule pin is identical. This is therefore a clean single-variable
comparison, as the DRIVE one also was.

The exact commits are recorded in the report header and in every row of `summary.csv`, so
a stale report can always be told apart from a current one.

### Read section 2 of the report with suspicion

The report's section 2 prose is generated from templates that fire regardless of sample
size, and with a single comparison two of its statements are vacuous and read as far more
encouraging than the data warrants:

- *"Centre and tail agree, which is the straightforward case."* With one comparison the
  median and the 90th percentile are **the same single number**, so agreement is
  arithmetic, not evidence.
- *"Both strata point the same way (1 and 1 comparisons), so the result does not depend on
  whether QC-dropped pairs are included. That is a meaningfully stronger finding than
  either stratum alone."* It is not a stronger finding; it is the same lone serum counted
  twice. (Note that in DRIVE, with 33 comparisons, the two strata genuinely *disagreed*
  and the sign of the result followed the stratum — so this is exactly the check that
  cannot be done here.)

The report does also flag the stratum as `__UNDERPOWERED` and prints the "only 1 serum x
plate-pair comparisons ... should not be interpreted" warning. Believe the warning, not
the prose.

**It is worth re-running as more data accumulates.** The reason this comparison cannot
settle the question is entirely the repeat-sera shortage: **1 of 307 sera** was measured
on more than one plate, versus 11 of 62 in DRIVE. That kills sections 2 (reproducibility),
3 (adjudication of run-unique titers) and the between-plate half of 4a, leaving the
comparison resting on 1, 4b, 5, 6 and 7 — none of which observes the final titers'
precision directly. Every additional serum run on multiple plates tightens this directly,
and a handful would restore the decisive metric.

Re-running is also worthwhile if any of the following changes, each of which could move
the answer rather than merely sharpen it:

- **More repeat sera**, for the reason above. This is overwhelmingly the main one here.
- **A different viral library**, especially one with more barcodes per strain or a more
  uneven barcode abundance distribution. The case for summing counts gets stronger as
  barcode counts get shallower or more skewed, so a conclusion drawn on this library does
  not transfer automatically.
- **Changed count or goodness-of-fit thresholds.** The comparison checked in here carries
  a known confound: the collapsed run kept the same count thresholds, which apply to
  summed counts and so run roughly 2x looser. A run with those rescaled would separate
  "collapsing helps" from "collapsing loosened QC". This was deliberately not done here.
- **A new `seqneut-pipeline` version** that changes QC stages or how replicates are
  combined.

Because `summary.csv` uses a schema that does not vary between comparisons, successive
re-runs stack into one table and can be compared directly — see the section on that below.
Keeping the old comparison directory alongside a new one is the intended workflow.

## Known artifacts

Two things to know before reading the report's numbers literally. Both follow from keeping
the script byte-identical to the DRIVE copy, and both were judged not worth diverging the
script over; neither changes any conclusion.

- **Section 7 reports the *designed* library, not the one behind the titers.** The script
  merges every library listed under `viral_libraries` in `config.yml`, and here the
  designed library is a superset of the actual one, so `library_structure.csv` says 235
  barcodes / 114 strains where the neutralization plates actually used 186 barcodes / 91
  strains. Median barcodes per strain is 2 either way, so the section's argument is
  unaffected, but the distribution table hides the two strains that carry a single barcode
  in the actual library (for which collapsing is a no-op). The same merge makes the
  `virus_absent_from_run` check in section 1a a superset test, which is harmless — it
  cannot produce a false "absent" attribution.
- **`median_between_plate_sd_log2` in section 4a comes from one serum.** It is reported
  with n = 91, but those are 91 strains from the *single* repeat serum's single plate pair,
  not 91 independent estimates. It is also uncentered, unlike section 2's
  `sd_log2_centered`, so the gap between the runs (0.434 vs 0.475) is mostly the differing
  plate offset rather than differing precision — which is why it points the opposite way
  from the RMSD on the very same plate pair (0.914 vs 0.910). None of it is meaningful at
  this sample size. The within-plate figure (0.284, n = 25,480) is solidly estimated and is
  the half of 4a worth trusting.

## Why this is not just a diff

Four ways to compare two runs are all confounded, and the script is shaped mostly by
avoiding them. This matters because each one has an intuitive reading that points the
wrong way.

**Drop counts are not a quality signal on their own.** The per-serum QC check
`max_fold_change_from_median` compares each replicate against the median over
replicates. A titer with a single replicate has a fold-change of exactly 1 and therefore
*cannot fail the check*. Collapsing strain barcodes takes 99.7% of titers here down to one
replicate, so its per-serum drop count of 0 (against 21 without collapsing) is not
evidence of cleaner data — it is the check being unable to fire. The script always reports
drop counts next to the fraction of single-replicate titers that makes them interpretable.

**Statistics computed on each run's retained titers are not comparable between runs,**
because the runs retain different sets, and a run that filters harder scores better on
what survives. Every quality metric is computed on a *matched set* with each run's
per-serum QC ignored. The "as shipped" numbers appear too, labelled as not comparable.

**Restricting to the two runs' shared *final* titers is itself biased** — the subtlest of
the four. Conditioning on both runs' per-serum QC having passed deletes the pairs whose
replicates disagreed, which are exactly the pairs a replicate-rich run rejects and a
replicate-poor run keeps. So it favours whichever run filters more effectively, by removing
the cases where the other approach might have won. Section 2 therefore reports
**`pre_qc_union`** (every pair with per-plate data in both runs, no QC applied) as its
primary stratum and keeps `final_matched` as secondary. Here `pre_qc_union` adds back 331
pairs that `final_matched` excludes — but with only one serum contributing, the comparison
between the strata is not informative on this dataset.

**Count-based QC thresholds change meaning** when `collapse_strain_barcodes` differs:
`avg_barcode_counts_per_well` and friends apply to summed counts once barcodes are
collapsed, so they are effectively looser by roughly the barcodes-per-strain factor —
about 2x here. The script detects this from the two configs and prints a standing caveat.
It cannot correct for it; that needs a third run with rescaled thresholds, which was not
done.

## What the report contains

1. **Coverage and concordance** — which serum-virus titers exist in which run, why each
   missing one is missing, and how similar the shared titers are. On this dataset there
   are **no scope differences at all**: both runs cover the same sera, strains and plate
   sets, so the entire coverage difference (63 titers only without collapsing, 40 only
   with it, out of a 27,904 union) is attributable to QC.
2. **Plate-to-plate reproducibility** — normally the decisive precision metric.
   **Effectively unavailable here**: 1 serum, 1 plate pair. See the warning above.
3. **Adjudication of run-unique titers** — do the titers only one run reports actually
   replicate? Also unavailable: only 1 unique titer has multi-plate data.
4. **Variance components and curve-fit quality** — the mechanism behind any difference.
   4b (curve-fit quality) is well powered here and is the clearest signal in the report;
   4a is half-usable, see "Known artifacts".
5. **Retention accounting by QC stage.**
6. **Titer-matrix coherence** — a replicate-free noise proxy that uses every serum, not
   only the repeat ones. Relatively more important here than in DRIVE, precisely because
   it does not need repeat sera.
7. **Dataset context** — library structure and count depth, which bound the possible
   benefit of collapsing.

### Two definitions worth knowing

**The matched set** is serum-virus pairs where both runs measured the same serum on the
*same set of plates*, per-serum QC ignored. Requiring an identical plate set matters as
soon as the two commits differ in which plates they include: a serum measured on three
plates in one run and two in the other has different amounts of averaging behind its
titer, and comparing them measures that instead of the change under test. Here the two
commits have identical plate coverage, so nothing is lost to this requirement.

**Reproducibility is reported as `sd_log2_centered`**, the spread of log2 titer
differences between two plates after removing their median offset — not as Pearson *r*.
A plate-wide dilution offset is a nuisance shared by every strain on the plate and is not
the noise under study, so it is removed. And *r* depends on the spread of true titers
across strains, so a run that slightly compresses dynamic range scores worse without
being any less precise. `pearson_r` is still reported, for continuity with the pipeline's
own `plate_to_plate_corr` output.

## `summary.csv` and comparing across datasets

`summary.csv` has one row per metric in a schema that does not vary between
comparisons:

```
comparison_id, run_a, run_b, run_a_commit, run_b_commit, section, metric, stratum,
value_a, value_b, delta, n, test, statistic, p_value, higher_is_better, note
```

Because the schema is fixed, summary files from several comparisons stack into one table
and can be analysed together:

```python
import pandas as pd, glob
all_comparisons = pd.concat(
    [pd.read_csv(f) for f in glob.glob("non-pipeline_analyses/*comparison*/**/summary.csv", recursive=True)]
)
```

**Check `stratum` for a `__UNDERPOWERED` suffix and exclude those rows.** Restricted
strata can shrink to a handful of comparisons where the apparent difference is driven by
which sera happened to survive the restriction; they are kept for completeness, not for
interpretation. On this dataset that flag is doing real work — it is what marks section 2
as uninterpretable.

### Stacked against the DRIVE comparison

This is the second dataset the same question has been asked on, so the cross-dataset read
the schema exists for is now possible. Stacking this `summary.csv` with the one in
`flu-seqneut-DRIVEstudy` (schemas verified identical) gives, for A = no_collapse and
B = collapse:

| metric | DRIVE (A -> B) | this dataset (A -> B) | consistent? |
| --- | --- | --- | --- |
| median R² of curve fits | 0.969 -> 0.979 | 0.947 -> 0.961 | yes, collapsing better |
| fraction of curves R² < 0.9 | 24.3% -> 17.9% | 31.4% -> 23.0% | yes, collapsing better |
| median RMSD of curve fits | 0.066 -> 0.056 | 0.072 -> 0.061 | yes, collapsing better |
| PC1-5 cumulative variance | 0.938 -> 0.932 | 0.900 -> 0.898 | yes, no_collapse better, negligibly |
| fraction single-replicate titers | 4.0% -> 82.3% | 4.8% -> 99.7% | yes, worse here |
| predicted max gain from averaging | 5.1% | 7.8% | yes, single-digit % both |
| median \|log2 FC\| between settings | 0.129 (1.09x) | 0.081 (1.06x) | agreement closer here |
| Pearson r between settings | 0.975 | 0.989 | agreement closer here |
| plate-to-plate `sd_log2_centered` | 0.750 -> 0.718, n=33, p=0.57 | n=1, no test | **only DRIVE informative** |
| barcodes per strain | 3 | 2 | less redundancy here |

**Every direction that can be compared replicates.** Collapsing consistently fits better
curves, and the headroom for replicate averaging is single-digit percent in both. The
coherence rows also share a sign, but read nothing into that: the gaps (0.16 and 0.64
percentage points) carry no significance test, sign agreement across two datasets is
half-likely by chance, and the metric structurally tilts toward the non-collapsing run
anyway — see the caveat below. The one metric that
actually measures the precision of the final titers is informative only in DRIVE, where it
was a coin flip (collapsing won 16 of 33 comparisons, p = 0.57).

So the two datasets together say the same thing more confidently than either alone: the
mechanism-level advantage of collapsing is real and reproducible, and it does not visibly
reach the titers. Notably this holds despite this library being the less favourable case
for collapsing on paper — 2 barcodes per strain rather than 3, and counts deep enough
(median ~4,000 per barcode-well, 6.5% below 50) that there is little Poisson noise for
summing to remove.

## Reusing this on another comparison

The script takes any two commit-ishes and does not assume the difference between them is
`collapse_strain_barcodes`. Group names, plate names, sera and strains all come from the
data. Each analysis is isolated, so an impossible or failing one becomes a visible note
in the report rather than a crash, and the report lists what was skipped and why.

Grep the script for `ADAPT:` for every tunable knob and `ASSUMPTION:` for every place
something about the current dataset is taken for granted. All of them were checked against
this repository before the run above; the ones that mattered:

- **Repeat sera.** Sections 2, 3 and 4a need sera measured on more than one plate. Check
  `n_repeat_sera` in the report header before trusting any precision claim. **This is what
  went wrong here**, and it is the first thing to check on any new dataset.
- **QC-drop stage keys.** `BARCODE_STAGE_DROP_KEYS` and `CURVEFIT_STAGE_DROP_KEYS` list
  the `qc_drops.yml` keys the retention accounting knows about. All five keys written at
  both commits are covered, so section 5 misses nothing. A new pipeline QC stage would
  need adding there.
- **Threshold uniformity.** `_find_threshold` returns the first matching value and so
  assumes plates share defaults via YAML anchors. Verified: all 28 plates carry identical
  count thresholds, so the caveat it prints quotes the right numbers.
- **Titer definition.** The script refuses to compare titers across a change in
  `titer_as` or `titer_units`. Both runs use the same values, so no issue.
- **Pipeline schema.** Both commits are pipeline 8.0.0 and write the `plate` and
  `dropped_by_qc` columns and all three `qc_drops` YAMLs, so none of the
  older-pipeline shims fire and no quantity is reported as unknown.
- **Groups.** This dataset has a single group (`human`), where DRIVE had two, so the
  per-group breakouts collapse to one row. Nothing assumes more than one group.
- **Barcodes per strain.** With one barcode per strain, collapsing is a no-op and the
  comparison cannot show anything; the script warns. Here it is 2-3, which is enough to
  compare but is the low end — see "Known artifacts" on how the library is reported.
- **Censored titers.** Titers at a bound are excluded from log-scale statistics and
  counted separately. Collapsing produces somewhat fewer of them (722 upper / 287 lower
  vs 970 / 335), consistent with pooled counts slightly extending the measurable range.

## Interpreting the results

**Expect a trade-off, not a winner.** Two runs of the same assay on the same reads rarely
differ enough for one to dominate. What differs is *which failure mode each is exposed
to*. Report the trade-off rather than picking a side the data does not support.

Trust the metrics in roughly this order: plate-to-plate reproducibility on
`pre_qc_union` (2), adjudication of unique titers (3), variance components and curve-fit
quality (4), titer-matrix coherence (6). **On this dataset the first two are unavailable**,
so the ordering effectively starts at (4) and (6) — which is the same as saying this
dataset cannot answer the question, only bound how much it could possibly matter. Treat as
*not* quality signals on their own: raw counts of dropped titers (5), the number of curves
fit (4b), and anything computed on each run's own retained set.

Four checks before believing a section 2 result:

1. **Is it significant?** If not, the conclusion is "no difference detected", *not* "no
   difference". Here no test could even be run. Also read
   `n_pairs_where_b_more_reproducible`: a near coin-flip alongside a non-significant
   p-value is much stronger evidence of a genuine wash than the p-value alone.
2. **Do centre and tail agree?** Compare the median against the `_p75`/`_p90` rows. Two
   runs can share a median spread while one has clearly worse worst cases, and it is the
   worst cases that mislead a downstream analysis. With one comparison, as here, these are
   the same number and the check is vacuous.
3. **Do `pre_qc_union` and `final_matched` agree?** If yes, the finding is robust to how
   the set was chosen. If no, the metric definition is deciding the comparison and the
   right conclusion is that neither approach is clearly better. Also vacuous at n = 1.
4. **What does section 4 say the headroom is?** If the predicted gain from replicate
   averaging is a few percent, then a few-percent gap in reproducibility is fully
   explained by that alone and implies nothing about counting noise. Here it is 7.8%.

### Arguments that sound decisive but are not

Cutting in both directions, because it is easy to marshal only the ones supporting a
preferred conclusion:

- *"It drops fewer titers, so it is cleaner."* Fewer drops can mean less detection — see
  confound 1. The 21 -> 0 change in per-serum drops here is almost entirely the check
  becoming unable to fire.
- *"It loses the ability to detect replicate disagreement."* This partly begs the
  question. Under collapsing, barcode disagreement is not an *undetected* error but an
  *eliminated* one, absorbed into a single better-powered measurement. The objection only
  bites if barcode disagreement signals something else wrong — a mis-assigned barcode,
  contamination — which is what section 3 tests, and section 3 is unavailable here.
- *"Averaging replicates must be more precise."* Only for the part of the noise that is
  independent between replicates. Section 4a's predicted 7.8% gain is an **upper** bound:
  it models within-plate scatter as independent, so to the extent barcode differences are
  systematic, averaging helps even less.
- *"Summing counts must be more precise."* Only where counts limit precision. Section 7 is
  the check: at a median of ~4,000 counts per barcode-well with 6.5% below 50, Poisson
  noise is already negligible for the typical measurement and the gain is confined to the
  low-abundance tail.
- *"The titer matrix is more coherent without collapsing, so collapsing adds noise."* The
  gap is 0.16 percentage points with no significance test attached, which is not evidence
  of anything. Worse, the metric is not neutral between the two settings: it rewards
  concentrating variance in the leading components, and averaging a strain's barcode
  titers is itself a smoothing operation, so the non-collapsing run can score higher by
  smoothing away genuine strain-specific signal rather than by carrying less noise. The
  script's own section-6 note flags this failure mode. Treat section 6 here as finding
  nothing, not as a narrow win.
- *"Collapsing fits better curves, so it gives better titers."* The better fits are real
  and replicate across both datasets, but a curve fit is an intermediate, not the answer.
  Nothing here shows the improvement surviving into the titers, and section 2 — which
  would have shown it — could not be run.
- *"More curves fit, or more replicates, is more data."* Curve and replicate counts are
  design consequences, not quality measures. The 56,403 -> 27,932 drop in curve count is
  arithmetic, not loss.
- *"The looser count thresholds under collapsing inflate its retention."* True, but it
  cuts both ways: it also means collapsing was evaluated with *inappropriately lax* count
  QC and still produced titers indistinguishable from the alternative.

### The two asymmetries that do survive

- **No redundancy against a failed fit.** One curve per strain means a single bad fit
  removes the titer outright, where several replicates still yield a titer from the
  survivors. `lost_at_curvefit_qc` and `lost_at_barcode_or_well_qc` in section 1a size
  this; here all 63 titers missing from the collapsed run were lost before any curve was
  fit, and collapsing's worst individual fit is worse than uncollapsed's (min R² -1.22 vs
  -0.51).
- **Count-weighting is the efficient combination.** Summing counts is inverse-variance
  weighting; averaging titers across barcodes weights a shallow barcode equally with a
  deep one. This favours collapsing, grows as barcode abundance gets more uneven, and none
  of the metrics here isolates it directly — a gap worth closing by stratifying section 2
  on per-strain barcode evenness if it becomes decision-relevant.

Also note that a run with mostly single-replicate titers has an empty `titer_sem` and a
dead `min_replicates` lever. On this dataset that applies to 99.7% of the collapsed run's
titers, which is the most concrete consequence of the choice found anywhere in this
comparison — it matters if downstream analyses want per-titer uncertainty regardless of
how the precision comparison lands.
