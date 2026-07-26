#!/usr/bin/env python3
"""Compare the results of two ``seqneut-pipeline`` runs held at two git commits.

The pipeline's key result CSVs are committed to this repository (see ``.gitignore``,
which excludes ``results/**`` but re-includes the important files), so two runs can be
compared by reading those files out of two commits with ``git show``. Nothing is
re-run, and no working tree is touched.

    python <this-dir>/compare_pipeline_runs.py \
        --run-a 417764a --label-a no_collapse \
        --run-b 963ae9c --label-b collapse \
        --comparison-id my_comparison

Prefer explicit SHAs over branch names or ``HEAD`` for a comparison you intend to commit:
a branch silently re-points as work continues, so re-running the "same" command later can
compare something else.

Outputs, under ``--out`` (default: a ``<comparison-id>/`` subdirectory beside this script):

    report.html   self-contained report, interactive Altair charts
    summary.csv   one row per metric, stable schema, safe to stack across comparisons
    tables/*.csv  every table behind every number in the report

======================================================================================
WHY THE COMPARISON IS BUILT THE WAY IT IS  (read this before changing any metric)
======================================================================================

The three most obvious ways to compare two runs are all confounded, and this script is
shaped mostly by avoiding them:

1.  "Which run drops fewer serum-virus titers?" is not a quality signal on its own.
    The per-serum QC check ``max_fold_change_from_median`` compares each replicate to
    the median over replicates. A run with one replicate per titer has
    ``fc_from_median == 1`` by construction, so the check *cannot fire*. Collapsing
    strain barcodes takes ~82% of titers in this dataset down to one replicate, so it
    "drops fewer titers" mostly by losing the ability to detect bad ones. Any drop
    count is therefore reported next to the replicate distribution that makes it
    interpretable, never on its own.

2.  Any statistic computed on each run's *retained* titers is not comparable between
    runs, because the runs retain different sets. A run that filters more aggressively
    gets a better-looking correlation on what is left. Every quality metric here is
    therefore computed on a MATCHED set (defined below) with each run's per-serum QC
    ignored; the "as shipped" version is reported too, but always labelled secondary.

3.  Count-based QC thresholds change meaning between runs that differ in
    ``collapse_strain_barcodes``. ``avg_barcode_counts_per_well`` and
    ``min_no_serum_count_per_viral_barcode_well`` apply to summed counts once barcodes
    are collapsed, so they are effectively looser by roughly the number of barcodes per
    strain. This script cannot correct for that; it detects the situation from the two
    configs and prints a standing caveat. Removing the confound needs a third run with
    rescaled thresholds.

THE MATCHED SET. Quality metrics use serum-virus pairs where both runs measured the
same serum on the same set of plates, ignoring per-serum QC. Requiring the same plate
set matters as soon as the two commits differ in which plates they include: otherwise a
serum measured on 3 plates in one run and 2 in the other contributes a different amount
of averaging to each, and the comparison silently measures that instead of the change
under test.

4.  A fourth confound, found the hard way and the subtlest of the lot: restricting a
    quality metric to pairs both runs report in their FINAL titers conditions on both
    runs' per-serum QC having passed, which is not neutral. It deletes the pairs whose
    replicates disagreed -- exactly the pairs a replicate-rich run rejects and a
    replicate-poor run keeps -- and so favours whichever run filters more effectively by
    removing the cases where the other approach might have won. Section 2 therefore
    reports `pre_qc_union` (every pair with per-plate data in both runs, no QC applied)
    as its PRIMARY stratum and `final_matched` as secondary, and says so loudly when the
    two disagree. On this repository's collapse-vs-nocollapse comparison they did
    disagree, and the sign of the result followed the stratum rather than the data.

REASON HIERARCHY FOR A MISSING TITER. Absence is attributed to the first cause that
applies, because only the last three are quality signals at all:

    1. serum_absent_from_run       serum not in that run's groups_sera_by_plate.csv
    2. virus_absent_from_run       strain not in that run's viral library
    3. plate_set_differs           serum in both runs but measured on different plates
    4. lost_at_barcode_or_well_qc  no curve was ever fit
    5. lost_at_curvefit_qc         curve fit but failed goodness-of-fit
    6. dropped_by_per_serum_qc     titer existed, failed serum-level QC

Causes 1-3 are scope differences between the commits, not evidence about either
approach, and are excluded from every quality metric.

======================================================================================
INTERPRETING THE RESULTS
======================================================================================

Expect a trade-off, not a winner. Two runs of the same assay on the same reads rarely
differ enough for one to dominate; what differs is WHICH failure mode each is exposed to.
State the trade-off rather than picking a side the data does not support.

Read the sections in this order, and resist stopping early:

*   Section 2 `pre_qc_union` is the headline. If the paired p-value is not significant,
    the honest conclusion is "no difference detected", NOT "no difference". With a
    handful of repeat sera only a large effect is detectable at all.
*   Check the median against the p75/p90 rows. Two runs can share a median spread while
    one has clearly worse worst cases, and it is the worst cases that mislead a
    downstream analysis. When centre and tail disagree, which one matters is a judgement
    about downstream use, not something the data settles.
*   Check whether `pre_qc_union` and `final_matched` agree. If they do, the finding is
    robust. If not, the metric definition is deciding the comparison and the right
    conclusion is that neither approach is clearly better.
*   Read section 4 before believing section 2. If the predicted gain from replicate
    averaging is a few percent, then a few-percent gap in reproducibility is fully
    explained by that alone and implies nothing about counting noise.

Arguments that sound decisive but are not, in both directions:

*   "It drops fewer titers, so it is cleaner." No: see confound 1. Fewer drops can mean
    less detection.
*   "It loses the ability to detect replicate disagreement." This partly begs the
    question. Under collapsing, barcode disagreement is not an undetected error but an
    eliminated one, absorbed into a single better-powered measurement. The objection only
    bites if barcode disagreement signals something else wrong -- a mis-assigned barcode,
    contamination -- which section 3 tests and this script does not otherwise assume.
*   "Averaging replicates must be more precise." Only for the part of the noise that is
    independent between replicates. Section 4a's predicted gain is an UPPER bound: it
    models within-plate scatter as independent, so to the extent barcode differences are
    systematic, averaging helps even less than predicted.
*   "Summing counts must be more precise." Only where counts limit precision. Section 7
    is the check: with a median of thousands of counts per barcode-well, Poisson noise is
    already negligible for the typical measurement and the gain is confined to the
    low-abundance tail.
*   "More curves fit / more replicates is more data." Curve counts are a design
    consequence, not a quality measure.

One genuinely asymmetric consideration, easy to miss: a run that fits one curve per
strain has no redundancy against a failed fit, so one bad fit removes the titer outright,
where a run with several replicates still reports a titer from the survivors. Look at
`lost_at_curvefit_qc` and `lost_at_barcode_or_well_qc` in section 1a for its size. Against
that, summing counts is inverse-variance weighting, which is the statistically efficient
combination, whereas averaging titers across barcodes weights a shallow barcode equally
with a deep one -- a real advantage for collapsing that grows as barcode abundance gets
more uneven, and that none of the metrics here isolates directly.

======================================================================================
ADAPTING THIS TO ANOTHER COMPARISON OR DATASET
======================================================================================

Grep for ``ADAPT:`` for every knob, and ``ASSUMPTION:`` for every place where something
about this dataset is being taken for granted. The script is written to degrade
gracefully: each analysis is wrapped so that a failure or an impossible analysis
becomes a note in the report rather than a crash, and the report states which sections
were skipped and why.

The things most likely to need attention on a new dataset:

*   REPEAT SERA. Sections 2, 3 and 5 need sera measured on more than one plate. This
    dataset has 11 of 62 sera on 3 plates each; everything else is single-plate. With
    zero repeat sera those sections are skipped entirely and the comparison rests on
    sections 1, 4, 6 and 7. Check ``n_repeat_sera`` in the report header before
    trusting any reproducibility claim, and note that the repeat sera may not be
    representative of the rest.

*   DIFFERENT PLATE SETS between commits. Handled by the reason hierarchy and the
    matched-set definition, but read the coverage table: if ``plate_set_differs`` or
    ``serum_absent_from_run`` are large, the two runs are not really measuring the same
    experiment and the headline concordance number means little.

*   BARCODES PER STRAIN. If the library has one barcode per strain, collapsing is a
    no-op and this whole comparison is meaningless; the script warns. The
    noise-reduction argument for collapsing scales with how *uneven* the barcode
    distribution is, not just how many barcodes there are, so section 7 reports
    evenness rather than only the count.

*   QC DROP KEYS ARE NOT THE SAME OBJECT IN THE TWO RUNS. With barcodes collapsed the
    pipeline names the collapsed barcode after its strain, so ``qc_drops.yml`` keys are
    strain names; without collapsing they are barcode sequences. ``_normalise_to_strain``
    maps both onto strains using the viral library. If a future comparison changes
    something else about barcode naming, that function is where to fix it.

*   GROUPS. Group names (here ``G4``, ``TestPlate``) are never hardcoded; they come from
    the data. Groups differ wildly in design in this dataset -- the ``TestPlate`` sera
    are all single-plate -- so most tables are also broken out per group.

*   TITER DEFINITION. ``ASSUMPTION:`` the two runs use the same ``titer_as``
    (midpoint vs nt50) and the same ``titer_units``. The script checks and refuses to
    compare titers across a change in either, since the numbers would not be on the
    same scale.

*   CENSORED TITERS. Titers at a bound (``titer_bound`` of ``upper``/``lower``) are
    excluded from all log-scale statistics and counted separately. If a future dataset
    has many censored titers, consider a survival-style treatment instead; a simple
    exclusion biases toward the measurable middle of the range.
"""

from __future__ import annotations

import argparse
import base64
import functools
import io
import json
import subprocess
import sys
import textwrap
import traceback
from dataclasses import dataclass, field
from pathlib import Path

import altair as alt
import numpy
import pandas as pd
import yaml
from jinja2 import Environment
from scipy import stats

# ======================================================================================
# Knobs. ADAPT: these are the thresholds and choices most worth revisiting per dataset.
# ======================================================================================

# Fold-change bins used when summarising how much two runs' titers disagree. Reported as
# fractions of the matched set, so they are comparable across datasets of any size.
FOLD_CHANGE_BINS = [2, 4, 10]

# Number of principal components reported in the titer-matrix coherence section.
N_PCS = 5

# A serum needs at least this many plates to contribute to any reproducibility metric.
MIN_PLATES_FOR_REPRODUCIBILITY = 2

# A table with at least this many rows is treated as bulk per-observation output and is
# written only under ``--full-tables``. Such tables are typically near-verbatim copies of
# pipeline output that is already committed at both commits, so they dominate the on-disk
# size while adding little that ``git show`` cannot reproduce. ADAPT: raise or lower this
# if a dataset's useful tables sit either side of it; anything omitted is always named in
# the console output and the report footer, never dropped silently.
BULK_TABLE_MIN_ROWS = 10_000

# Minimum number of serum x plate-pair comparisons before a reproducibility statistic is
# reported at all. Restricted strata (notably the discordant tail) can shrink to a
# handful of comparisons, where the difference between two runs is dominated by which
# sera happened to survive the restriction. Such strata are written to the tables but
# their summary rows are flagged as underpowered so a meta-analysis can exclude them.
MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM = 8

# ADAPT: the pipeline writes plate-level QC drops under these keys. The first group
# happens before any curve is fit, the second is the curve-fitting QC. If the pipeline
# gains a new QC stage, add its key here or it will be silently unattributed.
BARCODE_STAGE_DROP_KEYS = ("wells", "barcodes", "barcode_wells")
CURVEFIT_STAGE_DROP_KEYS = ("barcode_serum_replicates", "serum_replicates")

# Titers not at a bound. ASSUMPTION: the pipeline writes exactly this string for
# uncensored titers.
UNCENSORED = "interpolated"

KEY = ["group", "serum", "virus"]


# ======================================================================================
# Reading a run out of a commit
# ======================================================================================


class Git:
    """Read files out of a commit without touching the working tree."""

    def __init__(self, repo: Path):
        self.repo = repo

    def _run(self, *args: str) -> str:
        proc = subprocess.run(
            ["git", "-C", str(self.repo), *args],
            capture_output=True,
            text=True,
        )
        if proc.returncode:
            raise RuntimeError(f"git {' '.join(args)} failed:\n{proc.stderr}")
        return proc.stdout

    @functools.lru_cache(maxsize=None)
    def resolve(self, commit: str) -> tuple[str, str]:
        """Return the (abbreviated sha, subject line) of a commit-ish, for provenance."""
        sha = self._run("rev-parse", "--short", commit).strip()
        subject = self._run("log", "-1", "--format=%s", commit).strip()
        return sha, subject

    @functools.lru_cache(maxsize=None)
    def ls_tree(self, commit: str) -> tuple[str, ...]:
        return tuple(self._run("ls-tree", "-r", "--name-only", commit).splitlines())

    @functools.lru_cache(maxsize=None)
    def show(self, commit: str, path: str) -> str | None:
        """File contents at a commit, or None if the file is not there."""
        try:
            return self._run("show", f"{commit}:{path}")
        except RuntimeError:
            return None

    def read_csv(self, commit: str, path: str) -> pd.DataFrame | None:
        text = self.show(commit, path)
        if text is None:
            return None
        return pd.read_csv(io.StringIO(text))

    def read_yaml(self, commit: str, path: str):
        text = self.show(commit, path)
        if text is None:
            return None
        return yaml.safe_load(io.StringIO(text))

    def concat_csvs(self, commit: str, paths, **extra_cols) -> pd.DataFrame:
        """Read and concatenate many CSVs, skipping any that are absent."""
        frames = []
        for path in paths:
            frame = self.read_csv(commit, path)
            if frame is None or frame.empty:
                continue
            for col, func in extra_cols.items():
                frame[col] = func(path)
            frames.append(frame)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)


@dataclass
class Run:
    """Everything read out of one commit."""

    label: str
    commit: str
    sha: str
    subject: str
    config: dict
    titers: pd.DataFrame  # aggregated, post-QC: the run's final answer
    per_replicate: pd.DataFrame  # one row per replicate, pre-per-serum-QC
    curvefits: pd.DataFrame  # one row per fitted curve, with r2/rmsd
    plate_to_plate: pd.DataFrame  # the pipeline's own correlation CSVs
    sera_by_plate: pd.DataFrame
    plate_qc_drops: dict
    groups_sera_qc_drops: dict
    barcode_to_strain: dict
    notes: list[str] = field(default_factory=list)

    @property
    def collapsed(self) -> bool:
        value = self.config.get("collapse_strain_barcodes", False)
        if isinstance(value, str):
            return value.strip().lower() in ("true", "yes", "1")
        return bool(value)

    @property
    def strains(self) -> set[str]:
        return set(self.barcode_to_strain.values())

    @property
    def sera(self) -> set[tuple[str, str]]:
        if self.sera_by_plate.empty:
            return set()
        return set(map(tuple, self.sera_by_plate[["group", "serum"]].values))

    @property
    def plates_by_serum(self) -> dict[tuple[str, str], frozenset]:
        """(group, serum) -> the set of plates that serum was measured on."""
        if self.sera_by_plate.empty:
            return {}
        return {
            (row.group, row.serum): frozenset(str(row.plates).split(";"))
            for row in self.sera_by_plate.itertuples()
        }

    @property
    def n_barcodes_per_strain(self) -> pd.Series:
        if not self.barcode_to_strain:
            return pd.Series(dtype=int)
        return pd.Series(self.barcode_to_strain).value_counts()


def load_run(git: Git, commit: str, label: str) -> Run:
    """Read every file this comparison needs out of one commit."""
    sha, subject = git.resolve(commit)
    tree = git.ls_tree(commit)

    config = yaml.safe_load(git.show(commit, "config.yml") or "{}") or {}

    # The aggregated titers are the run's final answer: one file per group.
    titers = git.concat_csvs(
        commit, [p for p in tree if p.startswith("results/aggregated_titers/titers")]
    )

    # Per-replicate titers carry the plate and replicate structure, and include rows
    # that per-serum QC later dropped (flagged by `dropped_by_qc`). This is the only
    # place the pre-QC picture survives, so most quality metrics are built from it.
    per_replicate = git.concat_csvs(
        commit,
        [p for p in tree if p.endswith("/titers_per_replicate.csv")],
    )

    curvefits = git.concat_csvs(
        commit, [p for p in tree if p.endswith("/curvefits.csv")]
    )

    plate_to_plate = git.concat_csvs(
        commit,
        [p for p in tree if p.startswith("results/plate_to_plate_corrs/")],
    )

    sera_by_plate = git.read_csv(commit, "results/sera/groups_sera_by_plate.csv")
    if sera_by_plate is None:
        sera_by_plate = pd.DataFrame(columns=["group", "serum", "plates"])

    plate_qc_drops = git.read_yaml(commit, "results/qc_drops/plate_qc_drops.yml") or {}
    groups_sera_qc_drops = (
        git.read_yaml(commit, "results/qc_drops/groups_sera_qc_drops.yml") or {}
    )

    # The viral library maps barcodes to strains. Needed both to know which strains a
    # run could have measured at all, and to normalise QC drop keys, which are barcode
    # sequences without collapsing and strain names with it.
    barcode_to_strain: dict[str, str] = {}
    notes: list[str] = []
    for lib_name, lib_path in (config.get("viral_libraries") or {}).items():
        lib = git.read_csv(commit, lib_path)
        if lib is None:
            notes.append(f"viral library {lib_name!r} ({lib_path}) not found at {label}")
            continue
        if not {"barcode", "strain"} <= set(lib.columns):
            notes.append(f"viral library {lib_name!r} lacks barcode/strain columns")
            continue
        barcode_to_strain.update(dict(zip(lib["barcode"], lib["strain"])))

    # Older pipeline versions wrote a narrower schema. Normalise before any analysis so
    # that the analyses can assume one shape. See `_add_plate_column`.
    known_plates = set((config.get("plates") or {}).keys())
    for name, frame in (("titers_per_replicate", per_replicate), ("curvefits", curvefits)):
        added = _add_plate_column(frame, known_plates)
        if added:
            notes.append(
                f"`{name}` at {label} has no `plate` column (an older pipeline "
                "version); it was derived from the `replicate` column's plate prefix."
            )

    return Run(
        label=label,
        commit=commit,
        sha=sha,
        subject=subject,
        config=config,
        titers=titers,
        per_replicate=per_replicate,
        curvefits=curvefits,
        plate_to_plate=plate_to_plate,
        sera_by_plate=sera_by_plate,
        plate_qc_drops=plate_qc_drops,
        groups_sera_qc_drops=groups_sera_qc_drops,
        barcode_to_strain=barcode_to_strain,
        notes=notes,
    )


def _add_plate_column(frame: pd.DataFrame, known_plates: set[str]) -> bool:
    """Add a `plate` column derived from `replicate`, in place. Returns whether it did.

    Pipeline versions before the `plate` column existed name each replicate
    ``<plate>-<barcode>`` (or ``<plate>-<strain>`` when barcodes are collapsed). Matching
    against the plate names in the config rather than splitting on the first ``-`` keeps
    this correct for plate names that themselves contain a hyphen, and for strain names
    that contain one (e.g. ``A/StPetersburg/RII-MH144113/2023``).

    ADAPT: if a future pipeline version changes the replicate naming scheme, this is the
    single place to teach it the new one.
    """
    if frame.empty or "plate" in frame.columns or "replicate" not in frame.columns:
        return False

    # Longest match first, so `plate7.5` wins over a hypothetical `plate7`.
    ordered = sorted(known_plates, key=len, reverse=True)

    def derive(replicate: str) -> str | None:
        text = str(replicate)
        for plate in ordered:
            if text.startswith(f"{plate}-"):
                return plate
        return text.split("-", 1)[0] if "-" in text else None

    frame["plate"] = frame["replicate"].map(derive)
    return True


# ======================================================================================
# Small shared helpers
# ======================================================================================


def _log2(series: pd.Series) -> pd.Series:
    """log2 of a titer, with non-positive values as NaN rather than -inf."""
    return numpy.log2(series.where(series > 0))


def uncensored(frame: pd.DataFrame) -> pd.DataFrame:
    """Drop titers sitting at a bound. See the docstring note on censoring."""
    if "titer_bound" not in frame.columns:
        return frame
    return frame[frame["titer_bound"] == UNCENSORED]


def _normalise_to_strain(key: str, run: Run) -> str:
    """Map a QC-drop key onto a strain name.

    Without collapsing, the pipeline keys barcode-level drops by barcode sequence; with
    collapsing the collapsed barcode is named after its strain. Both have to land on the
    strain for the two runs' drop tallies to be comparable at all.
    """
    first = str(key).split()[0] if str(key).split() else str(key)
    if first in run.barcode_to_strain:
        return run.barcode_to_strain[first]
    return first


def per_serum_qc_drops(run: Run) -> pd.DataFrame:
    """Serum-virus titers that existed but failed per-serum QC.

    Preferred source is the `dropped_by_qc` flag on the replicate rows, since that is
    what the pipeline itself acts on. Pipeline versions predating that column need the
    `groups_sera_qc_drops.yml` fallback, without which such a run would appear to drop
    nothing at all -- a silent and badly misleading result.
    """
    if not run.per_replicate.empty and "dropped_by_qc" in run.per_replicate.columns:
        grouped = (
            run.per_replicate.groupby(KEY, as_index=False)["dropped_by_qc"]
            .all()
            .rename(columns={"dropped_by_qc": "all_replicates_dropped"})
        )
        return grouped[grouped["all_replicates_dropped"]][KEY]

    records = [
        {"group": group, "serum": serum, "virus": virus}
        for group, sera in (run.groups_sera_qc_drops or {}).items()
        if isinstance(sera, dict)
        for serum, viruses in sera.items()
        if isinstance(viruses, dict)
        for virus in viruses
    ]
    return pd.DataFrame(records, columns=KEY)


def per_serum_qc_drops_available(run: Run) -> bool:
    """Whether per-serum QC drops can be determined for this run at all.

    A run with neither the `dropped_by_qc` column nor `groups_sera_qc_drops.yml` yields
    an empty drop set, which is indistinguishable from "dropped nothing" unless it is
    reported as unavailable. Since drop counts are the most misread number in this whole
    comparison, a missing value must never be rendered as a zero.
    """
    if not run.per_replicate.empty and "dropped_by_qc" in run.per_replicate.columns:
        return True
    return bool(run.groups_sera_qc_drops)


def plate_titers(run: Run) -> pd.DataFrame:
    """One titer per (group, serum, virus, plate), pre-per-serum-QC.

    Replicates within a plate are combined by median on the log scale, matching how the
    pipeline aggregates replicates. With barcodes collapsed there is normally one
    replicate per plate, so this is a no-op for such runs -- which is exactly the point:
    it puts both runs on the same footing of "one number per plate".
    """
    if run.per_replicate.empty:
        return pd.DataFrame(columns=KEY + ["plate", "titer", "n_replicates_on_plate"])
    frame = uncensored(run.per_replicate).copy()
    frame["log2_titer"] = _log2(frame["titer"])
    frame = frame.dropna(subset=["log2_titer"])
    if "dropped_by_qc" not in frame.columns:
        # An older run without the flag: treat every replicate as retained, so the
        # "as shipped" view degrades to the full set rather than to nothing.
        frame["dropped_by_qc"] = False
    out = frame.groupby(KEY + ["plate"], as_index=False).agg(
        log2_titer=("log2_titer", "median"),
        n_replicates_on_plate=("log2_titer", "size"),
        any_retained=("dropped_by_qc", lambda s: bool((~s.astype(bool)).any())),
    )
    out["titer"] = 2 ** out["log2_titer"]
    return out


# ======================================================================================
# Summary rows: the stable, stackable schema
# ======================================================================================

SUMMARY_COLUMNS = [
    "comparison_id",
    "run_a",
    "run_b",
    "run_a_commit",
    "run_b_commit",
    "section",
    "metric",
    "stratum",
    "value_a",
    "value_b",
    "delta",
    "n",
    "test",
    "statistic",
    "p_value",
    "higher_is_better",
    "note",
]


class Summary:
    """Accumulates the report's headline numbers in one fixed schema.

    The schema is deliberately dataset-agnostic and identical for every comparison, so
    that `summary.csv` files from many comparisons stack into one table. That is what
    makes it possible to ask across datasets whether an effect is consistent, which no
    amount of reading individual HTML reports can answer. `higher_is_better` is carried
    per metric so a meta-analysis can score "which run won" without having to remember
    the polarity of each metric.
    """

    def __init__(self, comparison_id: str, run_a: Run, run_b: Run):
        self.base = {
            "comparison_id": comparison_id,
            "run_a": run_a.label,
            "run_b": run_b.label,
            "run_a_commit": run_a.sha,
            "run_b_commit": run_b.sha,
        }
        self.rows: list[dict] = []

    def add(
        self,
        section: str,
        metric: str,
        *,
        stratum: str = "all",
        value_a=None,
        value_b=None,
        n=None,
        test=None,
        statistic=None,
        p_value=None,
        higher_is_better=None,
        note=None,
    ) -> None:
        delta = None
        if value_a is not None and value_b is not None:
            try:
                delta = float(value_b) - float(value_a)
            except (TypeError, ValueError):
                delta = None
        self.rows.append(
            {
                **self.base,
                "section": section,
                "metric": metric,
                "stratum": stratum,
                "value_a": value_a,
                "value_b": value_b,
                "delta": delta,
                "n": n,
                "test": test,
                "statistic": statistic,
                "p_value": p_value,
                "higher_is_better": higher_is_better,
                "note": note,
            }
        )

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.rows, columns=SUMMARY_COLUMNS)


@dataclass
class Section:
    """One section of the report."""

    key: str
    title: str
    lead: str = ""
    notes: list[str] = field(default_factory=list)
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    charts: list[tuple[str, object]] = field(default_factory=list)
    skipped: str | None = None

    def note(self, text: str) -> None:
        self.notes.append(textwrap.dedent(text).strip())


# ======================================================================================
# Section 1: coverage and concordance of the final titers
# ======================================================================================


def section_coverage(run_a: Run, run_b: Run, summary: Summary) -> Section:
    """Which serum-virus titers exist in which run, and why they are missing.

    This is the first output deliberately: it establishes how much of the two runs is
    even comparable before any quality claim is made, and the reason breakdown is what
    keeps scope differences between the commits from being mistaken for QC behaviour.
    """
    sec = Section(
        "coverage",
        "1a. Coverage: which serum-virus titers exist in which run",
        lead=(
            "Partition of every serum-virus pair either run reports. Missing titers are "
            "attributed to the first applicable cause; only the last three causes say "
            "anything about the approaches being compared."
        ),
    )

    keys_a = set(map(tuple, run_a.titers[KEY].values)) if not run_a.titers.empty else set()
    keys_b = set(map(tuple, run_b.titers[KEY].values)) if not run_b.titers.empty else set()

    plates_a, plates_b = run_a.plates_by_serum, run_b.plates_by_serum
    qc_a = set(map(tuple, per_serum_qc_drops(run_a).values))
    qc_b = set(map(tuple, per_serum_qc_drops(run_b).values))

    # A curve was fit at all iff the pair appears in per_replicate.
    fit_a = (
        set(map(tuple, run_a.per_replicate[KEY].drop_duplicates().values))
        if not run_a.per_replicate.empty
        else set()
    )
    fit_b = (
        set(map(tuple, run_b.per_replicate[KEY].drop_duplicates().values))
        if not run_b.per_replicate.empty
        else set()
    )
    # Curves that were fit on some plate, from the curvefits files: lets a titer lost
    # before any fit be told apart from one whose every fit failed goodness-of-fit.
    curve_a = (
        set(map(tuple, run_a.curvefits[KEY].drop_duplicates().values))
        if not run_a.curvefits.empty and set(KEY) <= set(run_a.curvefits.columns)
        else set()
    )
    curve_b = (
        set(map(tuple, run_b.curvefits[KEY].drop_duplicates().values))
        if not run_b.curvefits.empty and set(KEY) <= set(run_b.curvefits.columns)
        else set()
    )

    def reason(key, run, other_keys, plates_self, plates_other, qc_self, fit_self, curve_self):
        group, serum, virus = key
        if (group, serum) not in {(g, s) for g, s in plates_self}:
            return "serum_absent_from_run"
        if run.strains and virus not in run.strains:
            return "virus_absent_from_run"
        if plates_self.get((group, serum)) != plates_other.get((group, serum)):
            return "plate_set_differs"
        if key in qc_self:
            return "dropped_by_per_serum_qc"
        if key in fit_self:
            return "dropped_by_per_serum_qc"
        if key in curve_self:
            return "lost_at_curvefit_qc"
        return "lost_at_barcode_or_well_qc"

    records = []
    for key in sorted(keys_a | keys_b):
        in_a, in_b = key in keys_a, key in keys_b
        records.append(
            {
                "group": key[0],
                "serum": key[1],
                "virus": key[2],
                "in_a": in_a,
                "in_b": in_b,
                "status": "shared"
                if in_a and in_b
                else (f"only_{run_a.label}" if in_a else f"only_{run_b.label}"),
                "missing_from": None
                if in_a and in_b
                else (run_b.label if in_a else run_a.label),
                "reason_missing": None
                if in_a and in_b
                else (
                    reason(key, run_b, keys_a, plates_b, plates_a, qc_b, fit_b, curve_b)
                    if in_a
                    else reason(
                        key, run_a, keys_b, plates_a, plates_b, qc_a, fit_a, curve_a
                    )
                ),
            }
        )
    coverage = pd.DataFrame(records)
    # Write only the rows that differ. The shared rows are 98% of the frame and carry no
    # information beyond "present in both", which `coverage_overview` already counts and
    # `concordance_per_titer` already enumerates. Among the differing rows `missing_from`
    # determines `in_a`, `in_b` and `status`, so those three are dropped as derived.
    sec.tables["coverage_differences"] = coverage[coverage["status"] != "shared"][
        [*KEY, "missing_from", "reason_missing"]
    ]

    n_union = len(coverage)
    n_shared = int((coverage["status"] == "shared").sum())
    n_only_a = int((coverage["status"] == f"only_{run_a.label}").sum())
    n_only_b = int((coverage["status"] == f"only_{run_b.label}").sum())

    overview = pd.DataFrame(
        [
            {"partition": "shared (both runs)", "n": n_shared},
            {"partition": f"only in {run_a.label}", "n": n_only_a},
            {"partition": f"only in {run_b.label}", "n": n_only_b},
        ]
    ).assign(pct_of_union=lambda x: (100 * x["n"] / max(n_union, 1)).round(2))
    sec.tables["coverage_overview"] = overview

    reasons = (
        coverage.dropna(subset=["reason_missing"])
        .groupby(["missing_from", "reason_missing"], as_index=False)
        .size()
        .rename(columns={"size": "n"})
        .sort_values(["missing_from", "n"], ascending=[True, False])
    )
    sec.tables["coverage_reasons"] = reasons

    # Which strains the unique titers pile up on. Concentration on a few strains is the
    # signal worth chasing: it means one run systematically accepts or rejects a strain.
    by_virus = (
        coverage[coverage["status"] != "shared"]
        .groupby(["missing_from", "virus"], as_index=False)
        .size()
        .rename(columns={"size": "n_titers_missing"})
        .sort_values("n_titers_missing", ascending=False)
    )
    sec.tables["unique_titers_by_virus"] = by_virus

    summary.add("coverage", "n_titers", value_a=len(keys_a), value_b=len(keys_b))
    summary.add("coverage", "n_shared_titers", value_a=n_shared, value_b=n_shared, n=n_union)
    summary.add(
        "coverage",
        "n_titers_unique_to_run",
        value_a=n_only_a,
        value_b=n_only_b,
        n=n_union,
        note="unique titers are not automatically bad; section 3 adjudicates them",
    )
    for (missing_from, why), grp in reasons.groupby(["missing_from", "reason_missing"]):
        n = int(grp["n"].iloc[0])
        summary.add(
            "coverage",
            f"n_missing__{why}",
            stratum=f"missing_from={missing_from}",
            value_a=n if missing_from == run_a.label else None,
            value_b=n if missing_from == run_b.label else None,
        )

    scope = reasons[
        reasons["reason_missing"].isin(
            ["serum_absent_from_run", "virus_absent_from_run", "plate_set_differs"]
        )
    ]["n"].sum()
    if scope:
        sec.note(
            f"""
            {int(scope)} missing titers are scope differences between the commits
            (different sera, strains, or plate sets), not QC behaviour. They are
            excluded from every quality metric. Treat the headline concordance below
            with care: the two runs are not measuring quite the same experiment.
            """
        )
    else:
        sec.note(
            "No scope differences: both runs cover the same sera, strains and plate "
            "sets, so the whole coverage difference is attributable to QC."
        )

    if not by_virus.empty:
        worst = by_virus.iloc[0]
        sec.note(
            f"""
            Unique titers concentrate on particular strains -- the most affected is
            {worst['virus']}, missing {int(worst['n_titers_missing'])} titers from
            {worst['missing_from']}. A strain that one run systematically rejects is
            the most informative thing in this table: section 3 tests whether those
            rejections were justified.
            """
        )

    if not coverage.empty:
        chart_data = overview.rename(columns={"partition": "Partition"})
        sec.charts.append(
            (
                "Coverage partition",
                alt.Chart(chart_data)
                .mark_bar()
                .encode(
                    x=alt.X("n:Q", title="number of serum-virus titers"),
                    y=alt.Y("Partition:N", sort="-x", title=None),
                    color=alt.Color("Partition:N", legend=None),
                    tooltip=["Partition", "n", "pct_of_union"],
                )
                .properties(height=100, width=520),
            )
        )
    return sec


def section_concordance(
    run_a: Run, run_b: Run, matched: pd.DataFrame, summary: Summary
) -> Section:
    """How similar the final titers are, on the matched set."""
    sec = Section(
        "concordance",
        "1b. Concordance: how similar the final titers are",
        lead=(
            "Agreement between the two runs' final titers on the matched set (same "
            "serum, same strain, same plate set), excluding censored titers. This is "
            "descriptive, not a quality metric: it says how much the two approaches "
            "actually differ, and therefore how much any quality difference can matter."
        ),
    )
    if matched.empty:
        sec.skipped = "no matched, uncensored serum-virus pairs in common"
        return sec

    frame = matched.copy()
    frame["log2_fc"] = frame["log2_titer_b"] - frame["log2_titer_a"]
    frame["abs_log2_fc"] = frame["log2_fc"].abs()
    frame["mean_log2_titer"] = (frame["log2_titer_a"] + frame["log2_titer_b"]) / 2
    # Write a narrow, rounded version: these floats come from arithmetic and default to 17
    # significant figures, which costs ~18 characters per value for precision far below
    # anything meaningful. `titer_bound_*` is constant here (the frame is uncensored by
    # construction) and `abs_log2_fc` and `mean_log2_titer` are functions of the columns
    # kept, so all four are dropped from the CSV rather than stored redundantly.
    sec.tables["concordance_per_titer"] = frame[
        [*KEY, "titer_a", "titer_b", "log2_fc"]
    ].assign(log2_fc=lambda x: x["log2_fc"].round(4))

    pearson = frame[["log2_titer_a", "log2_titer_b"]].corr().iloc[0, 1]
    spearman = frame[["log2_titer_a", "log2_titer_b"]].corr(method="spearman").iloc[0, 1]
    median_fc = frame["log2_fc"].median()
    median_abs = frame["abs_log2_fc"].median()

    # A paired test for systematic bias. Wilcoxon rather than a t-test because log
    # fold-changes here have heavy tails.
    bias_test = bias_p = None
    if len(frame) > 10 and frame["log2_fc"].abs().sum() > 0:
        bias_test, bias_p = stats.wilcoxon(frame["log2_fc"], zero_method="zsplit")

    rows = [
        {"metric": "n matched uncensored pairs", "value": len(frame)},
        {"metric": "Pearson r (log2 titer)", "value": round(pearson, 4)},
        {"metric": "Spearman rho", "value": round(spearman, 4)},
        {
            "metric": "median log2 FC (b/a): systematic bias",
            "value": round(median_fc, 4),
        },
        {"metric": "  as fold-change", "value": round(2**median_fc, 4)},
        {"metric": "median |log2 FC|: typical disagreement", "value": round(median_abs, 4)},
        {"metric": "  as fold-change", "value": round(2**median_abs, 4)},
    ]
    for fold in FOLD_CHANGE_BINS:
        frac = float((frame["abs_log2_fc"] > numpy.log2(fold)).mean())
        rows.append({"metric": f"fraction differing > {fold}-fold", "value": round(frac, 4)})
        summary.add(
            "concordance",
            f"frac_differing_gt_{fold}fold",
            value_a=round(frac, 6),
            value_b=round(frac, 6),
            n=len(frame),
            note="symmetric measure of disagreement, not per-run",
        )
    sec.tables["concordance_summary"] = pd.DataFrame(rows)

    summary.add("concordance", "pearson_r_log2", value_a=round(pearson, 6), value_b=round(pearson, 6), n=len(frame))
    summary.add("concordance", "spearman_rho", value_a=round(spearman, 6), value_b=round(spearman, 6), n=len(frame))
    summary.add(
        "concordance",
        "median_log2_fc_b_minus_a",
        value_a=0.0,
        value_b=round(median_fc, 6),
        n=len(frame),
        test="wilcoxon_signed_rank" if bias_p is not None else None,
        statistic=bias_test,
        p_value=bias_p,
        note="positive means run b reports higher titers",
    )
    summary.add(
        "concordance",
        "median_abs_log2_fc",
        value_a=round(median_abs, 6),
        value_b=round(median_abs, 6),
        n=len(frame),
        higher_is_better=False,
        note="symmetric: how much the two runs disagree at all",
    )

    sec.note(
        f"""
        The two runs agree closely for most titers: typical disagreement is
        {2**median_abs:.2f}-fold and
        {100 * (frame['abs_log2_fc'] > 1).mean():.1f}% of titers differ by more than
        2-fold. The decision between the approaches therefore rests on that tail plus
        the coverage difference in 1a, not on the bulk of the data. Sections 2 and 3 are
        run on the tail as well as on everything, since the tail is where the two
        approaches actually make different claims.
        """
    )
    if bias_p is not None and bias_p < 0.05 and abs(median_fc) > 0.05:
        sec.note(
            f"""
            There is a systematic offset: run {run_b.label} reads
            {2**median_fc:.3f}x run {run_a.label} at the median
            (Wilcoxon p={bias_p:.2e}). Check the Bland-Altman panel for whether the
            offset depends on titer magnitude before treating it as a simple scale
            factor.
            """
        )

    per_virus = (
        frame.groupby("virus", as_index=False)
        .agg(
            n=("log2_fc", "size"),
            median_log2_fc=("log2_fc", "median"),
            median_abs_log2_fc=("abs_log2_fc", "median"),
        )
        .sort_values("median_abs_log2_fc", ascending=False)
    )
    sec.tables["concordance_per_virus"] = per_virus

    # Ship only titers and surrogate keys: `log2_fc` and `mean_log2_titer` are functions of
    # `titer_a`/`titer_b`, so `vega` recomputes them client-side rather than costing two
    # full-precision floats on every row, and the serum/strain names come from lookups.
    sample, lookups = _lookup_encode(
        _downsample(frame, 6000)[["group", "serum", "virus", "titer_a", "titer_b"]],
        ["serum", "virus"],
        carry={"group": "serum"},
    )
    tooltip = [
        alt.Tooltip("group:N"),
        alt.Tooltip("serum:N"),
        alt.Tooltip("virus:N"),
        alt.Tooltip("titer_a:Q", title=f"titer, {run_a.label}"),
        alt.Tooltip("titer_b:Q", title=f"titer, {run_b.label}"),
        alt.Tooltip("log2_fc:Q", format=".3f"),
    ]
    concordance_base = _apply_lookups(
        alt.Chart(sample).mark_circle(size=18, opacity=0.4), lookups
    ).transform_calculate(log2_fc="log(datum.titer_b / datum.titer_a) / log(2)")
    sec.charts.append(
        (
            "Titer concordance (log-log). Hover a point to identify the serum and strain.",
            concordance_base.encode(
                x=alt.X("titer_a:Q", scale=alt.Scale(type="log"), title=f"titer, {run_a.label}"),
                y=alt.Y("titer_b:Q", scale=alt.Scale(type="log"), title=f"titer, {run_b.label}"),
                color=alt.Color("group:N"),
                tooltip=tooltip,
            )
            .properties(width=430, height=430)
            .interactive(),
        )
    )
    sec.charts.append(
        (
            "Bland-Altman: does the difference depend on titer magnitude?",
            concordance_base.transform_calculate(
                mean_log2_titer="log(datum.titer_a * datum.titer_b) / (2 * log(2))"
            )
            .encode(
                x=alt.X("mean_log2_titer:Q", title="mean log2 titer of the two runs"),
                y=alt.Y("log2_fc:Q", title=f"log2({run_b.label} / {run_a.label})"),
                color=alt.Color("group:N"),
                tooltip=tooltip,
            )
            .properties(width=430, height=300)
            .interactive(),
        )
    )
    top_virus = per_virus.head(30)
    sec.charts.append(
        (
            "30 strains where the runs disagree most",
            alt.Chart(top_virus)
            .mark_bar()
            .encode(
                x=alt.X("median_abs_log2_fc:Q", title="median |log2 FC|"),
                y=alt.Y("virus:N", sort="-x", title=None),
                tooltip=["virus", "n", "median_log2_fc", "median_abs_log2_fc"],
            )
            .properties(width=430, height=520),
        )
    )
    return sec


def build_matched_set(run_a: Run, run_b: Run) -> tuple[pd.DataFrame, list[str]]:
    """Serum-virus pairs comparable between the runs, with each run's final titer.

    "Comparable" means both runs measured that serum on the same set of plates. See the
    docstring: requiring an identical plate set is what stops a difference in plate
    coverage between the commits from masquerading as a difference in method.
    """
    notes: list[str] = []
    if run_a.titers.empty or run_b.titers.empty:
        return pd.DataFrame(), ["one or both runs have no aggregated titers"]

    cols = KEY + ["titer", "titer_bound", "n_replicates"]
    a = uncensored(run_a.titers)[[c for c in cols if c in run_a.titers.columns]]
    b = uncensored(run_b.titers)[[c for c in cols if c in run_b.titers.columns]]
    merged = a.merge(b, on=KEY, suffixes=("_a", "_b"), validate="one_to_one")

    plates_a, plates_b = run_a.plates_by_serum, run_b.plates_by_serum
    same_plates = merged.apply(
        lambda r: plates_a.get((r["group"], r["serum"])) == plates_b.get((r["group"], r["serum"])),
        axis=1,
    )
    n_dropped = int((~same_plates).sum())
    if n_dropped:
        notes.append(
            f"{n_dropped} shared serum-virus pairs excluded from the matched set "
            "because the serum was measured on different plates in the two runs"
        )
    merged = merged[same_plates].copy()
    merged["log2_titer_a"] = _log2(merged["titer_a"])
    merged["log2_titer_b"] = _log2(merged["titer_b"])
    merged = merged.dropna(subset=["log2_titer_a", "log2_titer_b"])
    return merged, notes


# ======================================================================================
# Section 2: plate-to-plate reproducibility -- the decisive precision metric
# ======================================================================================


def section_reproducibility(
    run_a: Run, run_b: Run, matched: pd.DataFrame, summary: Summary
) -> Section:
    """Agreement between independent remeasurements of the same serum on other plates.

    This is the closest thing to ground truth for precision, and the only metric here
    that does not depend on either run's internal replicate structure. Two statistics
    are reported per serum and plate pair:

      sd_log2_centered  the spread of log2 titer differences after removing the median
                        plate-to-plate offset. This is the precision measure. Removing
                        the offset matters because a plate-wide dilution error is a
                        nuisance shared by every strain on the plate and is not the
                        noise under study.
      pearson_r         retained for continuity with the pipeline's own output, but it
                        depends on the spread of true titers across strains, so a run
                        that slightly compresses dynamic range scores worse without
                        being less precise. Prefer sd_log2_centered.
    """
    sec = Section(
        "reproducibility",
        "2. Plate-to-plate reproducibility (decisive precision metric)",
        lead=(
            "For sera measured on more than one plate, how well do independent "
            "remeasurements agree? Computed on the matched set with each run's "
            "per-serum QC ignored, so neither run is rewarded for filtering harder."
        ),
    )

    pt_a, pt_b = plate_titers(run_a), plate_titers(run_b)
    if pt_a.empty or pt_b.empty:
        sec.skipped = "per-replicate titers unavailable in one or both runs"
        return sec

    repeat_sera = {
        key
        for key, plates in run_a.plates_by_serum.items()
        if len(plates) >= MIN_PLATES_FOR_REPRODUCIBILITY
        and run_b.plates_by_serum.get(key) == plates
    }
    if not repeat_sera:
        sec.skipped = (
            f"no serum is measured on >= {MIN_PLATES_FOR_REPRODUCIBILITY} plates with "
            "the same plate set in both runs, so no independent remeasurement exists"
        )
        sec.note(
            """
            Without repeat sera the comparison cannot measure precision directly and
            must rest on sections 1, 4, 6 and 7. This is the single biggest limitation
            to check when reusing this script on a new dataset.
            """
        )
        return sec

    matched_keys = set(map(tuple, matched[KEY].values)) if not matched.empty else set()
    discordant_keys = (
        set(map(tuple, matched[matched["log2_titer_b"].sub(matched["log2_titer_a"]).abs() > 1][KEY].values))
        if not matched.empty
        else set()
    )

    # The PRIMARY stratum. `matched_keys` comes from the two runs' *final* titers, so it
    # silently conditions on both runs' per-serum QC having passed. That is not neutral:
    # it removes precisely the pairs whose replicates disagreed, which are the pairs a
    # run with more replicates rejects and a run with fewer keeps. Restricting to it
    # therefore biases in favour of whichever run filters more effectively, by deleting
    # the cases where the other run's approach might have been the better one.
    #
    # On this repository's collapse-vs-nocollapse comparison the bias was large enough to
    # reverse the sign of the result: median spread was 0.688 (uncollapsed) vs 0.697
    # (collapsed) on the final-matched set, but 0.750 vs 0.718 -- favouring collapsing --
    # on the pre-QC union. Neither gap was significant, but the metric choice, not the
    # data, decided which approach appeared better. Prefer `pre_qc_union`.
    pre_qc_keys = set(map(tuple, pt_a[KEY].values)) & set(map(tuple, pt_b[KEY].values))
    n_qc_excluded = len(pre_qc_keys - matched_keys)
    winners: dict[str, tuple[float, float, int]] = {}

    def pairwise(pt: pd.DataFrame, restrict: set, label: str) -> pd.DataFrame:
        """Per (serum, plate pair) reproducibility statistics for one run."""
        frame = pt.copy()
        frame["key"] = list(map(tuple, frame[KEY].values))
        frame = frame[frame["key"].isin(restrict)]
        rows = []
        for (group, serum), grp in frame.groupby(["group", "serum"]):
            if (group, serum) not in repeat_sera:
                continue
            plates = sorted(grp["plate"].unique())
            for i, p1 in enumerate(plates):
                for p2 in plates[i + 1 :]:
                    left = grp[grp["plate"] == p1].set_index("virus")["log2_titer"]
                    right = grp[grp["plate"] == p2].set_index("virus")["log2_titer"]
                    common = left.index.intersection(right.index)
                    if len(common) < 3:
                        continue
                    diff = (right[common] - left[common]).astype(float)
                    offset = float(diff.median())
                    centered = diff - offset
                    r = (
                        float(numpy.corrcoef(left[common], right[common])[0, 1])
                        if len(common) > 2 and left[common].std() > 0 and right[common].std() > 0
                        else numpy.nan
                    )
                    rows.append(
                        {
                            "run": label,
                            "group": group,
                            "serum": serum,
                            "plate_1": p1,
                            "plate_2": p2,
                            "n_strains": int(len(common)),
                            "median_log2_offset": round(offset, 4),
                            "sd_log2_centered": round(float(centered.std(ddof=1)), 4),
                            "rmsd_log2": round(float(numpy.sqrt((diff**2).mean())), 4),
                            "pearson_r": round(r, 4) if r == r else numpy.nan,
                        }
                    )
        return pd.DataFrame(rows)

    for scope_label, restrict in [
        ("pre_qc_union", pre_qc_keys),
        ("final_matched", matched_keys),
        ("discordant_tail", discordant_keys),
    ]:
        if not restrict:
            continue
        rep_a = pairwise(pt_a, restrict, run_a.label)
        rep_b = pairwise(pt_b, restrict, run_b.label)
        if rep_a.empty or rep_b.empty:
            continue
        both = pd.concat([rep_a, rep_b], ignore_index=True)
        sec.tables[f"reproducibility_{scope_label}"] = both

        pair_cols = ["group", "serum", "plate_1", "plate_2"]
        paired = rep_a.merge(rep_b, on=pair_cols, suffixes=("_a", "_b"))
        sec.tables[f"reproducibility_{scope_label}_paired"] = paired

        for metric, higher_better in [
            ("sd_log2_centered", False),
            ("rmsd_log2", False),
            ("pearson_r", True),
        ]:
            va = paired[f"{metric}_a"].astype(float)
            vb = paired[f"{metric}_b"].astype(float)
            ok = va.notna() & vb.notna()
            if ok.sum() < 3:
                continue
            test = p = None
            if ok.sum() >= 6 and float((vb[ok] - va[ok]).abs().sum()) > 0:
                test, p = stats.wilcoxon(va[ok], vb[ok], zero_method="zsplit")
            underpowered = int(ok.sum()) < MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM
            note = "median over serum x plate-pair; paired test across those pairs"
            if underpowered:
                note = (
                    f"UNDERPOWERED: only {int(ok.sum())} serum x plate-pair comparisons "
                    f"(< {MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM}); do not interpret"
                )
            summary.add(
                "reproducibility",
                metric,
                stratum=scope_label
                + ("__UNDERPOWERED" if underpowered else ""),
                value_a=round(float(va[ok].median()), 6),
                value_b=round(float(vb[ok].median()), 6),
                n=int(ok.sum()),
                test="wilcoxon_signed_rank" if p is not None else None,
                statistic=test,
                p_value=p,
                higher_is_better=higher_better,
                note=note,
            )
            # Tail as well as centre. Two runs can have the same median spread while one
            # has appreciably worse worst cases, and for titers it is the worst cases
            # that mislead a downstream analysis. On this repository's comparison the
            # median favoured keeping barcodes separate while the 90th percentile
            # favoured collapsing -- a disagreement invisible from medians alone.
            if metric == "sd_log2_centered" and not underpowered:
                for quantile in (0.75, 0.90):
                    summary.add(
                        "reproducibility",
                        f"sd_log2_centered_p{int(quantile * 100)}",
                        stratum=scope_label,
                        value_a=round(float(va[ok].quantile(quantile)), 6),
                        value_b=round(float(vb[ok].quantile(quantile)), 6),
                        n=int(ok.sum()),
                        higher_is_better=False,
                        note="tail of the spread distribution: worst-case reproducibility",
                    )
                summary.add(
                    "reproducibility",
                    "n_pairs_where_b_more_reproducible",
                    stratum=scope_label,
                    value_a=int((va[ok] < vb[ok]).sum()),
                    value_b=int((vb[ok] < va[ok]).sum()),
                    n=int(ok.sum()),
                    note="sign test view: how often each run wins, irrespective of size",
                )
        if len(paired) < MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM:
            sec.note(
                f"""
                The `{scope_label}` stratum reduces to only {len(paired)} serum x
                plate-pair comparisons, below the {MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM}
                needed for a meaningful claim. Its numbers are in the tables and are
                flagged `__UNDERPOWERED` in summary.csv, but should not be interpreted:
                with this few comparisons the apparent difference is driven by which
                sera happened to survive the restriction.
                """
            )

        med_a = float(paired["sd_log2_centered_a"].median())
        med_b = float(paired["sd_log2_centered_b"].median())
        winners[scope_label] = (med_a, med_b, len(paired))
        if scope_label == "pre_qc_union":
            better = run_b.label if med_b < med_a else run_a.label
            p90_a = float(paired["sd_log2_centered_a"].quantile(0.90))
            p90_b = float(paired["sd_log2_centered_b"].quantile(0.90))
            better_tail = run_b.label if p90_b < p90_a else run_a.label
            sec.note(
                f"""
                On the pre-QC union (the primary stratum), median plate-to-plate spread is
                {med_a:.3f} log2 ({2**med_a:.2f}-fold) for {run_a.label} and
                {med_b:.3f} ({2**med_b:.2f}-fold) for {run_b.label} over {len(paired)}
                serum x plate-pair comparisons, so **{better}** is more reproducible at
                the centre. At the 90th percentile the figures are {p90_a:.3f} vs
                {p90_b:.3f}, so **{better_tail}** has the better worst cases.
                """
                + (
                    " Centre and tail agree, which is the straightforward case."
                    if better == better_tail
                    else " **Centre and tail disagree**, so one approach is better "
                    "typically and "
                    "the other is better in the tail. Which matters more depends on "
                    "whether a few badly wrong titers or a slightly worse average is "
                    "the bigger problem downstream -- this is a judgement call, not "
                    "something the data settles."
                )
            )
            sec.note(
                f"""
                Based on {len(repeat_sera)} repeat sera. Every serum-virus titer outside
                those sera has no independent remeasurement, so this section speaks for
                only part of the dataset and cannot be broken down per strain with
                confidence. Read the paired p-value in summary.csv before treating any
                small gap as real: with this many comparisons only a large effect is
                detectable, so a non-significant result means "not shown", not "no
                difference".
                """
            )

    # Whether the choice of stratum changed the answer. If it did, the metric definition
    # rather than the data is deciding the comparison, and that has to be said out loud.
    if "pre_qc_union" in winners and "final_matched" in winners:
        (pre_a, pre_b, n_pre) = winners["pre_qc_union"]
        (fin_a, fin_b, n_fin) = winners["final_matched"]
        if (pre_b < pre_a) != (fin_b < fin_a):
            sec.note(
                f"""
                **The two strata disagree, so the choice of stratum decides the answer.**
                Pre-QC union: {pre_a:.3f} vs {pre_b:.3f} favours
                {run_b.label if pre_b < pre_a else run_a.label}. Final-matched:
                {fin_a:.3f} vs {fin_b:.3f} favours
                {run_b.label if fin_b < fin_a else run_a.label}. `final_matched`
                conditions on both runs' per-serum QC passing, which removes
                {n_qc_excluded} pairs (per-serum QC drops, plus any censored in either
                run) -- among them exactly the pairs whose replicates disagreed, which one
                run rejects and the other keeps. That makes it the biased view and
                `pre_qc_union` the one to trust. Conclude that neither approach is
                clearly better on precision rather than quoting whichever stratum
                supports a preferred answer.
                """
            )
        else:
            sec.note(
                f"""
                Both strata point the same way ({n_pre} and {n_fin} comparisons), so the
                result does not depend on whether QC-dropped pairs are included. That is
                a meaningfully stronger finding than either stratum alone.
                """
            )
    if n_qc_excluded:
        sec.note(
            f"""
            `pre_qc_union` includes {n_qc_excluded} serum-virus pairs that `final_matched`
            excludes because one run's per-serum QC dropped them or one run censored them.
            Comparing the two
            answers a distinct question: whether the pairs a run discards are ones the
            other run measures well. If a run's spread barely worsens when they are added
            back, its QC was discarding data it could actually handle.
            """
        )

    # As-shipped variant: each run's own retained titers, on its own set. Reported for
    # completeness and explicitly not comparable -- the sets differ by construction.
    shipped = []
    for run, pt in [(run_a, pt_a), (run_b, pt_b)]:
        if "any_retained" not in pt.columns:
            continue
        kept = pt[pt["any_retained"].astype(bool)]
        keys = set(map(tuple, kept[KEY].values))
        shipped.append(pairwise(kept, keys, run.label))
    if shipped:
        shipped_frame = pd.concat([s for s in shipped if not s.empty], ignore_index=True)
        if not shipped_frame.empty:
            sec.tables["reproducibility_as_shipped"] = shipped_frame
            agg = shipped_frame.groupby("run", as_index=False).agg(
                n_pairs=("sd_log2_centered", "size"),
                median_sd_log2_centered=("sd_log2_centered", "median"),
                median_pearson_r=("pearson_r", "median"),
                mean_n_strains=("n_strains", "mean"),
            )
            sec.tables["reproducibility_as_shipped_summary"] = agg
            sec.note(
                """
                The "as shipped" table uses each run's own retained titers and is NOT a
                fair comparison: a run that drops more discordant titers scores better
                on the ones it keeps. Compare `mean_n_strains` between the runs to see
                how different the underlying sets are. It is included only because it
                is what the pipeline's own plate_to_plate_corr output reflects.
                """
            )

    # Chart the primary stratum, not the biased one.
    key_table = sec.tables.get("reproducibility_pre_qc_union")
    if key_table is not None and not key_table.empty:
        sec.charts.append(
            (
                "Plate-to-plate spread per serum and plate pair (lower is more precise)",
                alt.Chart(key_table.assign(pair=lambda x: x["plate_1"] + " vs " + x["plate_2"]))
                .mark_point(size=55, filled=True, opacity=0.75)
                .encode(
                    x=alt.X("sd_log2_centered:Q", title="SD of log2 titer difference (offset removed)"),
                    y=alt.Y("serum:N", title=None),
                    color=alt.Color("run:N", title="run"),
                    shape=alt.Shape("pair:N", title="plate pair"),
                    tooltip=["run", "serum", "pair", "n_strains", "sd_log2_centered", "rmsd_log2", "pearson_r", "median_log2_offset"],
                )
                .properties(width=460, height=alt.Step(16)),
            )
        )
        paired = sec.tables.get("reproducibility_pre_qc_union_paired")
        if paired is not None and not paired.empty:
            lim = float(
                pd.concat([paired["sd_log2_centered_a"], paired["sd_log2_centered_b"]]).max()
            )
            scatter = (
                alt.Chart(paired.assign(pair=lambda x: x["plate_1"] + " vs " + x["plate_2"]))
                .mark_circle(size=70, opacity=0.75)
                .encode(
                    x=alt.X("sd_log2_centered_a:Q", scale=alt.Scale(domain=[0, lim * 1.05]), title=f"spread, {run_a.label}"),
                    y=alt.Y("sd_log2_centered_b:Q", scale=alt.Scale(domain=[0, lim * 1.05]), title=f"spread, {run_b.label}"),
                    # The merge suffixes every column shared by the two runs, so the
                    # strain count is n_strains_a / n_strains_b here, not n_strains.
                    tooltip=[
                        "serum",
                        "pair",
                        "n_strains_a",
                        "n_strains_b",
                        "sd_log2_centered_a",
                        "sd_log2_centered_b",
                    ],
                )
                .properties(width=380, height=380)
            )
            diag = (
                alt.Chart(pd.DataFrame({"x": [0, lim * 1.05]}))
                .mark_line(strokeDash=[4, 4], color="grey")
                .encode(x="x:Q", y="x:Q")
            )
            sec.charts.append(
                (
                    "Paired: points below the diagonal favour "
                    f"{run_b.label}, above favour {run_a.label}",
                    scatter + diag,
                )
            )
    return sec


# ======================================================================================
# Section 3: are the run-unique titers any good?
# ======================================================================================


def section_adjudicate(run_a: Run, run_b: Run, summary: Summary) -> Section:
    """Test whether titers only one run reports actually replicate across plates.

    This converts the confounded "who drops fewer titers" question into a real one. If
    the titers unique to run X replicate badly across plates, the other run was right
    to discard them and X is silently keeping bad data. If they replicate as well as
    shared titers, the other run is over-filtering.
    """
    sec = Section(
        "adjudication",
        "3. Are the run-unique titers any good?",
        lead=(
            "For serum-virus titers that only one run reports, how consistent are that "
            "run's own per-plate measurements? Compared against the shared titers in "
            "the same run as a baseline."
        ),
    )

    keys_a = set(map(tuple, run_a.titers[KEY].values)) if not run_a.titers.empty else set()
    keys_b = set(map(tuple, run_b.titers[KEY].values)) if not run_b.titers.empty else set()
    shared = keys_a & keys_b

    rows = []
    for run, own, other in [(run_a, keys_a, keys_b), (run_b, keys_b, keys_a)]:
        pt = plate_titers(run)
        if pt.empty:
            continue
        multi = pt.groupby(KEY, as_index=False).agg(
            n_plates=("plate", "nunique"), sd_log2=("log2_titer", "std")
        )
        multi = multi[
            (multi["n_plates"] >= MIN_PLATES_FOR_REPRODUCIBILITY) & multi["sd_log2"].notna()
        ]
        if multi.empty:
            continue
        multi["key"] = list(map(tuple, multi[KEY].values))
        # Only titers this run actually reports can be "unique to" it. `plate_titers` is
        # pre-per-serum-QC and so also contains pairs this run itself dropped; leaving
        # those in would label them unique to this run when neither run reports them.
        multi = multi[multi["key"].isin(own)]
        if multi.empty:
            continue
        multi["status"] = numpy.where(
            multi["key"].isin(other), "shared", f"unique_to_{run.label}"
        )
        multi["run"] = run.label
        rows.append(multi.drop(columns=["key"]))

    if not rows:
        sec.skipped = (
            "no serum-virus titer has measurements on multiple plates, so consistency "
            "of the unique titers cannot be assessed"
        )
        return sec

    consistency = pd.concat(rows, ignore_index=True)
    # `sd_log2` is arithmetic output and defaults to 17 significant figures; `group` is
    # determined by `serum`.
    sec.tables["unique_titer_consistency"] = consistency.drop(columns=["group"]).assign(
        sd_log2=lambda x: x["sd_log2"].round(4)
    )

    agg = (
        consistency.groupby(["run", "status"], as_index=False)
        .agg(n=("sd_log2", "size"), median_sd_log2=("sd_log2", "median"))
        .sort_values(["run", "status"])
    )
    sec.tables["unique_titer_consistency_summary"] = agg
    verdicts = []
    for run in (run_a, run_b):
        sub = agg[agg["run"] == run.label]
        uniq = sub[sub["status"] == f"unique_to_{run.label}"]
        shr = sub[sub["status"] == "shared"]
        if uniq.empty or shr.empty:
            continue
        u, s = float(uniq["median_sd_log2"].iloc[0]), float(shr["median_sd_log2"].iloc[0])
        n_u = int(uniq["n"].iloc[0])
        # Unique titers often belong to single-plate sera, leaving very few with the
        # multi-plate data this test needs. Report the number rather than a verdict.
        underpowered = n_u < MIN_PAIRS_FOR_REPRODUCIBILITY_CLAIM
        summary.add(
            "adjudication",
            "median_sd_log2_unique_vs_shared",
            stratum=f"run={run.label}" + ("__UNDERPOWERED" if underpowered else ""),
            value_a=round(s, 6),
            value_b=round(u, 6),
            n=n_u,
            higher_is_better=False,
            note=(
                f"UNDERPOWERED: only {n_u} unique titers have multi-plate data"
                if underpowered
                else "value_a is this run's shared-titer baseline, value_b its unique "
                "titers; unique much worse means the other run was right to drop them"
            ),
        )
        ratio = u / s if s else float("nan")
        if underpowered:
            verdicts.append(
                f"For **{run.label}**, only {n_u} of its unique titers have multi-plate "
                f"data, too few to adjudicate (they scatter {u:.3f} log2 vs {s:.3f} for "
                "its shared titers, but do not read anything into that)."
            )
            continue
        verdicts.append(
            f"For **{run.label}**, the {n_u} unique titers with multi-plate data scatter "
            f"{u:.3f} log2 vs {s:.3f} for its shared titers ({ratio:.2f}x). "
            + (
                "They are appreciably noisier, so the other run was plausibly right to "
                "discard them."
                if ratio > 1.3
                else "They are about as consistent as the rest, so the other run may be "
                "over-filtering."
            )
        )
    for verdict in verdicts:
        sec.note(verdict)
    if not verdicts:
        sec.note(
            """
            Not enough unique titers have multi-plate measurements to adjudicate them.
            This is common when the unique titers belong to single-plate sera; it is a
            real limit, not a null result.
            """
        )

    # A boxplot aggregates, so it needs only the value and the two grouping columns --
    # shipping serum/virus/group per row would cost far more than the plot uses.
    sec.charts.append(
        (
            "Per-plate scatter of unique vs shared titers (lower is more consistent)",
            alt.Chart(
                consistency[["run", "status", "sd_log2"]].assign(
                    sd_log2=lambda x: x["sd_log2"].round(4)
                )
            )
            .mark_boxplot(extent="min-max")
            .encode(
                x=alt.X("sd_log2:Q", title="SD of log2 titer across plates"),
                y=alt.Y("status:N", title=None),
                color=alt.Color("status:N", legend=None),
                row=alt.Row("run:N", title="run"),
            )
            .properties(width=430, height=90),
        )
    )
    return sec


# ======================================================================================
# Section 4: mechanism -- where does the noise live, and do the curves fit better?
# ======================================================================================


def section_variance_components(run_a: Run, run_b: Run, summary: Summary) -> Section:
    """Split log2 titer noise into within-plate (between replicate) and between-plate.

    This is the analysis that explains *why* one approach wins. If replicates within a
    plate scatter about as much as plates do, then a strain's barcodes are genuine
    replicates and averaging them buys real precision. If within-plate scatter is small
    relative to between-plate, the barcodes are pseudo-replicates sharing plate-level
    error, the effective replicate count is ~1, and collapsing them costs little.
    """
    sec = Section(
        "variance",
        "4a. Variance components: are replicates real or pseudo-replicates?",
        lead=(
            "Within-plate scatter is between replicates of the same strain on one "
            "plate (for an uncollapsed run, its barcodes). Between-plate scatter is "
            "between independent remeasurements. Their ratio says how much averaging "
            "replicates can possibly help."
        ),
    )
    rows = []
    for run in (run_a, run_b):
        if run.per_replicate.empty:
            continue
        frame = uncensored(run.per_replicate).copy()
        frame["log2_titer"] = _log2(frame["titer"])
        frame = frame.dropna(subset=["log2_titer"])
        if frame.empty:
            continue
        within = (
            frame.groupby(KEY + ["plate"], as_index=False)
            .agg(n=("log2_titer", "size"), sd=("log2_titer", "std"))
            .dropna(subset=["sd"])
        )
        plate_means = frame.groupby(KEY + ["plate"], as_index=False)["log2_titer"].median()
        between = (
            plate_means.groupby(KEY, as_index=False)
            .agg(n_plates=("plate", "nunique"), sd=("log2_titer", "std"))
            .dropna(subset=["sd"])
        )
        rows.append(
            {
                "run": run.label,
                "collapsed": run.collapsed,
                "n_within_plate_estimates": len(within),
                "median_within_plate_sd_log2": round(float(within["sd"].median()), 4)
                if not within.empty
                else numpy.nan,
                "n_between_plate_estimates": len(between),
                "median_between_plate_sd_log2": round(float(between["sd"].median()), 4)
                if not between.empty
                else numpy.nan,
                "mean_replicates_per_plate": round(
                    float(frame.groupby(KEY + ["plate"]).size().mean()), 3
                ),
            }
        )
    if not rows:
        sec.skipped = "no per-replicate titers available"
        return sec

    table = pd.DataFrame(rows)
    sec.tables["variance_components"] = table
    for row in table.itertuples():
        summary.add(
            "variance",
            "median_within_plate_sd_log2",
            stratum=f"run={row.run}",
            value_a=row.median_within_plate_sd_log2 if row.run == run_a.label else None,
            value_b=row.median_within_plate_sd_log2 if row.run == run_b.label else None,
            n=row.n_within_plate_estimates,
            note="NaN/absent when a run has one replicate per plate, e.g. collapsed",
        )
        summary.add(
            "variance",
            "median_between_plate_sd_log2",
            stratum=f"run={row.run}",
            value_a=row.median_between_plate_sd_log2 if row.run == run_a.label else None,
            value_b=row.median_between_plate_sd_log2 if row.run == run_b.label else None,
            n=row.n_between_plate_estimates,
            higher_is_better=False,
        )

    # Turn the two scatter estimates into the quantity that actually matters: how much
    # precision averaging k replicates can buy. Between-plate scatter is measured from
    # plate medians, so it is already net of replicate averaging and behaves as the
    # irreducible plate-level term; the within-plate term is the part averaging shrinks.
    uncollapsed = table[~table["collapsed"].astype(bool)]
    if not uncollapsed.empty:
        row = uncollapsed.iloc[0]
        w, b, k = (
            row["median_within_plate_sd_log2"],
            row["median_between_plate_sd_log2"],
            row["mean_replicates_per_plate"],
        )
        if w == w and b == b and b > 0 and k and k > 1:
            plate_var, within_var = float(b) ** 2, float(w) ** 2
            sd_single = float(numpy.sqrt(plate_var + within_var))
            sd_averaged = float(numpy.sqrt(plate_var + within_var / float(k)))
            gain = 100 * (1 - sd_averaged / sd_single)
            sec.tables["replicate_averaging_gain"] = pd.DataFrame(
                [
                    {
                        "run": row["run"],
                        "within_plate_sd_log2": round(float(w), 4),
                        "between_plate_sd_log2": round(float(b), 4),
                        "mean_replicates_per_plate": round(float(k), 3),
                        "predicted_sd_single_replicate": round(sd_single, 4),
                        "predicted_sd_averaging_k_replicates": round(sd_averaged, 4),
                        "predicted_pct_sd_reduction_from_averaging": round(gain, 2),
                    }
                ]
            )
            summary.add(
                "variance",
                "predicted_pct_sd_reduction_from_replicate_averaging",
                stratum=f"run={row['run']}",
                value_a=round(gain, 4),
                value_b=None,
                n=int(row["n_within_plate_estimates"]),
                note=(
                    "upper bound on what keeping replicates separate can buy; compare "
                    "against the observed reproducibility gap in section 2"
                ),
            )
            sec.note(
                f"""
                In {row['run']} (replicates kept separate), replicates of a strain on one
                plate scatter {w:.3f} log2 ({2**float(w):.2f}-fold) while plates scatter
                {b:.3f} ({2**float(b):.2f}-fold). Because the plate-level term is much the
                larger and is *not* reduced by averaging, averaging the
                {k:.1f} replicates per plate is predicted to cut total scatter only from
                {sd_single:.3f} to {sd_averaged:.3f} log2, a **{gain:.1f}% reduction**.
                """
                + (
                    "That is a small headroom, so losing replicate averaging costs "
                    "little and can plausibly be offset by any reduction in counting "
                    "noise from summing counts. Compare this predicted gain against the "
                    "observed reproducibility gap in section 2: if they are the same "
                    "size, replicate averaging explains the whole difference."
                    if gain < 15
                    else "That is substantial headroom, so keeping replicates separate "
                    "buys real precision that collapsing gives up. Section 2 should show "
                    "a gap of roughly this size if that is the dominant effect."
                )
            )
    sec.note(
        """
        A run with one replicate per plate has no within-plate estimate at all, which is
        not a sign of low noise -- it is the absence of a measurement. That is also why
        such a run's `titer_sem` is mostly empty and why its `max_fold_change_from_median`
        QC check cannot fire.
        """
    )
    return sec


def section_curvefit_quality(run_a: Run, run_b: Run, summary: Summary) -> Section:
    """Compare goodness of fit of the individual neutralization curves.

    Direct evidence on the mechanism: if summing counts before fitting reduces noise,
    the fits should be measurably better. Unlike the titer comparisons this needs no
    matching, because every curve is its own observation -- but note the two runs fit
    *different numbers* of curves, so only the distributions are comparable, not counts.
    """
    sec = Section(
        "curvefits",
        "4b. Curve-fit quality",
        lead=(
            "Distribution of R^2 and RMSD over all fitted curves. Better fits mean less "
            "noise in the underlying fraction-infectivity measurements."
        ),
    )
    rows, per_curve = [], []
    for run in (run_a, run_b):
        cf = run.curvefits
        if cf.empty or not {"r2", "rmsd"} <= set(cf.columns):
            continue
        rows.append(
            {
                "run": run.label,
                "n_curves": len(cf),
                "mean_r2": round(float(cf["r2"].mean()), 4),
                "median_r2": round(float(cf["r2"].median()), 4),
                "frac_r2_below_0.9": round(float((cf["r2"] < 0.9).mean()), 4),
                "min_r2": round(float(cf["r2"].min()), 4),
                "mean_rmsd": round(float(cf["rmsd"].mean()), 4),
                "median_rmsd": round(float(cf["rmsd"].median()), 4),
                "max_rmsd": round(float(cf["rmsd"].max()), 4),
            }
        )
        per_curve.append(cf.assign(run=run.label)[["run", "plate", "serum", "virus", "r2", "rmsd"]])
    if not rows:
        sec.skipped = "no curvefits.csv with r2/rmsd columns in either run"
        return sec

    table = pd.DataFrame(rows)
    sec.tables["curvefit_quality"] = table
    a = table[table["run"] == run_a.label]
    b = table[table["run"] == run_b.label]
    if not a.empty and not b.empty:
        for metric, higher_better in [
            ("median_r2", True),
            ("mean_r2", True),
            ("frac_r2_below_0.9", False),
            ("median_rmsd", False),
        ]:
            summary.add(
                "curvefits",
                metric,
                value_a=float(a[metric].iloc[0]),
                value_b=float(b[metric].iloc[0]),
                higher_is_better=higher_better,
                note="distributions over all fitted curves; curve counts differ by design",
            )
        summary.add(
            "curvefits",
            "n_curves",
            value_a=int(a["n_curves"].iloc[0]),
            value_b=int(b["n_curves"].iloc[0]),
            note="not a quality metric: collapsing fits one curve per strain, not per barcode",
        )
        better = run_b.label if float(b["median_r2"].iloc[0]) > float(a["median_r2"].iloc[0]) else run_a.label
        sec.note(
            f"""
            Median R^2 is {float(a['median_r2'].iloc[0]):.4f} for {run_a.label} and
            {float(b['median_r2'].iloc[0]):.4f} for {run_b.label}; median RMSD
            {float(a['median_rmsd'].iloc[0]):.4f} vs {float(b['median_rmsd'].iloc[0]):.4f}.
            **{better}** fits the better curves. This is evidence about the mechanism,
            not about the final titers: better individual fits can still be outweighed
            by the loss of replicate averaging and of replicate-based error detection.
            """
        )

    curves = pd.concat(per_curve, ignore_index=True)
    sec.tables["curvefit_per_curve"] = curves
    # These two charts only bin and count, so bin in `pandas` and ship the ~60 bins per run
    # instead of tens of thousands of raw values. That is both far smaller and strictly more
    # faithful than plotting a downsample, since every curve contributes to the counts.
    for column, title in [("r2", "R^2"), ("rmsd", "RMSD")]:
        sec.charts.append(
            (
                f"Distribution of curve-fit {title} (all {len(curves):,} curves)",
                alt.Chart(_histogram(curves, column, "run"))
                .mark_area(opacity=0.5, interpolate="step")
                .encode(
                    x=alt.X("bin_start:Q", title=title),
                    x2=alt.X2("bin_end:Q"),
                    y=alt.Y("count:Q", stack=None, title="curves"),
                    color=alt.Color("run:N"),
                    tooltip=["run:N", "bin_start:Q", "bin_end:Q", "count:Q"],
                )
                .properties(width=430, height=230),
            )
        )
    return sec


# ======================================================================================
# Section 5: retention accounting by stage
# ======================================================================================


def section_retention(run_a: Run, run_b: Run, summary: Summary) -> Section:
    """Where each run loses data, split by QC stage.

    Reported by stage rather than as one total because the stages mean different things,
    and because the per-serum stage is the one whose comparison is confounded.
    """
    sec = Section(
        "retention",
        "5. Retention accounting by QC stage",
        lead=(
            "Drops at each QC stage, with barcode-level keys normalised to strains so "
            "the two runs are counted in the same units."
        ),
    )
    rows = []
    for run in (run_a, run_b):
        counts = {"run": run.label, "collapsed": run.collapsed}
        barcode_stage = curvefit_stage = 0
        for plate, drops in (run.plate_qc_drops or {}).items():
            if not isinstance(drops, dict):
                continue
            for stage_key, entries in drops.items():
                if not isinstance(entries, dict):
                    continue
                if stage_key in BARCODE_STAGE_DROP_KEYS:
                    barcode_stage += len(entries)
                elif stage_key in CURVEFIT_STAGE_DROP_KEYS:
                    curvefit_stage += len(entries)
        counts["barcode_or_well_stage_drops"] = (
            barcode_stage if run.plate_qc_drops else numpy.nan
        )
        counts["curvefit_stage_drops"] = curvefit_stage if run.plate_qc_drops else numpy.nan
        available = per_serum_qc_drops_available(run)
        counts["per_serum_qc_titer_drops"] = (
            len(per_serum_qc_drops(run)) if available else numpy.nan
        )
        if not available:
            sec.note(
                f"""
                Per-serum QC drops cannot be determined for **{run.label}**: it has
                neither a `dropped_by_qc` column nor
                `results/qc_drops/groups_sera_qc_drops.yml`. Reported as blank rather
                than zero, since zero would read as "dropped nothing".
                """
            )
        if not run.plate_qc_drops:
            sec.note(
                f"""
                Plate-stage QC drops cannot be determined for **{run.label}**:
                `results/qc_drops/plate_qc_drops.yml` is not committed at that commit.
                Reported as blank rather than zero.
                """
            )
        counts["n_final_titers"] = len(run.titers)
        if not run.titers.empty and "n_replicates" in run.titers.columns:
            counts["frac_titers_with_1_replicate"] = round(
                float((run.titers["n_replicates"] <= 1).mean()), 4
            )
            counts["median_n_replicates"] = float(run.titers["n_replicates"].median())
        rows.append(counts)

    table = pd.DataFrame(rows)
    sec.tables["retention_by_stage"] = table
    a = table[table["run"] == run_a.label]
    b = table[table["run"] == run_b.label]
    if not a.empty and not b.empty:
        for metric in [
            "barcode_or_well_stage_drops",
            "curvefit_stage_drops",
            "per_serum_qc_titer_drops",
            "n_final_titers",
        ]:
            summary.add(
                "retention",
                metric,
                value_a=a[metric].iloc[0],
                value_b=b[metric].iloc[0],
                note=(
                    "CONFOUNDED: see the note in this section"
                    if metric == "per_serum_qc_titer_drops"
                    else None
                ),
            )
        if "frac_titers_with_1_replicate" in table.columns:
            fa = float(a["frac_titers_with_1_replicate"].iloc[0])
            fb = float(b["frac_titers_with_1_replicate"].iloc[0])
            summary.add(
                "retention",
                "frac_titers_with_1_replicate",
                value_a=fa,
                value_b=fb,
                higher_is_better=False,
                note="a titer with one replicate cannot fail max_fold_change_from_median",
            )
            sec.note(
                f"""
                **The per-serum drop counts are not a quality signal on their own.**
                `max_fold_change_from_median` compares replicates to their median, so a
                titer with a single replicate can never fail it. Single-replicate
                titers: {100 * fa:.1f}% in {run_a.label}, {100 * fb:.1f}% in
                {run_b.label}. Whichever run has more of them will appear to drop fewer
                titers regardless of data quality. Section 3 is the substitute that
                actually tests whether the retained titers are sound.
                """
            )
    sec.note(
        """
        A run that fits one curve per strain has no redundancy against curve-fit
        failure: a single failed fit removes the titer outright, whereas a run with
        several replicates per strain still reports a titer from the survivors. Check
        `lost_at_curvefit_qc` in the section 1a reason table for the size of this
        effect -- it is a robustness advantage independent of precision.
        """
    )
    return sec


# ======================================================================================
# Section 6: titer-matrix coherence, a replicate-free noise proxy over all sera
# ======================================================================================


def section_coherence(run_a: Run, run_b: Run, matched: pd.DataFrame, summary: Summary) -> Section:
    """Variance explained by the leading principal components of the titer matrix.

    Antigenic structure is genuinely low-rank: sera respond to strains in a few
    correlated patterns. Cleaner data concentrates more variance in the leading
    components. Unlike the reproducibility sections this uses every serum, not just the
    repeat ones, so it is the main quality signal available for single-plate sera.

    ASSUMPTION: both runs are compared on the identical submatrix, so the comparison is
    not affected by one run simply having more entries.
    """
    sec = Section(
        "coherence",
        "6. Titer-matrix coherence (all sera, no replicates needed)",
        lead=(
            "Fraction of variance captured by the leading principal components of the "
            "serum-by-strain log2 titer matrix, on the identical complete submatrix for "
            "both runs. More variance in the early components means less noise."
        ),
    )
    if matched.empty:
        sec.skipped = "no matched titers to build a matrix from"
        return sec

    wide_a = matched.pivot_table(index=["group", "serum"], columns="virus", values="log2_titer_a")
    wide_b = matched.pivot_table(index=["group", "serum"], columns="virus", values="log2_titer_b")
    # Restrict to a complete submatrix common to both runs: drop the strains, then the
    # sera, with any missing entry. Reported dimensions make the cost of this explicit.
    both_ok = wide_a.notna() & wide_b.notna()
    keep_virus = both_ok.columns[both_ok.all(axis=0)]
    if len(keep_virus) < 3:
        keep_virus = both_ok.columns[both_ok.mean(axis=0) > 0.95]
        both_ok = both_ok[keep_virus]
        keep_serum = both_ok.index[both_ok.all(axis=1)]
    else:
        keep_serum = both_ok.index
    wide_a = wide_a.loc[keep_serum, keep_virus].dropna(axis=0, how="any")
    wide_b = wide_b.loc[wide_a.index, keep_virus]
    if wide_a.shape[0] < 4 or wide_a.shape[1] < 4:
        sec.skipped = (
            f"complete submatrix too small ({wide_a.shape[0]} sera x "
            f"{wide_a.shape[1]} strains) for a meaningful decomposition"
        )
        return sec

    rows = []
    for run, wide in [(run_a, wide_a), (run_b, wide_b)]:
        centered = wide.to_numpy(dtype=float)
        centered = centered - centered.mean(axis=0, keepdims=True)
        singular = numpy.linalg.svd(centered, compute_uv=False)
        var = singular**2
        frac = var / var.sum()
        row = {"run": run.label}
        cumulative = 0.0
        for i in range(min(N_PCS, len(frac))):
            cumulative += float(frac[i])
            row[f"pc{i + 1}_frac_var"] = round(float(frac[i]), 4)
            row[f"pc1_to_{i + 1}_cum_frac_var"] = round(cumulative, 4)
        rows.append(row)
    table = pd.DataFrame(rows)
    sec.tables["titer_matrix_coherence"] = table

    a = table[table["run"] == run_a.label]
    b = table[table["run"] == run_b.label]
    cum_col = f"pc1_to_{min(N_PCS, wide_a.shape[1])}_cum_frac_var"
    if cum_col in table.columns and not a.empty and not b.empty:
        va, vb = float(a[cum_col].iloc[0]), float(b[cum_col].iloc[0])
        summary.add(
            "coherence",
            cum_col,
            value_a=va,
            value_b=vb,
            n=int(wide_a.shape[0] * wide_a.shape[1]),
            higher_is_better=True,
            note=f"identical {wide_a.shape[0]} sera x {wide_a.shape[1]} strains submatrix",
        )
        better = run_b.label if vb > va else run_a.label
        sec.note(
            f"""
            On an identical {wide_a.shape[0]} x {wide_a.shape[1]} submatrix, the first
            {min(N_PCS, wide_a.shape[1])} components explain {100 * va:.1f}% of variance
            in {run_a.label} and {100 * vb:.1f}% in {run_b.label}: **{better}** gives
            the more coherent titer matrix. The gap is usually small; treat it as
            supporting evidence, and note that a run could in principle score well here
            by over-smoothing genuine strain-specific signal.
            """
        )
    return sec


# ======================================================================================
# Section 7: dataset context -- how much room is there for collapsing to help?
# ======================================================================================


def section_dataset_context(run_a: Run, run_b: Run, git: Git, summary: Summary) -> Section:
    """Library structure and count depth: the prior on how much collapsing can help.

    The case for collapsing barcodes is that summing counts raises the number of
    virions behind each measurement and so cuts counting noise. That gain scales with
    how many barcode-wells are actually count-limited, which is a property of the
    dataset, not of either run. A library whose barcodes all have deep counts has
    little to gain; one with a long low-abundance tail has a lot.
    """
    sec = Section(
        "context",
        "7. Dataset context: how much room is there for collapsing to help?",
        lead=(
            "Library structure and sequencing depth. These are properties of the data, "
            "identical in both runs, that set an upper bound on the possible benefit."
        ),
    )
    rows = []
    for run in (run_a, run_b):
        per_strain = run.n_barcodes_per_strain
        if per_strain.empty:
            continue
        rows.append(
            {
                "run": run.label,
                "collapse_strain_barcodes": run.collapsed,
                "n_barcodes": len(run.barcode_to_strain),
                "n_strains": int(per_strain.size),
                "median_barcodes_per_strain": float(per_strain.median()),
                "min_barcodes_per_strain": int(per_strain.min()),
                "max_barcodes_per_strain": int(per_strain.max()),
            }
        )
    if rows:
        table = pd.DataFrame(rows)
        sec.tables["library_structure"] = table
        sec.tables["barcodes_per_strain_distribution"] = (
            run_a.n_barcodes_per_strain.value_counts()
            .sort_index()
            .rename_axis("n_barcodes")
            .reset_index(name="n_strains")
        )
        if float(table["max_barcodes_per_strain"].max()) <= 1:
            sec.note(
                "**Every strain has a single barcode, so collapsing is a no-op and this "
                "comparison cannot show a difference.** Check that the right commits "
                "and library are being compared."
            )
        summary.add(
            "context",
            "median_barcodes_per_strain",
            value_a=float(table["median_barcodes_per_strain"].iloc[0]),
            value_b=float(table["median_barcodes_per_strain"].iloc[-1]),
            note="library property, identical in both runs unless the library changed",
        )
    return sec


def add_count_depth(sec: Section, run: Run, git: Git, summary: Summary, max_files: int) -> None:
    """Optional: read the raw barcode counts to quantify the count-limited tail.

    Opt-in behind --count-depth because these are thousands of small per-well CSVs and
    reading them all is by far the slowest part of the script. Counts are pre-collapse
    and therefore identical in both runs, so they are read from one run only.
    """
    paths = [p for p in git.ls_tree(run.commit) if p.startswith("results/barcode_counts/")]
    if not paths:
        sec.note("`--count-depth` requested but no results/barcode_counts/ files are committed.")
        return
    sampled = paths[:: max(1, len(paths) // max_files)][:max_files]
    per_barcode, per_strain = [], []
    for path in sampled:
        frame = git.read_csv(run.commit, path)
        if frame is None or "count" not in frame.columns or "barcode" not in frame.columns:
            continue
        frame = frame.copy()
        frame["strain"] = frame["barcode"].map(run.barcode_to_strain)
        frame = frame.dropna(subset=["strain"])
        if frame.empty:
            continue
        per_barcode.append(frame["count"])
        per_strain.append(frame.groupby("strain")["count"].sum())
    if not per_barcode:
        sec.note("`--count-depth` requested but no readable barcode-count files were found.")
        return
    bc = pd.concat(per_barcode, ignore_index=True)
    st = pd.concat(per_strain, ignore_index=True)

    thresholds = [10, 50, 100, 500]
    table = pd.DataFrame(
        [
            {
                "level": "per barcode (uncollapsed)",
                "n_observations": len(bc),
                "median_count": float(bc.median()),
                **{f"frac_below_{t}": round(float((bc < t).mean()), 4) for t in thresholds},
            },
            {
                "level": "per strain, summed (collapsed)",
                "n_observations": len(st),
                "median_count": float(st.median()),
                **{f"frac_below_{t}": round(float((st < t).mean()), 4) for t in thresholds},
            },
        ]
    )
    sec.tables["count_depth"] = table
    frac_low = float((bc < 50).mean())
    summary.add(
        "context",
        "frac_barcode_wells_below_50_counts",
        value_a=round(frac_low, 6),
        value_b=round(frac_low, 6),
        n=len(bc),
        note=f"sampled {len(sampled)} of {len(paths)} count files; dataset property",
    )
    sec.note(
        f"""
        Sampled {len(sampled)} of {len(paths)} count files. Median count per barcode-well
        is {float(bc.median()):.0f} and {100 * frac_low:.1f}% of barcode-wells fall below
        50 counts. """
        + (
            "Most barcodes are already deep enough that Poisson noise is small, so the "
            "counting-noise argument for collapsing applies mainly to the low-abundance "
            "tail rather than to the typical measurement."
            if float(bc.median()) > 500
            else "A large share of barcode-wells are count-limited, so collapsing has "
            "substantial room to reduce counting noise."
        )
    )


# ======================================================================================
# Report rendering
# ======================================================================================


def _histogram(
    frame: pd.DataFrame, column: str, by: str, bins: int = 60
) -> pd.DataFrame:
    """Bin `column` within each group of `by`, returning one row per bin.

    A binned-and-counted chart needs only the counts, so binning here rather than letting
    `vega-lite` do it means the embedded data is a few dozen rows instead of every
    observation -- and unlike plotting a downsample, no observation is discarded. Shared bin
    edges across groups so the distributions are directly comparable.
    """
    values = frame[column].dropna()
    if values.empty:
        return pd.DataFrame(columns=[by, "bin_start", "bin_end", "count"])
    edges = numpy.linspace(float(values.min()), float(values.max()), bins + 1)
    rows = []
    for key, grp in frame.groupby(by):
        counts, _ = numpy.histogram(grp[column].dropna(), bins=edges)
        for i, count in enumerate(counts):
            if count:
                rows.append(
                    {
                        by: key,
                        "bin_start": round(float(edges[i]), 4),
                        "bin_end": round(float(edges[i + 1]), 4),
                        "count": int(count),
                    }
                )
    return pd.DataFrame(rows)


def _lookup_encode(
    frame: pd.DataFrame, columns: list[str], carry: dict[str, str] | None = None
) -> tuple[pd.DataFrame, list]:
    """Replace repeated string columns with integer codes plus small lookup tables.

    `altair` embeds the plotted frame as inline JSON, so a long strain or serum name costs
    its full length on every row. Each named column is replaced by an integer surrogate and
    re-joined client-side with `transform_lookup`, which is the pattern
    `seqneut-pipeline/CLAUDE.md` ("Keeping Altair plots small") prescribes for a column
    functionally determined by a key. Here the surrogate *is* the key, so the join is 1:1
    by construction.

    `carry` maps an extra column to the column that determines it -- `group` is determined
    by `serum`, so it rides along in the serum lookup instead of being repeated per row.
    Returns the narrowed frame and the `transform_lookup` arguments to apply to the chart.
    A tooltip or encoding referring to a looked-up field must state its type explicitly
    (`alt.Tooltip("virus:N")`), since the field is absent from the embedded data.
    """
    carry = carry or {}
    out = frame.copy()
    lookups = []
    for column in columns:
        if column not in out.columns:
            continue
        extra = [c for c, src in carry.items() if src == column and c in out.columns]
        table = (
            out[[column, *extra]]
            .drop_duplicates(subset=[column])
            .reset_index(drop=True)
            .reset_index(names=f"{column}_i")
        )
        assert table[column].is_unique, f"{column} lookup is not 1:1"
        codes = dict(zip(table[column], table[f"{column}_i"]))
        out[f"{column}_i"] = out[column].map(codes)
        out = out.drop(columns=[column, *extra])
        lookups.append(
            {
                "lookup": f"{column}_i",
                "from_": alt.LookupData(
                    table, key=f"{column}_i", fields=[column, *extra]
                ),
            }
        )
    return out, lookups


def _apply_lookups(chart, lookups: list):
    for kwargs in lookups:
        chart = chart.transform_lookup(**kwargs)
    return chart


def _downsample(frame: pd.DataFrame, limit: int) -> pd.DataFrame:
    """Cap rows sent to a chart, deterministically. ADAPT: raise `limit` if a chart
    looks sparse; Altair embeds the data in the HTML so this trades file size for
    completeness. Any truncation is reported in the chart caption."""
    if len(frame) <= limit:
        return frame
    step = len(frame) // limit + 1
    return frame.iloc[::step].copy()


REPORT_TEMPLATE_SOURCE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<script src="https://cdn.jsdelivr.net/npm/vega@{{ vega_version }}"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@{{ vegalite_version }}"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         max-width: 1080px; margin: 2rem auto; padding: 0 1.5rem; color: #1a1a1a; line-height: 1.55; }
  h1 { border-bottom: 3px solid #333; padding-bottom: .4rem; }
  h2 { margin-top: 2.6rem; border-bottom: 1px solid #ddd; padding-bottom: .3rem; }
  h3 { margin-top: 1.6rem; font-size: 1.02rem; color: #444; }
  table { border-collapse: collapse; margin: 1rem 0; font-size: .86rem; }
  th, td { border: 1px solid #ddd; padding: .32rem .6rem; text-align: right; }
  th { background: #f4f4f4; text-align: left; }
  td:first-child, th:first-child { text-align: left; }
  .note { background: #fbf7e8; border-left: 4px solid #d8b64a; padding: .7rem 1rem; margin: .9rem 0; }
  .skip { background: #f0f0f0; border-left: 4px solid #999; padding: .7rem 1rem; color: #555; }
  .caveat { background: #fdecea; border-left: 4px solid #d9534f; padding: .7rem 1rem; margin: .9rem 0; }
  .meta { background: #f6f8fa; padding: .9rem 1.1rem; font-size: .88rem; border-radius: 5px; }
  .meta code { background: #e8eaed; padding: .1rem .3rem; border-radius: 3px; }
  .chart { margin: 1.2rem 0; }
  .caption { font-size: .85rem; color: #555; margin-bottom: .35rem; }
  .toc a { display: block; margin: .18rem 0; }
  footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ddd; font-size: .82rem; color: #666; }
</style>
</head>
<body>
<h1>{{ title }}</h1>

<div class="meta">
<strong>Run A ({{ run_a.label }})</strong>: <code>{{ run_a.commit }}</code>
 = <code>{{ run_a.sha }}</code> &mdash; {{ run_a.subject }}<br>
<strong>Run B ({{ run_b.label }})</strong>: <code>{{ run_b.commit }}</code>
 = <code>{{ run_b.sha }}</code> &mdash; {{ run_b.subject }}<br>
<strong>collapse_strain_barcodes</strong>: {{ run_a.label }} = {{ run_a.collapsed }},
 {{ run_b.label }} = {{ run_b.collapsed }}<br>
<strong>Repeat sera</strong> (measured on &ge;{{ min_plates }} plates, same plate set in
 both runs): {{ n_repeat_sera }} of {{ n_sera }}<br>
<strong>Comparison id</strong>: <code>{{ comparison_id }}</code>
</div>

{% for caveat in caveats %}<div class="caveat">{{ caveat }}</div>{% endfor %}

<h2>Contents</h2>
<div class="toc">
{% for sec in sections %}<a href="#{{ sec.key }}">{{ sec.title }}</a>{% endfor %}
<a href="#how-to-read">How to read this report</a>
</div>

{% for sec in sections %}
<h2 id="{{ sec.key }}">{{ sec.title }}</h2>
{% if sec.lead %}<p>{{ sec.lead }}</p>{% endif %}
{% if sec.skipped %}
  <div class="skip"><strong>Skipped:</strong> {{ sec.skipped }}</div>
{% endif %}
{% for note in sec.notes %}<div class="note">{{ note | md }}</div>{% endfor %}
{% for name, html in sec.rendered_tables %}
  <h3>{{ name }}</h3>
  {{ html }}
{% endfor %}
{% for caption, div_id, spec in sec.rendered_charts %}
  <div class="chart">
    <div class="caption">{{ caption }}</div>
    <div id="{{ div_id }}"></div>
    <script>vegaEmbed('#{{ div_id }}', {{ spec }}, {actions: {export: true, source: false, editor: false}});</script>
  </div>
{% endfor %}
{% endfor %}

<h2 id="how-to-read">How to read this report</h2>
<p>Metrics worth trusting most, in order:</p>
<ol>
<li><strong>Plate-to-plate reproducibility on the matched set</strong> (section 2). The only
    metric based on genuinely independent remeasurement. Prefer
    <code>sd_log2_centered</code> over <code>pearson_r</code>.</li>
<li><strong>Adjudication of run-unique titers</strong> (section 3). Turns the confounded
    "who drops fewer titers" question into a testable one.</li>
<li><strong>Variance components and curve-fit quality</strong> (section 4). Explains the
    mechanism behind any difference.</li>
<li><strong>Titer-matrix coherence</strong> (section 6). Weaker, but uses every serum
    rather than only the repeat ones.</li>
</ol>
<p>Metrics that are <em>not</em> quality signals on their own: raw counts of dropped
titers (section 5), number of curves fit (section 4b), and anything computed on each
run's own retained set.</p>
<p>Every number above comes from a CSV in <code>tables/</code>, and the headline metrics
are in <code>summary.csv</code>, whose schema is fixed across comparisons so several can
be stacked and compared.</p>

{% if omitted_tables %}
<p><strong>Bulk tables omitted from <code>tables/</code></strong> to keep the output small.
Re-run with <code>--full-tables</code> to write them; each is largely re-derivable from the
pipeline output committed at the two commits.</p>
<ul>
{% for name, rows in omitted_tables %}<li><code>{{ name }}.csv</code> ({{ "{:,}".format(rows) }} rows)</li>{% endfor %}
</ul>
{% endif %}

<footer>
Generated by <code>{{ script_name }}</code>. Charts render via the vega CDN, so a network
connection is needed to view them; the underlying data is embedded in this file.
</footer>
</body>
</html>
"""


def _md_inline(text: str) -> str:
    """Minimal markdown: **bold** and `code`, with everything else escaped."""
    import html
    import re

    out = html.escape(str(text))
    out = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", out, flags=re.S)
    out = re.sub(r"`(.+?)`", r"<code>\1</code>", out, flags=re.S)
    return out.replace("\n", " ")


def render_report(
    path: Path,
    *,
    comparison_id: str,
    run_a: Run,
    run_b: Run,
    sections: list[Section],
    caveats: list[str],
    n_repeat_sera: int,
    n_sera: int,
    omitted_tables: list | None = None,
) -> None:
    schema = getattr(alt, "SCHEMA_VERSION", "v5.0.0").lstrip("v")
    vegalite_major = schema.split(".")[0]
    vega_version = {"5": "5", "6": "6"}.get(vegalite_major, "5")

    prepared = []
    for i, sec in enumerate(sections):
        sec.rendered_tables = [
            (
                name.replace("_", " "),
                frame.to_html(index=False, na_rep="", float_format=lambda v: f"{v:g}"),
            )
            for name, frame in sec.tables.items()
            # Very long per-observation tables go to CSV only; the HTML would be unusable.
            if len(frame) <= 60
        ]
        sec.rendered_charts = []
        for j, (caption, chart) in enumerate(sec.charts):
            try:
                spec = chart.to_json(indent=None)
            except Exception as err:  # a bad chart must not lose the whole report
                sec.notes.append(f"chart {caption!r} could not be rendered: {err}")
                continue
            sec.rendered_charts.append((caption, f"chart-{i}-{j}", spec))
        prepared.append(sec)

    # The `md` filter has to exist on the environment before the template is compiled.
    env = Environment(autoescape=False)
    env.filters["md"] = _md_inline
    html = env.from_string(REPORT_TEMPLATE_SOURCE).render(
        title=f"seqneut-pipeline run comparison: {run_a.label} vs {run_b.label}",
        comparison_id=comparison_id,
        run_a=run_a,
        run_b=run_b,
        sections=prepared,
        caveats=[_md_inline(c) for c in caveats],
        n_repeat_sera=n_repeat_sera,
        n_sera=n_sera,
        omitted_tables=omitted_tables or [],
        min_plates=MIN_PLATES_FOR_REPRODUCIBILITY,
        vega_version=vega_version,
        vegalite_version=schema,
        script_name=Path(__file__).name,
    )
    path.write_text(html)


# ======================================================================================
# Cross-run sanity checks
# ======================================================================================


def cross_run_caveats(run_a: Run, run_b: Run) -> list[str]:
    """Standing caveats detected by comparing the two runs' configs.

    These are things the script cannot correct for and that a reader must know about.
    """
    caveats: list[str] = []

    # Titers must be on the same scale to be comparable at all.
    for column in ("titer_as", "titer_units"):
        vals = {}
        for run in (run_a, run_b):
            if not run.titers.empty and column in run.titers.columns:
                vals[run.label] = set(run.titers[column].dropna().unique())
        if len(vals) == 2 and list(vals.values())[0] != list(vals.values())[1]:
            caveats.append(
                f"**The two runs report different `{column}` ({vals}).** Their titers "
                "are not on the same scale and the concordance and reproducibility "
                "numbers below are not meaningful. Fix this before drawing conclusions."
            )

    # The count-threshold confound: collapsing changes what count thresholds mean.
    if run_a.collapsed != run_b.collapsed:
        collapsed_run = run_b if run_b.collapsed else run_a
        other = run_a if run_b.collapsed else run_b
        count_keys = [
            "avg_barcode_counts_per_well",
            "min_no_serum_count_per_viral_barcode_well",
            "min_neut_standard_count_per_well",
        ]
        same = []
        for key in count_keys:
            va = _find_threshold(collapsed_run.config, key)
            vb = _find_threshold(other.config, key)
            if va is not None and va == vb:
                same.append(f"{key}={va}")
        if same:
            n_bc = collapsed_run.n_barcodes_per_strain
            factor = float(n_bc.median()) if not n_bc.empty else float("nan")
            caveats.append(
                f"**Count-threshold confound.** `{collapsed_run.label}` collapses strain "
                f"barcodes but keeps the same count thresholds as `{other.label}` "
                f"({', '.join(same)}). Those thresholds apply to summed counts once "
                f"barcodes are collapsed, so they are effectively about {factor:.1f}x "
                "looser in the collapsed run. Part of any retention difference is this, "
                "not the collapsing itself. Separating the two needs a third run with "
                "the thresholds rescaled."
            )

    if run_a.collapsed == run_b.collapsed:
        caveats.append(
            f"Both runs have `collapse_strain_barcodes = {run_a.collapsed}`, so they "
            "differ in some other way. The section leads still apply, but wording that "
            "refers to collapsing may not describe what is actually being compared."
        )

    lib_a, lib_b = run_a.strains, run_b.strains
    if lib_a and lib_b and lib_a != lib_b:
        caveats.append(
            f"**The viral libraries differ**: {len(lib_a - lib_b)} strains only in "
            f"{run_a.label}, {len(lib_b - lib_a)} only in {run_b.label}. These are "
            "attributed as `virus_absent_from_run` in section 1a and excluded from "
            "quality metrics."
        )

    plates_a, plates_b = run_a.plates_by_serum, run_b.plates_by_serum
    differing = [k for k in set(plates_a) & set(plates_b) if plates_a[k] != plates_b[k]]
    only_a, only_b = set(plates_a) - set(plates_b), set(plates_b) - set(plates_a)
    if differing or only_a or only_b:
        caveats.append(
            f"**The runs differ in plate coverage**: {len(differing)} sera measured on "
            f"different plate sets, {len(only_a)} sera only in {run_a.label}, "
            f"{len(only_b)} only in {run_b.label}. Sera with differing plate sets are "
            "excluded from the matched set, since the amount of averaging behind their "
            "titers is not comparable."
        )
    return caveats


def _find_threshold(config: dict, key: str):
    """Search a config recursively for a threshold value, returning the first found.

    Thresholds live under per-plate blocks built from YAML anchors, so the same key
    appears many times. ADAPT: this returns the first value found and so assumes the
    plates share defaults, which is true here; a config with genuinely per-plate
    thresholds would need this to return the full set instead.
    """
    if isinstance(config, dict):
        for k, v in config.items():
            if k == key and not isinstance(v, (dict, list)):
                return v
            found = _find_threshold(v, key)
            if found is not None:
                return found
    elif isinstance(config, list):
        for item in config:
            found = _find_threshold(item, key)
            if found is not None:
                return found
    return None


# ======================================================================================
# Main
# ======================================================================================


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="See the module docstring for why each comparison is built the way it is.",
    )
    parser.add_argument("--run-a", required=True, help="commit-ish for the baseline run")
    parser.add_argument("--run-b", required=True, help="commit-ish for the run under test")
    parser.add_argument("--label-a", default=None, help="short name for run A")
    parser.add_argument("--label-b", default=None, help="short name for run B")
    parser.add_argument("--comparison-id", default=None, help="identifier used in summary.csv")
    parser.add_argument("--repo", default=".", help="path to the git repository")
    parser.add_argument("--out", default=None, help="output directory")
    parser.add_argument(
        "--count-depth",
        action="store_true",
        help="also read results/barcode_counts/ to quantify the count-limited tail "
        "(slow: thousands of small files)",
    )
    parser.add_argument(
        "--count-depth-max-files",
        type=int,
        default=400,
        help="how many count files to sample when --count-depth is given",
    )
    parser.add_argument(
        "--full-tables",
        action="store_true",
        help=f"also write bulk per-observation tables (>= {BULK_TABLE_MIN_ROWS:,} rows). "
        "These are mostly re-derivable from the committed pipeline output at the two "
        "commits, so they are omitted by default to keep the output small; any omitted "
        "table is named in the console output and the report footer",
    )
    args = parser.parse_args(argv)

    label_a = args.label_a or args.run_a.replace("/", "_")
    label_b = args.label_b or args.run_b.replace("/", "_")
    if label_a == label_b:
        label_a, label_b = f"{label_a}_a", f"{label_b}_b"
    comparison_id = args.comparison_id or f"{label_a}__vs__{label_b}"
    out_dir = Path(args.out) if args.out else Path(__file__).parent / comparison_id
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)

    try:
        alt.data_transformers.disable_max_rows()
    except Exception:
        pass

    git = Git(Path(args.repo).resolve())
    print(f"reading run A ({label_a}) from {args.run_a} ...", file=sys.stderr)
    run_a = load_run(git, args.run_a, label_a)
    print(f"reading run B ({label_b}) from {args.run_b} ...", file=sys.stderr)
    run_b = load_run(git, args.run_b, label_b)

    summary = Summary(comparison_id, run_a, run_b)
    caveats = cross_run_caveats(run_a, run_b) + run_a.notes + run_b.notes

    matched, matched_notes = build_matched_set(run_a, run_b)
    print(f"matched set: {len(matched)} serum-virus pairs", file=sys.stderr)

    repeat_sera = {
        k
        for k, plates in run_a.plates_by_serum.items()
        if len(plates) >= MIN_PLATES_FOR_REPRODUCIBILITY
        and run_b.plates_by_serum.get(k) == plates
    }

    # Each analysis is isolated: a failure becomes a visible note, never a lost report.
    sections: list[Section] = []
    builders = [
        ("coverage", lambda: section_coverage(run_a, run_b, summary)),
        ("concordance", lambda: section_concordance(run_a, run_b, matched, summary)),
        ("reproducibility", lambda: section_reproducibility(run_a, run_b, matched, summary)),
        ("adjudication", lambda: section_adjudicate(run_a, run_b, summary)),
        ("variance", lambda: section_variance_components(run_a, run_b, summary)),
        ("curvefits", lambda: section_curvefit_quality(run_a, run_b, summary)),
        ("retention", lambda: section_retention(run_a, run_b, summary)),
        ("coherence", lambda: section_coherence(run_a, run_b, matched, summary)),
        ("context", lambda: section_dataset_context(run_a, run_b, git, summary)),
    ]
    for key, builder in builders:
        print(f"  section {key} ...", file=sys.stderr)
        try:
            sections.append(builder())
        except Exception:
            failed = Section(key, f"{key} (failed)")
            failed.skipped = "this analysis raised an exception; see the traceback below"
            failed.note(f"`{traceback.format_exc(limit=6)}`")
            sections.append(failed)
            print(f"    FAILED: see report", file=sys.stderr)

    if args.count_depth:
        print("  reading barcode counts ...", file=sys.stderr)
        context = next((s for s in sections if s.key == "context"), None)
        if context is not None:
            try:
                add_count_depth(context, run_a, git, summary, args.count_depth_max_files)
            except Exception as err:
                context.note(f"count-depth analysis failed: `{err}`")

    for note in matched_notes:
        for sec in sections:
            if sec.key == "concordance":
                sec.note(note)

    # Write every table, including the long ones the HTML omits. Bulk per-observation
    # tables are held back unless --full-tables, but are always named so their absence is
    # visible rather than silent.
    omitted_tables: list[tuple[str, int]] = []
    for sec in sections:
        for name, frame in sec.tables.items():
            if len(frame) >= BULK_TABLE_MIN_ROWS and not args.full_tables:
                omitted_tables.append((name, len(frame)))
                continue
            frame.to_csv(out_dir / "tables" / f"{name}.csv", index=False)

    summary_frame = summary.frame()
    summary_frame.to_csv(out_dir / "summary.csv", index=False)

    render_report(
        out_dir / "report.html",
        comparison_id=comparison_id,
        run_a=run_a,
        run_b=run_b,
        sections=sections,
        caveats=caveats,
        n_repeat_sera=len(repeat_sera),
        n_sera=len(run_a.sera | run_b.sera),
        omitted_tables=omitted_tables,
    )

    # Console summary: enough to know the answer without opening the report.
    print(f"\n{'=' * 78}\n{comparison_id}\n{'=' * 78}")
    print(f"A = {label_a:24s} {run_a.sha}  collapse={run_a.collapsed}")
    print(f"B = {label_b:24s} {run_b.sha}  collapse={run_b.collapsed}")
    print(f"matched serum-virus pairs: {len(matched)};  repeat sera: {len(repeat_sera)}")
    for caveat in caveats:
        print(f"\n! {textwrap.shorten(caveat.replace('**', ''), 400)}")
    headline = summary_frame[
        summary_frame["metric"].isin(
            [
                "n_titers",
                "n_titers_unique_to_run",
                "median_abs_log2_fc",
                "sd_log2_centered",
                "median_r2",
                "per_serum_qc_titer_drops",
                "frac_titers_with_1_replicate",
            ]
        )
    ]
    if not headline.empty:
        print("\nheadline metrics:")
        print(
            headline[["section", "metric", "stratum", "value_a", "value_b", "p_value"]]
            .to_string(index=False)
        )
    if omitted_tables:
        print("\nbulk tables omitted (re-run with --full-tables to write them):")
        for name, rows in omitted_tables:
            print(f"  {name}.csv ({rows:,} rows)")
    skipped = [s for s in sections if s.skipped]
    if skipped:
        print("\nsections skipped:")
        for sec in skipped:
            print(f"  {sec.key}: {sec.skipped}")
    print(f"\nreport:  {out_dir / 'report.html'}")
    print(f"summary: {out_dir / 'summary.csv'}")
    print(f"tables:  {out_dir / 'tables'}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
