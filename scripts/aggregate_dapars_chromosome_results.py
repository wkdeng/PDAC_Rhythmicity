#!/usr/bin/env python3
"""Aggregate per-chromosome DaPars2 outputs into TCGA-style sample folders.

DaPars2_Multi_Sample_Multi_Chr.py writes one directory per sample/chromosome:

    {sample}_chr1/{sample}_result_result_temp.chr1.txt

This script combines those chromosome files into:

    results/{sample}/{sample}_result_All_Prediction_Results.txt

It preserves the chromosome-level outputs under ``result_by_chromosome`` and can
optionally move an older mixed layout into that structure.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


DEFAULT_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
MISSING = "NA"

CHROM_DIR_RE = re.compile(r"^(?P<sample>.+)_(?P<chrom>chr(?:[0-9]+|X|Y|M|MT))$")

TCGA_STYLE_COLUMNS = [
    "Gene",
    "fit_value",
    "Predicted_Proximal_APA",
    "Loci",
    "A_1_long_exp",
    "A_1_short_exp",
    "A_1_PDUI",
    "B_1_long_exp",
    "B_1_short_exp",
    "B_1_PDUI",
    "Group_A_Mean_PDUI",
    "Group_B_Mean_PDUI",
    "PDUI_Group_diff",
    "P_val",
    "adjusted.P_val",
    "Pass_Filter",
]


@dataclass
class SampleSummary:
    sample: str
    result_file: Path
    manifest_file: Path
    completed_file: Path
    rows_written: int
    chromosomes_expected: int
    chromosomes_observed: int
    missing_chromosomes: list[str]
    filled_expression_columns: bool
    pdui_source_columns: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate per-chromosome DaPars2 result_temp files into TCGA-style sample outputs."
    )
    parser.add_argument(
        "--cohort-root",
        type=Path,
        help="Cohort polyAPA root, e.g. data/polyAPA/GTEx. Used to derive defaults.",
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        help="Existing mixed root containing {sample}_chr* directories. Defaults are auto-detected.",
    )
    parser.add_argument(
        "--chrom-root",
        type=Path,
        help="Directory for chromosome-level folders. Default: {cohort-root}/result_by_chromosome.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Directory for final TCGA-style sample folders. Default: {cohort-root}/results.",
    )
    parser.add_argument(
        "--chromosomes",
        type=Path,
        help="File containing expected chromosomes, one per line. Default: chr1-22, chrX, chrY.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        dest="samples",
        help="Sample ID to aggregate. May be supplied multiple times. Defaults to all discovered samples.",
    )
    parser.add_argument(
        "--samples-file",
        type=Path,
        help="Optional file with one sample ID per line.",
    )
    parser.add_argument(
        "--move-chromosome-dirs",
        action="store_true",
        help="Move top-level {sample}_chr* directories from input-root into chrom-root before aggregation.",
    )
    parser.add_argument(
        "--allow-missing-chroms",
        action="store_true",
        help="Aggregate available chromosomes and report missing chromosomes instead of failing.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing aggregate result files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned moves and aggregations without writing outputs.",
    )
    parser.add_argument(
        "--write-summary",
        action="store_true",
        help="Write output-root/aggregation_summary.tsv after aggregating all requested samples.",
    )
    return parser.parse_args()


def require_path(name: str, value: Path | None) -> Path:
    if value is None:
        raise SystemExit(f"Missing required path: {name}. Provide --cohort-root or {name}.")
    return value


def has_chromosome_dirs(path: Path) -> bool:
    return path.exists() and any(CHROM_DIR_RE.match(child.name) for child in path.iterdir() if child.is_dir())


def resolve_roots(args: argparse.Namespace) -> tuple[Path | None, Path, Path, Path]:
    cohort_root = args.cohort_root.resolve() if args.cohort_root else None

    chrom_root = args.chrom_root
    if chrom_root is None and cohort_root is not None:
        chrom_root = cohort_root / "result_by_chromosome"
    chrom_root = require_path("--chrom-root", chrom_root).resolve()

    output_root = args.output_root
    if output_root is None and cohort_root is not None:
        output_root = cohort_root / "results"
    output_root = require_path("--output-root", output_root).resolve()

    input_root = args.input_root
    if input_root is None:
        if has_chromosome_dirs(chrom_root):
            input_root = chrom_root
        elif has_chromosome_dirs(output_root):
            input_root = output_root
        elif cohort_root is not None and has_chromosome_dirs(cohort_root):
            input_root = cohort_root
        else:
            input_root = chrom_root
    input_root = input_root.resolve()

    return cohort_root, input_root, chrom_root, output_root


def read_chromosomes(path: Path | None) -> list[str]:
    if path is None:
        return DEFAULT_CHROMS
    chroms = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not chroms:
        raise SystemExit(f"No chromosomes found in {path}")
    return chroms


def read_samples(samples: list[str] | None, samples_file: Path | None) -> list[str]:
    out: list[str] = []
    if samples:
        out.extend(samples)
    if samples_file:
        out.extend(line.strip() for line in samples_file.read_text().splitlines() if line.strip())
    return sorted(dict.fromkeys(out))


def discover_samples(chrom_root: Path) -> list[str]:
    samples = set()
    if not chrom_root.exists():
        return []
    for child in chrom_root.iterdir():
        if not child.is_dir():
            continue
        match = CHROM_DIR_RE.match(child.name)
        if match:
            samples.add(match.group("sample"))
    return sorted(samples)


def move_chromosome_dirs(input_root: Path, chrom_root: Path, dry_run: bool) -> int:
    if input_root == chrom_root:
        print(f"Chromosome directories already use chrom-root: {chrom_root}")
        return 0
    if not input_root.exists():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")

    moved = 0
    for child in sorted(input_root.iterdir()):
        if not child.is_dir() or not CHROM_DIR_RE.match(child.name):
            continue
        dest = chrom_root / child.name
        if dest.exists():
            print(f"SKIP existing destination: {dest}")
            continue
        print(f"MOVE {child} -> {dest}")
        moved += 1
        if not dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(child), str(dest))
    return moved


def find_chromosome_file(chrom_root: Path, sample: str, chrom: str) -> Path | None:
    chrom_dir = chrom_root / f"{sample}_{chrom}"
    exact = chrom_dir / f"{sample}_result_result_temp.{chrom}.txt"
    if exact.exists():
        return exact
    if chrom_dir.exists():
        matches = sorted(chrom_dir.glob(f"*result_temp.{chrom}.txt"))
        if matches:
            return matches[0]
    return None


def select_pdui_column(fieldnames: Iterable[str], sample: str) -> str | None:
    fields = list(fieldnames)
    preferred = ["A_1_PDUI", f"{sample}_PDUI"]
    for col in preferred:
        if col in fields:
            return col
    pdui_cols = [col for col in fields if col.endswith("_PDUI")]
    if pdui_cols:
        sample_hits = [col for col in pdui_cols if sample in col]
        return sample_hits[0] if sample_hits else pdui_cols[0]
    return None


def clean_value(value: str | None, default: str = MISSING) -> str:
    if value is None or value == "":
        return default
    return str(value)


def normalize_row(row: dict[str, str], pdui_col: str | None) -> tuple[list[str], bool]:
    pdui = clean_value(row.get("A_1_PDUI") or (row.get(pdui_col) if pdui_col else None))

    a_long = row.get("A_1_long_exp")
    a_short = row.get("A_1_short_exp")
    b_long = row.get("B_1_long_exp")
    b_short = row.get("B_1_short_exp")
    filled_expression = any(value in (None, "") for value in [a_long, a_short, b_long, b_short])

    out = {
        "Gene": row.get("Gene"),
        "fit_value": row.get("fit_value"),
        "Predicted_Proximal_APA": row.get("Predicted_Proximal_APA"),
        "Loci": row.get("Loci"),
        "A_1_long_exp": a_long,
        "A_1_short_exp": a_short,
        "A_1_PDUI": pdui,
        "B_1_long_exp": b_long,
        "B_1_short_exp": b_short,
        "B_1_PDUI": row.get("B_1_PDUI") or pdui,
        "Group_A_Mean_PDUI": row.get("Group_A_Mean_PDUI") or pdui,
        "Group_B_Mean_PDUI": row.get("Group_B_Mean_PDUI") or pdui,
        "PDUI_Group_diff": row.get("PDUI_Group_diff") or "0",
        "P_val": row.get("P_val") or "1",
        "adjusted.P_val": row.get("adjusted.P_val") or "1",
        "Pass_Filter": row.get("Pass_Filter") or MISSING,
    }
    return [clean_value(out[col]) for col in TCGA_STYLE_COLUMNS], filled_expression


def collect_chromosome_files(
    chrom_root: Path,
    sample: str,
    chromosomes: list[str],
    allow_missing: bool,
) -> tuple[list[tuple[str, Path]], list[str]]:
    observed: list[tuple[str, Path]] = []
    missing: list[str] = []
    for chrom in chromosomes:
        chrom_file = find_chromosome_file(chrom_root, sample, chrom)
        if chrom_file is None:
            missing.append(chrom)
        else:
            observed.append((chrom, chrom_file))
    if missing and not allow_missing:
        raise FileNotFoundError(
            f"{sample}: missing {len(missing)} chromosome result files: {', '.join(missing)}"
        )
    return observed, missing


def aggregate_sample(
    sample: str,
    chrom_root: Path,
    output_root: Path,
    chromosomes: list[str],
    allow_missing: bool,
    force: bool,
    dry_run: bool,
) -> SampleSummary:
    observed, missing = collect_chromosome_files(chrom_root, sample, chromosomes, allow_missing)
    sample_out_dir = output_root / sample
    result_file = sample_out_dir / f"{sample}_result_All_Prediction_Results.txt"
    manifest_file = sample_out_dir / f"{sample}_aggregation_manifest.tsv"
    completed_file = sample_out_dir / "dapars_completed.txt"

    if result_file.exists() and not force and not dry_run:
        raise FileExistsError(f"Refusing to overwrite existing result file without --force: {result_file}")

    print(
        f"AGGREGATE {sample}: {len(observed)}/{len(chromosomes)} chromosomes -> {result_file}"
    )
    if missing:
        print(f"  Missing chromosomes: {', '.join(missing)}")
    if dry_run:
        return SampleSummary(
            sample=sample,
            result_file=result_file,
            manifest_file=manifest_file,
            completed_file=completed_file,
            rows_written=0,
            chromosomes_expected=len(chromosomes),
            chromosomes_observed=len(observed),
            missing_chromosomes=missing,
            filled_expression_columns=False,
            pdui_source_columns=[],
        )

    sample_out_dir.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    filled_expression_columns = False
    pdui_source_columns: set[str] = set()

    with result_file.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(TCGA_STYLE_COLUMNS)

        for chrom, chrom_file in observed:
            with chrom_file.open(newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if not reader.fieldnames:
                    raise ValueError(f"{chrom_file} has no header")
                pdui_col = select_pdui_column(reader.fieldnames, sample)
                if pdui_col is None:
                    raise ValueError(f"{chrom_file} has no PDUI column in header: {reader.fieldnames}")
                pdui_source_columns.add(pdui_col)

                for row in reader:
                    normalized, filled_expr = normalize_row(row, pdui_col)
                    writer.writerow(normalized)
                    rows_written += 1
                    filled_expression_columns = filled_expression_columns or filled_expr

            print(f"  {chrom}: {chrom_file}")

    summary = SampleSummary(
        sample=sample,
        result_file=result_file,
        manifest_file=manifest_file,
        completed_file=completed_file,
        rows_written=rows_written,
        chromosomes_expected=len(chromosomes),
        chromosomes_observed=len(observed),
        missing_chromosomes=missing,
        filled_expression_columns=filled_expression_columns,
        pdui_source_columns=sorted(pdui_source_columns),
    )
    write_manifest(summary)
    completed_file.write_text(
        "\n".join(
            [
                f"Completed at {datetime.now().isoformat(timespec='seconds')}",
                f"Aggregated rows: {rows_written}",
                f"Observed chromosomes: {len(observed)}/{len(chromosomes)}",
                f"Result file: {result_file}",
                "",
            ]
        )
    )
    return summary


def write_manifest(summary: SampleSummary) -> None:
    rows = [
        ("sample", summary.sample),
        ("result_file", str(summary.result_file)),
        ("rows_written", str(summary.rows_written)),
        ("chromosomes_expected", str(summary.chromosomes_expected)),
        ("chromosomes_observed", str(summary.chromosomes_observed)),
        ("missing_chromosomes", ",".join(summary.missing_chromosomes)),
        ("filled_expression_columns", str(summary.filled_expression_columns)),
        ("pdui_source_columns", ",".join(summary.pdui_source_columns)),
    ]
    with summary.manifest_file.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(rows)


def write_summary(output_root: Path, summaries: list[SampleSummary]) -> Path:
    summary_path = output_root / "aggregation_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "result_file",
                "rows_written",
                "chromosomes_expected",
                "chromosomes_observed",
                "missing_chromosomes",
                "filled_expression_columns",
                "pdui_source_columns",
            ]
        )
        for summary in summaries:
            writer.writerow(
                [
                    summary.sample,
                    summary.result_file,
                    summary.rows_written,
                    summary.chromosomes_expected,
                    summary.chromosomes_observed,
                    ",".join(summary.missing_chromosomes),
                    summary.filled_expression_columns,
                    ",".join(summary.pdui_source_columns),
                ]
            )
    return summary_path


def main() -> int:
    args = parse_args()
    _, input_root, chrom_root, output_root = resolve_roots(args)
    chromosomes = read_chromosomes(args.chromosomes)

    print(f"Input root: {input_root}")
    print(f"Chromosome root: {chrom_root}")
    print(f"Output root: {output_root}")
    print(f"Expected chromosomes: {', '.join(chromosomes)}")

    if args.move_chromosome_dirs:
        moved = move_chromosome_dirs(input_root, chrom_root, args.dry_run)
        print(f"Chromosome directories moved/planned: {moved}")

    requested_samples = read_samples(args.samples, args.samples_file)
    samples = requested_samples or discover_samples(chrom_root)
    if not samples:
        raise SystemExit(f"No samples found in chromosome root: {chrom_root}")
    print(f"Samples to aggregate: {len(samples)}")

    summaries = [
        aggregate_sample(
            sample=sample,
            chrom_root=chrom_root,
            output_root=output_root,
            chromosomes=chromosomes,
            allow_missing=args.allow_missing_chroms,
            force=args.force,
            dry_run=args.dry_run,
        )
        for sample in samples
    ]

    if args.write_summary and not args.dry_run:
        summary_path = write_summary(output_root, summaries)
        print(f"Summary written: {summary_path}")

    n_rows = sum(summary.rows_written for summary in summaries)
    print(f"Done. Samples: {len(summaries)}. Rows written: {n_rows}.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
