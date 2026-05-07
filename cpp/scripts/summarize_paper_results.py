#!/usr/bin/env python3
"""Summarize benchmark CSVs into the paper-style experimental table."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List

APPROACH_ORDER = [
    "GNi-SPO-IGG",
    "GNi-SPO-LCG",
    "GNi-SPO",
    "GNi-IGG",
    "GNi-LCG",
    "GNi",
    "GN1",
]

REQUIRED_COLUMNS = [
    "normalized_chi2",
    "ate",
    "solve_flops",
    "full_solve_flops",
    "update_flops",
]


def parse_float(value: str) -> float:
    value = (value or "").strip()
    if not value:
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def finite_values(rows: List[Dict[str, str]], column: str) -> List[float]:
    out: List[float] = []
    for row in rows:
        v = parse_float(row.get(column, ""))
        if math.isfinite(v):
            out.append(v)
    return out


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else math.nan


def final_value(rows: List[Dict[str, str]], column: str) -> float:
    for row in reversed(rows):
        v = parse_float(row.get(column, ""))
        if math.isfinite(v):
            return v
    return math.nan


def fmt_float(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    if value == 0:
        return "0"
    av = abs(value)
    if av < 1e-3 or av >= 1e4:
        return f"{value:.6e}"
    return f"{value:.6g}"


def fmt_intlike(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{int(round(value))}"


def infer_approach(path: Path) -> str:
    stem = path.stem
    for approach in APPROACH_ORDER:
        if stem == approach or stem.endswith("_" + approach) or stem.endswith("-" + approach):
            return approach
    return stem


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path} is missing required columns: {', '.join(missing)}")
        return list(reader)


def summarize_csv(path: Path, dataset: str) -> Dict[str, object]:
    rows = read_rows(path)
    approach = infer_approach(path)
    out: Dict[str, object] = {
        "dataset": dataset,
        "approach": approach,
        "source_csv": str(path),
        "increments": len(rows),
        "final_normalized_chi2": final_value(rows, "normalized_chi2"),
        "mean_normalized_chi2": mean(finite_values(rows, "normalized_chi2")),
        "final_ate": final_value(rows, "ate"),
        "mean_ate": mean(finite_values(rows, "ate")),
        "mean_solve_flops": mean(finite_values(rows, "solve_flops")),
        "mean_full_solve_flops": mean(finite_values(rows, "full_solve_flops")),
        "mean_update_flops": mean(finite_values(rows, "update_flops")),
        "final_cumulative_solve_flops": final_value(rows, "cumulative_solve_flops"),
        "final_cumulative_full_solve_flops": final_value(rows, "cumulative_full_solve_flops"),
        "final_cumulative_update_flops": final_value(rows, "cumulative_update_flops"),
    }
    if approach.startswith("GNi-SPO"):
        out["paper_solve_flops"] = f"{fmt_intlike(float(out['mean_solve_flops']))} ({fmt_intlike(float(out['mean_full_solve_flops']))})"
    else:
        out["paper_solve_flops"] = fmt_intlike(float(out["mean_solve_flops"]))
    out["paper_update_flops"] = fmt_intlike(float(out["mean_update_flops"]))
    return out


def approach_sort_key(row: Dict[str, object]) -> tuple:
    approach = str(row.get("approach", ""))
    try:
        return (APPROACH_ORDER.index(approach), approach)
    except ValueError:
        return (len(APPROACH_ORDER), approach)


def write_summary_csv(rows: List[Dict[str, object]], path: Path) -> None:
    columns = [
        "dataset",
        "approach",
        "increments",
        "final_normalized_chi2",
        "mean_normalized_chi2",
        "final_ate",
        "mean_ate",
        "mean_solve_flops",
        "mean_full_solve_flops",
        "mean_update_flops",
        "paper_solve_flops",
        "paper_update_flops",
        "final_cumulative_solve_flops",
        "final_cumulative_full_solve_flops",
        "final_cumulative_update_flops",
        "source_csv",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def write_markdown(rows: List[Dict[str, object]], path: Path) -> None:
    headers = ["dataset", "approach", "Nchi2 final", "Nchi2 mean", "ATE final", "ATE mean", "solve FLOPs", "update FLOPs"]
    with path.open("w") as f:
        f.write("# Paper-style benchmark summary\n\n")
        f.write("| " + " | ".join(headers) + " |\n")
        f.write("|" + "|".join(["---"] * len(headers)) + "|\n")
        for row in rows:
            f.write(
                "| "
                + " | ".join([
                    str(row["dataset"]),
                    str(row["approach"]),
                    fmt_float(float(row["final_normalized_chi2"])),
                    fmt_float(float(row["mean_normalized_chi2"])),
                    fmt_float(float(row["final_ate"])),
                    fmt_float(float(row["mean_ate"])),
                    str(row["paper_solve_flops"]),
                    str(row["paper_update_flops"]),
                ])
                + " |\n"
            )
        f.write("\nFor SPO variants, solve FLOPs are shown as partial solve FLOPs with full-solve replacement FLOPs in parentheses.\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_dir", type=Path, help="Directory containing benchmark CSV files")
    parser.add_argument("--dataset", default=None, help="Dataset label to write in the summary")
    parser.add_argument("--out-csv", type=Path, default=None, help="Output summary CSV path")
    parser.add_argument("--out-md", type=Path, default=None, help="Output markdown summary path")
    args = parser.parse_args()

    if not args.csv_dir.is_dir():
        raise SystemExit(f"Not a directory: {args.csv_dir}")
    dataset = args.dataset or args.csv_dir.name
    csv_paths = sorted(p for p in args.csv_dir.glob("*.csv") if p.name != "paper_summary.csv")
    if not csv_paths:
        raise SystemExit(f"No benchmark CSV files found in {args.csv_dir}")

    rows = [summarize_csv(path, dataset) for path in csv_paths]
    rows.sort(key=approach_sort_key)
    out_csv = args.out_csv or (args.csv_dir / "paper_summary.csv")
    out_md = args.out_md or (args.csv_dir / "paper_summary.md")
    write_summary_csv(rows, out_csv)
    write_markdown(rows, out_md)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
