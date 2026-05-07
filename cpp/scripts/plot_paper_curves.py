#!/usr/bin/env python3
"""Plot paper-style ATE and cumulative-update-FLOP curves from benchmark CSVs."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import List, Tuple

APPROACH_ORDER = ["GNi-SPO-IGG", "GNi-SPO-LCG", "GNi-SPO", "GNi-IGG", "GNi-LCG", "GNi", "GN1"]


def infer_approach(path: Path) -> str:
    for name in APPROACH_ORDER:
        if path.stem == name or path.stem.endswith("_" + name) or path.stem.endswith("-" + name):
            return name
    return path.stem


def read_series(path: Path, ycol: str) -> Tuple[str, List[int], List[float]]:
    xs: List[int] = []
    ys: List[float] = []
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if "batch" not in (reader.fieldnames or []) or ycol not in (reader.fieldnames or []):
            raise SystemExit(f"{path} is missing batch/{ycol}")
        for row in reader:
            try:
                x = int(row["batch"])
                y = float(row[ycol])
            except ValueError:
                continue
            if math.isfinite(y):
                xs.append(x)
                ys.append(y)
    return infer_approach(path), xs, ys


def sorted_csvs(csv_dir: Path) -> List[Path]:
    files = [p for p in csv_dir.glob("*.csv") if p.name != "paper_summary.csv"]
    def key(p: Path) -> tuple:
        name = infer_approach(p)
        return (APPROACH_ORDER.index(name) if name in APPROACH_ORDER else len(APPROACH_ORDER), name)
    return sorted(files, key=key)


def plot_metric(csv_dir: Path, out_path: Path, ycol: str, ylabel: str, title: str) -> None:
    import matplotlib.pyplot as plt
    plt.figure()
    for path in sorted_csvs(csv_dir):
        label, xs, ys = read_series(path, ycol)
        if xs and ys:
            plt.plot(xs, ys, label=label)
    plt.xlabel("increment number")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_dir", type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--dataset", default=None)
    args = parser.parse_args()
    if not args.csv_dir.is_dir():
        raise SystemExit(f"Not a directory: {args.csv_dir}")
    out_dir = args.out_dir or args.csv_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    label = args.dataset or args.csv_dir.name
    plot_metric(args.csv_dir, out_dir / "ate_over_increments.png", "ate", "absolute trajectory error", f"ATE over increments: {label}")
    plot_metric(args.csv_dir, out_dir / "cumulative_update_flops.png", "cumulative_update_flops", "cumulative update FLOPs", f"Cumulative update FLOPs: {label}")
    print(f"Wrote {out_dir / 'ate_over_increments.png'}")
    print(f"Wrote {out_dir / 'cumulative_update_flops.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
