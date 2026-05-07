from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

DatasetName = Literal["MIT", "Intel", "CSAIL", "FR079", "FRH", "MITP"]
LinearBackend = Literal["auto", "legacy", "cholmod"]
GatingMode = Literal["igg", "lcg", "always"]

@dataclass(frozen=True)
class RunConfig:
    dataset: DatasetName = "MIT"
    data_dir: Path = Path("data")
    batch_size: int = 1
    max_gn_iter: int = 10
    dx_th: float | None = None
    use_reorder_edges: bool = True
    linear_backend: LinearBackend = "auto"
    cholmod_order: str = "amd"
    selective_alpha: float = 0.3
    gating_mode: GatingMode = "igg"

@dataclass(frozen=True)
class DatasetSpec:
    data_file: Path
    gt_file: Path | None
    dx_th: float
    lc_gap: int
    ent_th: float | None = None

def dataset_spec(cfg: RunConfig) -> DatasetSpec:
    base = cfg.data_dir
    table: dict[str, DatasetSpec] = {
        "MIT": DatasetSpec(base / "MITb_g2o.g2o", base / "MIT_ground_truth.mat", 1e-3, 4, 1.0),
        "Intel": DatasetSpec(base / "INTEL_g2o.g2o", base / "Intel_ground_truth.mat", 1e-6, 4, 0.72),
        "CSAIL": DatasetSpec(base / "CSAIL_P_toro.graph", base / "CSAIL_ground_truth.mat", 1e-5, 4, 0.95),
        "FR079": DatasetSpec(base / "FR079_P_toro.graph", base / "FR079_ground_truth.mat", 1e-4, 4, 0.6),
        "FRH": DatasetSpec(base / "FRH_P_toro.graph", base / "FRH_ground_truth.mat", 1e-7, 1, 0.45),
        "MITP": DatasetSpec(base / "MITP.mat", base / "MITP_ground_truth.mat", 1e-3, 4, 1.0),
    }
    spec = table[cfg.dataset]
    if cfg.dx_th is not None:
        spec = DatasetSpec(spec.data_file, spec.gt_file, cfg.dx_th, spec.lc_gap, spec.ent_th)
    return spec
