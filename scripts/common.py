"""Shared helpers for OPV-DB public release scripts."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
import zipfile

import pandas as pd


DEVICE_FILES = {
    "full": "data/opv_devices_full.csv",
    "strict_performance_benchmark": "data/opv_devices_strict_performance_benchmark.csv",
    "strict_molecular_benchmark": "data/opv_devices_strict_molecular_benchmark.csv",
}

MATERIALS_FILE = "data/materials_reference.csv"
SUMMARY_FILE = "validation/opv_quality_tier_summary.csv"
COVERAGE_FILE = "validation/opv_quality_tier_field_coverage.csv"
PCE_FILE = "validation/opv_pce_consistency_records.csv"
RULES_FILE = "validation/opv_quality_tier_rules.json"
MANIFEST_FILE = "metadata/release_manifest.json"
SHA256SUMS_FILE = "metadata/SHA256SUMS.txt"

EXPECTED_COUNTS = {
    "full": {"records": 38945, "unique_dois": 7886},
    "strict_performance_benchmark": {"records": 31403, "unique_dois": 7191},
    "strict_molecular_benchmark": {"records": 21763, "unique_dois": 5465},
    "materials_reference": {"records": 4548},
    "pce_consistency_records": {"records": 37511},
}

FORBIDDEN_PUBLIC_COLUMNS = {
    "source_table",
    "data_source_type",
    "is_baseline_ref",
    "data_quality_flag",
    "verified",
}

STRICT_REQUIRED_FIELDS = ["doi_norm", "donor", "acceptor", "pce", "voc", "jsc", "ff"]
STRICT_MOLECULAR_REQUIRED_FIELDS = STRICT_REQUIRED_FIELDS + ["donor_smiles", "acceptor_smiles"]
CORE_METRICS = ["voc", "jsc", "ff", "pce"]


@dataclass
class ReleaseReader:
    """Read files from a Zenodo ZIP or an extracted release directory."""

    root: Path
    zip_handle: zipfile.ZipFile | None = None

    @classmethod
    def open(cls, release: str | Path) -> "ReleaseReader":
        path = Path(release)
        if path.is_file():
            return cls(root=path, zip_handle=zipfile.ZipFile(path))
        if path.is_dir():
            return cls(root=path)
        raise FileNotFoundError(f"Release path does not exist: {path}")

    def close(self) -> None:
        if self.zip_handle is not None:
            self.zip_handle.close()

    def list_files(self) -> list[str]:
        if self.zip_handle is not None:
            return sorted(name for name in self.zip_handle.namelist() if not name.endswith("/"))
        return sorted(str(path.relative_to(self.root)).replace("\\", "/") for path in self.root.rglob("*") if path.is_file())

    def exists(self, path: str) -> bool:
        if self.zip_handle is not None:
            return path in self.zip_handle.namelist()
        return (self.root / path).is_file()

    def read_bytes(self, path: str) -> bytes:
        if self.zip_handle is not None:
            return self.zip_handle.read(path)
        return (self.root / path).read_bytes()

    def read_text(self, path: str) -> str:
        return self.read_bytes(path).decode("utf-8")

    def read_json(self, path: str) -> dict:
        return json.loads(self.read_text(path))

    def read_csv(self, path: str, **kwargs) -> pd.DataFrame:
        data = self.read_bytes(path)
        return pd.read_csv(BytesIO(data), **kwargs)


def ensure_out_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def parse_sha256sums(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        digest, path = parts
        result[path.strip()] = digest.strip()
    return result


def non_missing(series: pd.Series) -> pd.Series:
    if series.dtype.kind in "biufc":
        return series.notna()
    text = series.astype("string").str.strip()
    return series.notna() & text.notna() & text.ne("")


def filled_count(frame: pd.DataFrame, column: str) -> int:
    if column not in frame.columns:
        return 0
    return int(non_missing(frame[column]).sum())


def to_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(pd.NA, index=frame.index, dtype="Float64")
    return pd.to_numeric(frame[column], errors="coerce")


def normalized_ids(frame: pd.DataFrame) -> set[str]:
    if "id" not in frame.columns:
        return set()
    return set(frame["id"].dropna().astype(str))


def tier_display_name(tier: str) -> str:
    return {
        "full": "full",
        "strict_performance_benchmark": "strict performance",
        "strict_molecular_benchmark": "strict molecular",
    }.get(tier, tier)

