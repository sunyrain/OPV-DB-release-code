"""Validate the public OPV-DB Zenodo release package."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from common import (
    CORE_METRICS,
    COVERAGE_FILE,
    DEVICE_FILES,
    EXPECTED_COUNTS,
    FORBIDDEN_PUBLIC_COLUMNS,
    MANIFEST_FILE,
    MATERIALS_FILE,
    PCE_FILE,
    RULES_FILE,
    SHA256SUMS_FILE,
    STRICT_MOLECULAR_REQUIRED_FIELDS,
    STRICT_REQUIRED_FIELDS,
    SUMMARY_FILE,
    ReleaseReader,
    ensure_out_dir,
    filled_count,
    non_missing,
    normalized_ids,
    parse_sha256sums,
    sha256_bytes,
    to_numeric,
)


REQUIRED_FILES = [
    "CITATION.cff",
    "DATA_DICTIONARY.md",
    "LICENSE",
    "README.md",
    "THIRD_PARTY_ATTRIBUTION.md",
    MATERIALS_FILE,
    *DEVICE_FILES.values(),
    MANIFEST_FILE,
    SHA256SUMS_FILE,
    PCE_FILE,
    COVERAGE_FILE,
    RULES_FILE,
    SUMMARY_FILE,
]


def add_check(checks: list[dict], name: str, passed: bool, detail: str = "") -> None:
    checks.append({"check": name, "passed": bool(passed), "detail": detail})


def pct(numerator: int, denominator: int) -> float:
    return float(numerator) / float(denominator) * 100.0 if denominator else 0.0


def compare_float(actual: float, expected: float, tolerance: float = 0.05) -> bool:
    if pd.isna(actual) and pd.isna(expected):
        return True
    return abs(float(actual) - float(expected)) <= tolerance


def validate_checksums(reader: ReleaseReader, checks: list[dict]) -> None:
    manifest = reader.read_json(MANIFEST_FILE)
    sums = parse_sha256sums(reader.read_text(SHA256SUMS_FILE))

    manifest_paths = {item["path"]: item for item in manifest.get("files", [])}
    add_check(checks, "manifest_has_files", bool(manifest_paths), f"{len(manifest_paths)} files")
    add_check(checks, "sha256sums_has_entries", bool(sums), f"{len(sums)} entries")

    missing_from_zip = sorted(path for path in manifest_paths if not reader.exists(path))
    add_check(checks, "manifest_files_exist", not missing_from_zip, "; ".join(missing_from_zip[:5]))

    mismatches = []
    for path, item in manifest_paths.items():
        if not reader.exists(path):
            continue
        digest = sha256_bytes(reader.read_bytes(path))
        if digest != item.get("sha256"):
            mismatches.append(f"manifest:{path}")
        if path in sums and digest != sums[path]:
            mismatches.append(f"SHA256SUMS:{path}")
    add_check(checks, "checksums_match", not mismatches, "; ".join(mismatches[:10]))

    sums_not_in_manifest = sorted(set(sums) - set(manifest_paths))
    add_check(checks, "sha256sums_subset_of_manifest", not sums_not_in_manifest, "; ".join(sums_not_in_manifest[:5]))


def validate_device_tiers(reader: ReleaseReader, checks: list[dict]) -> dict[str, pd.DataFrame]:
    frames = {tier: reader.read_csv(path) for tier, path in DEVICE_FILES.items()}
    summary = reader.read_csv(SUMMARY_FILE).set_index("tier")
    coverage = reader.read_csv(COVERAGE_FILE)
    pce_records = reader.read_csv(PCE_FILE)
    materials = reader.read_csv(MATERIALS_FILE)

    for tier, frame in frames.items():
        expected = EXPECTED_COUNTS[tier]
        add_check(checks, f"{tier}_row_count", len(frame) == expected["records"], f"{len(frame)}")
        doi_count = int(frame["doi_norm"].dropna().nunique()) if "doi_norm" in frame.columns else 0
        add_check(checks, f"{tier}_doi_count", doi_count == expected["unique_dois"], f"{doi_count}")
        add_check(checks, f"{tier}_duplicate_ids", frame["id"].duplicated().sum() == 0, f"{int(frame['id'].duplicated().sum())}")
        forbidden = sorted(FORBIDDEN_PUBLIC_COLUMNS.intersection(frame.columns))
        add_check(checks, f"{tier}_no_internal_columns", not forbidden, ", ".join(forbidden))

        if tier != "full":
            required = STRICT_MOLECULAR_REQUIRED_FIELDS if tier == "strict_molecular_benchmark" else STRICT_REQUIRED_FIELDS
            missing_required = [column for column in required if column not in frame.columns or not bool(non_missing(frame[column]).all())]
            add_check(checks, f"{tier}_required_fields_complete", not missing_required, ", ".join(missing_required))

            rel_error = to_numeric(frame, "pce_relative_error_percent")
            add_check(checks, f"{tier}_pce_error_le_2", bool((rel_error <= 2.0).all()), f"max={rel_error.max()}")

            ranges = {
                "pce": (0, 30),
                "voc": (0, 2),
                "jsc": (0, 50),
                "ff": (0, 100),
            }
            range_issues = []
            for column, (lo, hi) in ranges.items():
                values = to_numeric(frame, column)
                bad = values.notna() & ~values.between(lo, hi)
                if bad.any():
                    range_issues.append(f"{column}:{int(bad.sum())}")
            add_check(checks, f"{tier}_physical_ranges", not range_issues, "; ".join(range_issues))

        pce_calc = to_numeric(frame, "voc") * to_numeric(frame, "jsc") * to_numeric(frame, "ff") / 100.0
        pce_stored = to_numeric(frame, "pce_recomputed")
        comparable = pce_calc.notna() & pce_stored.notna()
        max_delta = float((pce_calc[comparable] - pce_stored[comparable]).abs().max()) if comparable.any() else 0.0
        add_check(checks, f"{tier}_pce_recomputed_formula", max_delta <= 1e-6, f"max_delta={max_delta:.3g}")

        if tier in summary.index:
            row = summary.loc[tier]
            add_check(checks, f"{tier}_summary_records", int(row["records"]) == len(frame), f"{row['records']} vs {len(frame)}")
            add_check(checks, f"{tier}_summary_unique_dois", int(row["unique_dois"]) == doi_count, f"{row['unique_dois']} vs {doi_count}")
            for column in ["pce", "voc", "jsc", "ff", "donor", "acceptor", "donor_smiles", "acceptor_smiles"]:
                cov = coverage[(coverage["tier"] == tier) & (coverage["field"] == column)]
                if cov.empty:
                    continue
                actual = filled_count(frame, column)
                expected_filled = int(cov.iloc[0]["filled_records"])
                expected_pct = float(cov.iloc[0]["filled_percent"])
                add_check(checks, f"{tier}_{column}_coverage_count", actual == expected_filled, f"{actual} vs {expected_filled}")
                add_check(checks, f"{tier}_{column}_coverage_percent", compare_float(pct(actual, len(frame)), expected_pct), f"{pct(actual, len(frame)):.3f} vs {expected_pct:.3f}")

    full_ids = normalized_ids(frames["full"])
    perf_ids = normalized_ids(frames["strict_performance_benchmark"])
    mol_ids = normalized_ids(frames["strict_molecular_benchmark"])
    add_check(checks, "strict_performance_subset_of_full", perf_ids.issubset(full_ids), f"missing={len(perf_ids - full_ids)}")
    add_check(checks, "strict_molecular_subset_of_performance", mol_ids.issubset(perf_ids), f"missing={len(mol_ids - perf_ids)}")

    checkable_mask = pd.Series(True, index=frames["full"].index)
    for column in [*CORE_METRICS, "pce_recomputed", "pce_relative_error_percent"]:
        checkable_mask &= non_missing(frames["full"][column])
    checkable_ids = set(frames["full"].loc[checkable_mask, "id"].astype(str))
    pce_ids = set(pce_records["id"].dropna().astype(str))
    add_check(checks, "pce_consistency_row_count", len(pce_records) == EXPECTED_COUNTS["pce_consistency_records"]["records"], f"{len(pce_records)}")
    add_check(checks, "pce_consistency_matches_full_checkable_ids", pce_ids == checkable_ids, f"missing={len(checkable_ids - pce_ids)}, extra={len(pce_ids - checkable_ids)}")

    add_check(checks, "materials_row_count", len(materials) == EXPECTED_COUNTS["materials_reference"]["records"], f"{len(materials)}")
    add_check(checks, "materials_no_source_column", "source" not in materials.columns, "")
    add_check(checks, "materials_duplicate_names_allowed_documented", "name" in materials.columns, "reference table is supporting, not authoritative")

    return frames


def write_reports(checks: list[dict], out_dir: Path) -> bool:
    report = {
        "passed": all(item["passed"] for item in checks),
        "total_checks": len(checks),
        "failed_checks": [item for item in checks if not item["passed"]],
        "checks": checks,
    }
    (out_dir / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = ["# OPV-DB Release Validation", ""]
    lines.append(f"Overall result: {'PASS' if report['passed'] else 'FAIL'}")
    lines.append("")
    lines.append("| Check | Result | Detail |")
    lines.append("| --- | --- | --- |")
    for item in checks:
        result = "PASS" if item["passed"] else "FAIL"
        detail = str(item.get("detail", "")).replace("|", "\\|")
        lines.append(f"| `{item['check']}` | {result} | {detail} |")
    (out_dir / "validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return bool(report["passed"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="Path to Zenodo ZIP or extracted release directory.")
    parser.add_argument("--out", default="outputs/validation", help="Output directory for validation reports.")
    args = parser.parse_args()

    checks: list[dict] = []
    out_dir = ensure_out_dir(args.out)
    reader = ReleaseReader.open(args.release)
    try:
        available = set(reader.list_files())
        missing = sorted(path for path in REQUIRED_FILES if path not in available)
        add_check(checks, "required_files_present", not missing, "; ".join(missing))
        if not missing:
            validate_checksums(reader, checks)
            validate_device_tiers(reader, checks)
    finally:
        reader.close()

    passed = write_reports(checks, out_dir)
    print(f"Validation {'PASS' if passed else 'FAIL'}: {out_dir / 'validation_report.md'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())

