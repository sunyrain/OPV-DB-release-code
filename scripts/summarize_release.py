"""Recompute OPV-DB public release summary tables from archived CSV files."""

from __future__ import annotations

import argparse
import json

import pandas as pd

from common import (
    CORE_METRICS,
    DEVICE_FILES,
    MATERIALS_FILE,
    ReleaseReader,
    ensure_out_dir,
    filled_count,
    non_missing,
    to_numeric,
)


SUMMARY_FIELDS = [
    "doi_norm",
    "donor",
    "acceptor",
    "donor_canonical",
    "acceptor_canonical",
    "donor_smiles",
    "acceptor_smiles",
    "pce",
    "voc",
    "jsc",
    "ff",
    "pce_recomputed",
    "pce_relative_error_percent",
    "device_structure",
    "solvent",
    "additive",
    "active_layer_thickness",
    "annealing_temp",
]


def pct(numerator: int, denominator: int) -> float:
    return round(float(numerator) / float(denominator) * 100.0, 3) if denominator else 0.0


def label_series(frame: pd.DataFrame, canonical: str, raw: str) -> pd.Series:
    if canonical in frame.columns:
        result = frame[canonical].astype("string")
    else:
        result = pd.Series(pd.NA, index=frame.index, dtype="string")
    if raw in frame.columns:
        result = result.fillna(frame[raw].astype("string"))
    return result.str.strip().replace("", pd.NA)


def summarize_tier(name: str, frame: pd.DataFrame) -> dict:
    records = len(frame)
    all_four_mask = pd.Series(True, index=frame.index)
    for column in CORE_METRICS:
        all_four_mask &= non_missing(frame[column])

    pce_error = to_numeric(frame, "pce_relative_error_percent")
    pce_values = to_numeric(frame, "pce")
    donor_label = label_series(frame, "donor_canonical", "donor")
    acceptor_label = label_series(frame, "acceptor_canonical", "acceptor")

    checked = pce_error.notna()
    both_smiles = non_missing(frame["donor_smiles"]) & non_missing(frame["acceptor_smiles"])

    return {
        "tier": name,
        "records": records,
        "unique_dois": int(frame["doi_norm"].dropna().nunique()),
        "unique_donors": int(donor_label.dropna().nunique()),
        "unique_acceptors": int(acceptor_label.dropna().nunique()),
        "unique_donor_acceptor_pairs": int(pd.DataFrame({"d": donor_label, "a": acceptor_label}).dropna().drop_duplicates().shape[0]),
        "all_four_metrics_records": int(all_four_mask.sum()),
        "all_four_metrics_percent": pct(int(all_four_mask.sum()), records),
        "both_smiles_records": int(both_smiles.sum()),
        "both_smiles_percent": pct(int(both_smiles.sum()), records),
        "pce_records": filled_count(frame, "pce"),
        "pce_median": round(float(pce_values.median()), 4) if pce_values.notna().any() else None,
        "pce_mean": round(float(pce_values.mean()), 4) if pce_values.notna().any() else None,
        "pce_max": round(float(pce_values.max()), 4) if pce_values.notna().any() else None,
        "pce_error_checked_records": int(checked.sum()),
        "pce_error_le_2_records": int((pce_error <= 2).sum()),
        "pce_error_le_2_percent": pct(int((pce_error <= 2).sum()), int(checked.sum())),
        "pce_error_le_5_records": int((pce_error <= 5).sum()),
        "pce_error_le_5_percent": pct(int((pce_error <= 5).sum()), int(checked.sum())),
        "pce_error_le_10_records": int((pce_error <= 10).sum()),
        "pce_error_le_10_percent": pct(int((pce_error <= 10).sum()), int(checked.sum())),
        "pce_error_gt_10_records": int((pce_error > 10).sum()),
        "device_structure_records": filled_count(frame, "device_structure"),
        "device_structure_percent": pct(filled_count(frame, "device_structure"), records),
        "solvent_records": filled_count(frame, "solvent"),
        "solvent_percent": pct(filled_count(frame, "solvent"), records),
        "additive_records": filled_count(frame, "additive"),
        "additive_percent": pct(filled_count(frame, "additive"), records),
        "active_layer_thickness_records": filled_count(frame, "active_layer_thickness"),
        "active_layer_thickness_percent": pct(filled_count(frame, "active_layer_thickness"), records),
        "annealing_temp_records": filled_count(frame, "annealing_temp"),
        "annealing_temp_percent": pct(filled_count(frame, "annealing_temp"), records),
    }


def field_coverage(name: str, frame: pd.DataFrame) -> list[dict]:
    rows = []
    for field in SUMMARY_FIELDS:
        if field not in frame.columns:
            continue
        count = filled_count(frame, field)
        rows.append(
            {
                "tier": name,
                "field": field,
                "filled_records": count,
                "filled_percent": pct(count, len(frame)),
            }
        )
    return rows


def top_tables(full: pd.DataFrame, strict: pd.DataFrame) -> dict[str, pd.DataFrame]:
    full = full.copy()
    strict = strict.copy()
    full["donor_label"] = label_series(full, "donor_canonical", "donor")
    full["acceptor_label"] = label_series(full, "acceptor_canonical", "acceptor")
    strict["donor_label"] = label_series(strict, "donor_canonical", "donor")
    strict["acceptor_label"] = label_series(strict, "acceptor_canonical", "acceptor")
    strict["pce_num"] = to_numeric(strict, "pce")

    donors = full.groupby("donor_label", dropna=True).size().reset_index(name="records").sort_values("records", ascending=False)
    acceptors = full.groupby("acceptor_label", dropna=True).size().reset_index(name="records").sort_values("records", ascending=False)
    pairs_counts = full.groupby(["donor_label", "acceptor_label"], dropna=True).size().reset_index(name="records")
    pair_max = strict.groupby(["donor_label", "acceptor_label"], dropna=True)["pce_num"].max().reset_index(name="max_pce")
    pairs = pairs_counts.merge(pair_max, on=["donor_label", "acceptor_label"], how="left").sort_values(["records", "max_pce"], ascending=[False, False])
    return {
        "top_donors.csv": donors.head(100),
        "top_acceptors.csv": acceptors.head(100),
        "top_donor_acceptor_pairs.csv": pairs.head(100),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="Path to Zenodo ZIP or extracted release directory.")
    parser.add_argument("--out", default="outputs/summary", help="Output directory.")
    args = parser.parse_args()

    out_dir = ensure_out_dir(args.out)
    reader = ReleaseReader.open(args.release)
    try:
        frames = {tier: reader.read_csv(path) for tier, path in DEVICE_FILES.items()}
        materials = reader.read_csv(MATERIALS_FILE)
    finally:
        reader.close()

    summary = pd.DataFrame([summarize_tier(tier, frame) for tier, frame in frames.items()])
    coverage = pd.DataFrame([row for tier, frame in frames.items() for row in field_coverage(tier, frame)])
    summary.to_csv(out_dir / "tier_summary_recomputed.csv", index=False)
    coverage.to_csv(out_dir / "field_coverage_recomputed.csv", index=False)

    for filename, table in top_tables(frames["full"], frames["strict_performance_benchmark"]).items():
        table.to_csv(out_dir / filename, index=False)

    material_summary = {
        "records": int(len(materials)),
        "smiles_records": filled_count(materials, "smiles"),
        "unique_names": int(materials["name"].dropna().nunique()) if "name" in materials.columns else 0,
    }
    (out_dir / "material_summary.json").write_text(json.dumps(material_summary, indent=2), encoding="utf-8")

    key_metrics = {
        "tier_summary": summary.to_dict(orient="records"),
        "material_summary": material_summary,
    }
    (out_dir / "key_metrics.json").write_text(json.dumps(key_metrics, indent=2), encoding="utf-8")
    print(f"Wrote summary outputs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

