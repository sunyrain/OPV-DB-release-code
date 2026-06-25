"""Generate lightweight figures from the public OPV-DB release files."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import DEVICE_FILES, ReleaseReader, ensure_out_dir, filled_count, non_missing, to_numeric


COLORS = {
    "blue": "#1f5a8f",
    "teal": "#16807c",
    "green": "#3f8f5c",
    "purple": "#7a61a8",
    "orange": "#d98324",
    "red": "#b44a43",
    "gray": "#7e8790",
    "dark": "#17212b",
    "grid": "#d8dee6",
}


def style_axes(ax) -> None:
    ax.grid(color=COLORS["grid"], alpha=0.55)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.set_axisbelow(True)


def save(fig, out_dir: Path, name: str) -> None:
    fig.savefig(out_dir / f"{name}.png", dpi=220, bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def label_series(frame: pd.DataFrame, canonical: str, raw: str) -> pd.Series:
    result = frame[canonical].astype("string") if canonical in frame.columns else pd.Series(pd.NA, index=frame.index, dtype="string")
    if raw in frame.columns:
        result = result.fillna(frame[raw].astype("string"))
    return result.str.strip().replace("", pd.NA)


def acceptor_family(label: str) -> str:
    text = str(label).upper()
    fullerene_tokens = ["PCBM", "PC61", "PC71", "C60", "C70", "ICBA", "FULLERENE"]
    return "Fullerene" if any(token in text for token in fullerene_tokens) else "Non-fullerene"


def draw_field_coverage(frames: dict[str, pd.DataFrame], out_dir: Path) -> None:
    fields = [
        "pce",
        "voc",
        "jsc",
        "ff",
        "donor",
        "acceptor",
        "donor_smiles",
        "acceptor_smiles",
        "device_structure",
        "solvent",
        "additive",
        "active_layer_thickness",
        "annealing_temp",
    ]
    tiers = [
        ("full", "Full"),
        ("strict_performance_benchmark", "Strict performance"),
        ("strict_molecular_benchmark", "Strict molecular"),
    ]
    data = []
    for tier_key, tier_label in tiers:
        frame = frames[tier_key]
        for field in fields:
            data.append({"tier": tier_label, "field": field, "percent": filled_count(frame, field) / len(frame) * 100.0})
    coverage = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    x = np.arange(len(fields))
    width = 0.26
    for i, (_, tier_label) in enumerate(tiers):
        y = coverage.loc[coverage["tier"] == tier_label, "percent"].to_numpy()
        ax.bar(x + (i - 1) * width, y, width=width, label=tier_label)
    ax.set_xticks(x)
    ax.set_xticklabels(fields, rotation=40, ha="right")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Filled records (%)")
    ax.set_title("Field coverage by release layer", loc="left", weight="bold")
    ax.legend(frameon=False, ncol=3)
    style_axes(ax)
    save(fig, out_dir, "field_coverage")


def draw_performance_distributions(strict: pd.DataFrame, out_dir: Path) -> None:
    data = strict.copy()
    for column in ["pce", "voc", "jsc", "ff", "pce_relative_error_percent"]:
        data[column] = to_numeric(data, column)

    fig, axes = plt.subplots(2, 3, figsize=(12.2, 7.4))
    metrics = [
        ("pce", "PCE (%)", COLORS["blue"]),
        ("voc", "Voc (V)", COLORS["teal"]),
        ("jsc", "Jsc (mA cm$^{-2}$)", COLORS["green"]),
        ("ff", "FF (%)", COLORS["purple"]),
    ]
    for ax, (column, label, color) in zip(axes.flat[:4], metrics):
        values = data[column].dropna()
        ax.hist(values, bins=45, color=color, alpha=0.82)
        ax.set_xlabel(label)
        ax.set_ylabel("Records")
        ax.set_title(label, loc="left", weight="bold")
        style_axes(ax)

    ax = axes.flat[4]
    ax.scatter(data["jsc"], data["pce"], s=6, alpha=0.18, color=COLORS["orange"], edgecolors="none")
    ax.set_xlabel("Jsc (mA cm$^{-2}$)")
    ax.set_ylabel("PCE (%)")
    ax.set_title("Jsc-PCE coupling", loc="left", weight="bold")
    style_axes(ax)

    ax = axes.flat[5]
    err = data["pce_relative_error_percent"].dropna()
    ax.hist(err, bins=35, color=COLORS["red"], alpha=0.82)
    ax.set_xlabel("PCE relative error (%)")
    ax.set_ylabel("Records")
    ax.set_title("PCE recomputation error", loc="left", weight="bold")
    style_axes(ax)

    fig.suptitle("Strict performance benchmark", x=0.02, ha="left", weight="bold")
    fig.tight_layout()
    save(fig, out_dir, "performance_distributions")


def draw_material_landscape(full: pd.DataFrame, strict: pd.DataFrame, out_dir: Path) -> None:
    full = full.copy()
    strict = strict.copy()
    full["donor_label"] = label_series(full, "donor_canonical", "donor")
    full["acceptor_label"] = label_series(full, "acceptor_canonical", "acceptor")
    strict["donor_label"] = label_series(strict, "donor_canonical", "donor")
    strict["acceptor_label"] = label_series(strict, "acceptor_canonical", "acceptor")
    strict["pce_num"] = to_numeric(strict, "pce")
    strict["acceptor_family"] = strict["acceptor_label"].map(acceptor_family)

    donors = full.groupby("donor_label", dropna=True).size().sort_values(ascending=False).head(12).sort_values()
    acceptors = full.groupby("acceptor_label", dropna=True).size().sort_values(ascending=False).head(12).sort_values()
    pair_counts = full.groupby(["donor_label", "acceptor_label"], dropna=True).size().reset_index(name="records")
    pair_max = strict.groupby(["donor_label", "acceptor_label"], dropna=True)["pce_num"].max().reset_index(name="max_pce")
    pairs = pair_counts.merge(pair_max, on=["donor_label", "acceptor_label"], how="inner")
    pairs = pairs.sort_values(["records", "max_pce"], ascending=[False, False]).head(10).sort_values("records")
    pairs["label"] = pairs["donor_label"].astype(str) + " / " + pairs["acceptor_label"].astype(str)

    fig = plt.figure(figsize=(12.2, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.48)

    ax = fig.add_subplot(gs[0, 0])
    ax.barh(donors.index.astype(str), donors.values, color=COLORS["blue"], alpha=0.86)
    ax.set_xlabel("Records")
    ax.set_title("Frequent donor labels", loc="left", weight="bold")
    style_axes(ax)

    ax = fig.add_subplot(gs[0, 1])
    ax.barh(acceptors.index.astype(str), acceptors.values, color=COLORS["teal"], alpha=0.86)
    ax.set_xlabel("Records")
    ax.set_title("Frequent acceptor labels", loc="left", weight="bold")
    style_axes(ax)

    ax = fig.add_subplot(gs[1, 0])
    families = ["Fullerene", "Non-fullerene"]
    values = [strict.loc[strict["acceptor_family"].eq(family), "pce_num"].dropna() for family in families]
    parts = ax.violinplot(values, showmedians=True, widths=0.75)
    for body, color in zip(parts["bodies"], [COLORS["purple"], COLORS["green"]]):
        body.set_facecolor(color)
        body.set_alpha(0.58)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(families)
    ax.set_ylabel("PCE (%)")
    ax.set_title("Strict PCE by acceptor family", loc="left", weight="bold")
    style_axes(ax)

    ax = fig.add_subplot(gs[1, 1])
    ax.barh(pairs["label"], pairs["records"], color=COLORS["orange"], alpha=0.78)
    for y, (_, row) in enumerate(pairs.iterrows()):
        ax.text(row["records"] * 1.01, y, f"max {row['max_pce']:.2f}%", va="center", fontsize=8)
    ax.set_xlabel("Full archive records")
    ax.set_title("Top donor-acceptor pairs", loc="left", weight="bold")
    style_axes(ax)

    save(fig, out_dir, "material_landscape")


def draw_validation_summary(frames: dict[str, pd.DataFrame], out_dir: Path) -> None:
    full = frames["full"]
    strict = frames["strict_performance_benchmark"]
    mol = frames["strict_molecular_benchmark"]

    pce_error = to_numeric(full, "pce_relative_error_percent")
    thresholds = [2, 5, 10, 50]
    counts = [int((pce_error <= threshold).sum()) for threshold in thresholds]

    funnel_labels = ["Full", "Core metrics", "Strict perf.", "Strict molecular"]
    core_mask = pd.Series(True, index=full.index)
    for column in ["pce", "voc", "jsc", "ff"]:
        core_mask &= non_missing(full[column])
    funnel_counts = [len(full), int(core_mask.sum()), len(strict), len(mol)]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    ax = axes[0]
    ax.bar([f"<={t}%" for t in thresholds], counts, color=COLORS["blue"], alpha=0.85)
    ax.set_ylabel("Full-release records")
    ax.set_title("PCE consistency thresholds", loc="left", weight="bold")
    style_axes(ax)

    ax = axes[1]
    ax.bar(funnel_labels, funnel_counts, color=[COLORS["gray"], COLORS["teal"], COLORS["green"], COLORS["purple"]], alpha=0.86)
    ax.set_ylabel("Records")
    ax.set_title("Release-layer funnel", loc="left", weight="bold")
    for tick in ax.get_xticklabels():
        tick.set_rotation(18)
        tick.set_ha("right")
    style_axes(ax)

    fig.tight_layout()
    save(fig, out_dir, "validation_summary")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", required=True, help="Path to Zenodo ZIP or extracted release directory.")
    parser.add_argument("--out", default="outputs/figures", help="Output directory.")
    args = parser.parse_args()

    out_dir = ensure_out_dir(args.out)
    reader = ReleaseReader.open(args.release)
    try:
        frames = {tier: reader.read_csv(path) for tier, path in DEVICE_FILES.items()}
    finally:
        reader.close()

    draw_field_coverage(frames, out_dir)
    draw_performance_distributions(frames["strict_performance_benchmark"], out_dir)
    draw_material_landscape(frames["full"], frames["strict_performance_benchmark"], out_dir)
    draw_validation_summary(frames, out_dir)
    print(f"Wrote figures to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

