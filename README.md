# OPV-DB Release Code

This repository contains a minimal, redistribution-safe code package for the
OPV-DB Scientific Data release.

It is designed to work only from the public Zenodo data package:

https://doi.org/10.5281/zenodo.20841543

The scripts validate the archived CSV release, recompute summary tables, and
generate lightweight figures from the public files. They do not include
literature harvesting code, publisher full text, article PDFs, LLM extraction
prompts, agent workflows, API calls, or local intermediate databases.

## Recommended hosting

Use a standalone GitHub repository, for example `OPV-DB-release-code`. This is
cleaner for a paper Code Availability statement than mixing validation scripts
into the OPV2D website repository. The OPV2D repository can remain the static
viewer and website.

## Contents

- `scripts/common.py`: release ZIP/directory reader and shared helpers.
- `scripts/validate_release.py`: structural, checksum, tier, and PCE checks.
- `scripts/summarize_release.py`: recomputed tier summaries and field coverage.
- `scripts/make_figures.py`: publication-supporting overview figures.
- `requirements.txt`: minimal Python dependencies.

## Quick start

Create a Python environment and install dependencies:

```bash
python -m venv .venv
.venv/Scripts/activate
pip install -r requirements.txt
```

Run all checks against the Zenodo ZIP:

```bash
python scripts/validate_release.py --release opvdb_scidata_release_v1_20260625.zip --out outputs/validation
```

Recompute summary tables:

```bash
python scripts/summarize_release.py --release opvdb_scidata_release_v1_20260625.zip --out outputs/summary
```

Generate figures:

```bash
python scripts/make_figures.py --release opvdb_scidata_release_v1_20260625.zip --out outputs/figures
```

The `--release` argument can also point to an extracted release directory with
the same internal layout as the Zenodo ZIP.

## Scope

These scripts are for public-release verification and reuse. They intentionally
do not attempt to regenerate OPV-DB from raw literature assets. The public
release is a secondary literature-mined data package; the authoritative
redistributable inputs are the archived CSV files, validation tables, tier
rules, manifest, checksums, data dictionary, license, and attribution files.

