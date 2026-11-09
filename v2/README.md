# OPV-DB v2.0.0 analysis workflow

Input: the extracted, unmodified data archive from
https://doi.org/10.5281/zenodo.23951239.

```sh
python -m pip install -r v2/requirements.txt
python v2/reproduce.py --release /path/to/extracted/data-release --out /path/to/results
```

Run these commands from the repository root. Outputs are written separately.
No network access is required during analysis. Figure 2 uses DOI/year metadata
stored in `assets/figure2_article_metadata.csv`; the data release itself contains
no paper-reproduction files. Figure 1 uses the authored, editable SVG in assets.
Figures 2–4 and Supplementary Figure S1 are calculated from the published tables.

The binary table preserves 37 fields. Its supplement contains 45 descriptive
fields and id/doi (47 columns). Ternary, indoor and material-reference tables
have 68, 57 and 7 columns. Identified binary controls prepared without a third active material
are excluded from the ternary collection. Its literature coverage remains partial.

For each material role, HOMO/LUMO coverage requires both values. Each extra
information group in Figure 3d is combined independently with the four finite
performance values. Field presence is not an accuracy measurement.

Figure 4 uses nominal irradiance 100 mW/cm2; no observed irradiance is imputed.
FF <=1.5 is rescaled to percent only for the diagnostic. Strict performance
selection requires source/material fields, a complete quartet, the documented
physical ranges and relative PCE difference <=2%, followed by exact duplicate-key
selection. The molecular subset additionally requires donor and acceptor SMILES.
The program compares recomputed members and every cell against both published
37-column subsets. This verifies implementation of the selection rules, not
the chemical identity or extraction accuracy of every observation.

Outputs include inventories, field coverage, subset membership, numerical inputs,
figure files and captions, and a verification summary. Prompts and settings
document extraction, source-based completion and structure checks. They exclude
publisher text and credentials and do not constitute deterministic model replay.

The cumulative-curve export includes the record ID and denominator for each
subset, with tied differences ordered by record ID. This supports direct
comparison with the source data supplied with the manuscript.

Extraction templates document their respective processing stages. For reuse of
the released observations, the data dictionary defines the statistical meaning
of the principal PCE and the separate mean and best-value fields.

The 8 November 2026 file correction retains 32,778 binary observations and 27,863 strict-performance records; the strict molecular subset contains 19,299 records. Figure 2 groups recognised aliases using the DOI union and omits unspecified PCBM. The workflow and Figure 1 asset correspond to these corrected files under the same version DOI.
