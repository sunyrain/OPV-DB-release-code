# Changelog

## 2.0.0 corrected data and website sync — 2026-11-10

- Re-ran the offline workflow against the current archive at DOI
  10.5281/zenodo.23951239: 32,778 binary observations, 27,863 strict performance
  and 19,299 strict molecular records. Both subsets match the deposited tables
  in membership and every exported cell.
- Include record IDs and subset denominators in cumulative-curve source data;
  order tied differences by record ID for reproducible exports.
- Update code-availability wording and record the additional tested environment.
- The website at https://sunyrain.github.io/OPV2D/ uses the published
  v2.0.0 tables and links to the current DOI.

## 2.0.0 — 2026-11-04

- Target the published pure-data release, DOI 10.5281/zenodo.23951239.
- Recompute current figures and coverage from the 37-column binary table and
  47-column linked supplement, with the final ternary and indoor schemas.
- Verify strict-subset membership and all exported cells against independently
  recomputed selection; preserve legacy scripts for older releases.
- Include extraction templates, recorded settings and historical query definitions
  separately from the data archive.

Initial subset sizes before the data correction: 27,863 strict performance;
19,302 strict molecular. The current corrected molecular subset has 19,299.
