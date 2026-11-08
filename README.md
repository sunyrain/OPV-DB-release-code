# OPV-DB release code

Processing, validation and figure code for the OPV-DB organic photovoltaic
device performance dataset.

## Version 2.0.0

Download and extract [OPV-DB v2.0.0](https://doi.org/10.5281/zenodo.23951239).
The Zenodo archive is a pure-data release. Code, extraction templates and
recorded model settings are maintained separately in this repository.

```sh
python -m pip install -r v2/requirements.txt
python v2/reproduce.py --release /path/to/extracted/data-release --out /path/to/results
```

The release directory must contain `data/` and `validation/`. The offline
workflow reads the seven published tables, verifies their sizes and linked
identifiers, recomputes both strict binary subsets, compares their members and
all cells with the supplied subsets, and generates the main figures and
specialist-collection coverage figure. It does not modify the release.

The primary binary collection contains 32,778 records; the ternary and indoor
collections contain 2,588 and 1,678. The strict performance and molecular
subsets contain 27,863 and 19,299 records. The binary supplement has 45
descriptive fields plus two linking identifiers (47 columns).

See [workflow details](v2/README.md), [extraction settings](v2/prompts/execution_settings.json)
and the article-text-free templates in `v2/prompts/`. The six historical search
definitions are retained there with a description of their coverage and timestamp
meaning. Article full text, PDFs, source images and credentials are excluded.

## Earlier workflows

The existing `scripts/`, `prompts/` and root requirements are retained for earlier
release formats. They do not target the v2 pure-data layout. Data
[v1.1.0](https://doi.org/10.5281/zenodo.22177960) remains the revision comparison
baseline; [v1.0](https://doi.org/10.5281/zenodo.20841543) is a historical release.
The [OPV2D repository](https://github.com/sunyrain/OPV2D) contains the website.

Code is MIT licensed. The data are CC BY 4.0 with the attribution included in the
Zenodo archive. Model outputs cannot be reproduced deterministically from
prompts alone; the offline workflow reproduces published-table analyses.
