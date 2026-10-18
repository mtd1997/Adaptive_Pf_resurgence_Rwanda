# Metadata

## Files

| File | In git | What it is |
|---|---|---|
| `samples_metadata.tsv` | yes | de-identified sample metadata; every analysis script reads this |
| `samples_metadata_internal.tsv` | **no** | the full table, including exact age and health centre |
| `metadata_preparation_report.txt` | yes | what `prepare_metadata.R` did, and what it found |
| `amplicons.bed` | yes | the seven amplicon regions |
| `resistance_markers.tsv` | yes | which mutations count as resistance markers |

Both metadata files are produced by `scripts/00_metadata/prepare_metadata.R`
from the working spreadsheet. Do not edit them by hand — edit the source and
re-run, so the report stays an accurate record.

## De-identification

`samples_metadata.tsv` is the published file, and two columns are removed
from it.

**Exact age** is replaced by `age_band`. Age, district, sampling year and sex
together can single out one patient in a small district even with no name
attached. Banding age breaks that combination while keeping everything the
analyses need — the bands are the same ones the demographic table in the
supplement reports.

**Health centre** is removed. Health centre plus sampling date plus age
narrows a patient down more sharply than district does.

Nothing in the manuscript depends on either at full resolution. The published
analyses stratify by district and period.

The unbanded table stays in `samples_metadata_internal.tsv`, which
`.gitignore` excludes. If a reviewer or a reader asks for it, it is available
on request under a data-sharing agreement; the Data Availability statement
should say so.

## Columns in `samples_metadata.tsv`

| Column | Type | Description |
|---|---|---|
| `SAMPLE` | text | sample identifier; matches the sample names in the VCF headers |
| `period` | text | `Baseline-Kigali` or `Resurgence` |
| `region` | text | geographic grouping used in the regional comparisons |
| `district` | text | Rwandan district of collection |
| `year` | integer | year of collection |
| `age_band` | text | `<5`, `5-14`, `15-24`, `25-34`, `35-49`, `50+` |
| `sex` | text | as recorded at collection |
| `k13_status` | text | *kelch13* genotype grouping used in the PCoA panels |

`samples_metadata_internal.tsv` has the same columns plus `age` (exact),
`health_centre`, and `period_original`.

## About the period labels

The manuscript uses **Baseline-Kigali** and **Resurgence**. Earlier working
files used **Pre-Resurgence** and **Post-Resurgence** for the same two
groups. `prepare_metadata.R` accepts either and writes the published wording,
keeping the original value in `period_original` in the internal file so no
relabelling is silent.

Both labellings describe the same 257 and 919 samples.

## Sample count

1,176 samples.

The sequencing run produced 1,179. Three were dropped after the metadata
check established they were not of Rwandan origin. That exclusion happened
upstream of this directory and is not revisited by `prepare_metadata.R` —
the script does not drop samples, so a count other than 1,176 here means the
source spreadsheet changed rather than that the script filtered something.

## `amplicons.bed`

The seven targeted loci, 25,078 bp in total:

| Chromosome | Start | End | Gene |
|---|---|---|---|
| Pf3D7_02_v3 | 271575 | 274917 | MSP2 |
| Pf3D7_03_v3 | 221111 | 223145 | CSP |
| Pf3D7_04_v3 | 747897 | 750065 | DHFR |
| Pf3D7_05_v3 | 955955 | 963095 | MDR1 |
| Pf3D7_07_v3 | 402385 | 406341 | CRT |
| Pf3D7_08_v3 | 547896 | 551057 | DHPS |
| Pf3D7_13_v3 | 1724600 | 1727877 | K13 |

BED coordinates are zero-based and half-open; VCF positions are one-based.
Scripts that compare the two account for it.

This file drives the per-gene coverage matrix, the gene assignment in the
annotation step, and the gene ordering in the figures. Columns must be
tab-separated, with the gene name in column 4. Check that with:

    cat -A data/metadata/amplicons.bed | head -2

Tabs show as `^I`. Spaces instead of tabs will cause the gene column to be
read as part of the coordinates.

## `resistance_markers.tsv`

Which amino-acid changes are flagged as resistance-associated, and therefore
drawn as diamonds in the selection-scan figure. The original flagging rule
was inside a script that was not kept, so this list is a documented
reconstruction. `build_selection_scan_annotation.py --compare-to` verifies it
against the annotation table produced by the original analysis.
