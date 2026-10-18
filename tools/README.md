# Software

Every tool used in the analysis, the version used, and how to obtain it.

## Why the binaries are not here

GATK, SnpEff and the SnpEff genome database together come to several hundred
megabytes. This repository is archived to Zenodo and given a DOI, so
everything committed is preserved permanently — which is the right treatment
for code and the wrong treatment for redistributable third-party binaries.

What a reviewer needs in order to reproduce the analysis is the exact version
of each tool, where to get it, and the recipe for anything built locally. All
three are here. The SnpEff genome database in particular is not redistributed
but is fully rebuildable from `snpeff/build_pf3d7_database.sh`.

Install the tools wherever you like and point `config.sh` at them. The
pipeline reads every tool path from there, so nothing is tied to one machine.

## Versions used

`../environment/software_versions.tsv` is the authoritative list, and it is
generated from the machine the analysis ran on rather than written by hand —
see `../environment/capture_versions.sh`. That table is also what the Methods
software list should be checked against.

## A discrepancy worth knowing about

The gVCF headers written during variant calling record
`clair3_version=1.0.11`, while the Docker image used was tagged `v2.0.0`.

The explanation is that Clair3's Docker images carried the package's internal
version string, which lagged behind the image tag. The image actually used
was `hkubal/clair3:v2.0.0`.

This is recorded here rather than quietly reconciled, because a reader who
opens a gVCF will see `1.0.11` and should not have to wonder which is right.
The Methods should cite the image tag, since that is what pins the software.

## Container images

Clair3 was run from a pinned Docker image rather than a local installation:

    hkubal/clair3:v2.0.0

Pinning the tag matters. `hkubal/clair3:latest` pointed at v2.0.0 at the time
of the analysis and has moved since, so `latest` would no longer reproduce
these results.

### The basecalling model

Clair3 needs a model matched to the basecalling settings. The model used was:

    r1041_e82_400bps_sup_v500

This corresponds to R10.4.1 flow cells, 400 bps translocation speed, and
**super-accurate (SUP)** basecalling. It is not bundled with the Docker image
and was obtained from ONT's Rerio repository:

    https://github.com/nanoporetech/rerio

    python3 download_model.py --clair3

The models land in `rerio/clair3_models/`. Point `CLAIR3_MODEL_DIR` in
`config.sh` at that directory.

Using a model that does not match the basecalling accuracy changes the
variant calls. An earlier exploratory run in this project used a HAC model by
mistake; that run is not part of the published analysis, and the SUP model
named above is what produced the results in the paper.

## Tools used through a graphical interface

Three steps were run interactively and have no script:

| Tool | Used for |
|---|---|
| MinKNOW | sequencing run control and live basecalling |
| Qualimap | alignment quality inspection |
| IGV | manual inspection of variant calls |

`../scripts/00_basecalling/README.md` documents the basecalling settings in
full, including both passes and the barcoding kits. Qualimap and IGV were
used to look at the data, not to generate any published result, so no output
of theirs feeds a figure.
