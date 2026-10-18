# Stage 03 — Variant calling

`run_clair3.sh` calls SNPs per sample with Clair3 in haploid mode, producing a
VCF (used for drug-resistance analysis) and a gVCF (used for joint genotyping
in the population-genetics track).

## Clair3 version

The analyses were run with the `hkubal/clair3:latest` image, which at the time
of the analysis (November 2025) was **v2.0.0**. The script pins `v2.0.0`
explicitly so the call set can be reproduced; `latest` now resolves to a
different build.

## Model

`r1041_e82_400bps_sup_v500`, matching the R10.4.1 chemistry and the
super-accurate (SUP) basecalling described in `scripts/00_basecalling/`.

This model is **not** bundled in the Clair3 image. It is distributed through
ONT's rerio repository and must be downloaded first:

    git clone https://github.com/nanoporetech/rerio
    python3 rerio/download_model.py --clair3

Then set `CLAIR3_MODEL_DIR` in `config.sh` to the resulting
`rerio/clair3_models/` folder. The script mounts it into the container.

## Calling parameters

    --haploid_sensitive        retain positions supporting more than one allele
    --include_all_ctgs         call across all contigs
    --gvcf                     emit a gVCF alongside the VCF
    --var_pct_full=1           full-alignment model on all variant candidates
    --ref_pct_full=1           full-alignment model on all reference candidates
    --no_phasing_for_fa        no phasing during full alignment
    --use_whatshap_for_final_output_phasing
    --remove_intermediate_dir

`--haploid_sensitive` rather than `--haploid_precise` is deliberate. P.
falciparum is haploid in the human host, but natural infections are frequently
polyclonal, and mixed genotypes within a host are a subject of this study.
Forcing a single call per position would discard that signal.
