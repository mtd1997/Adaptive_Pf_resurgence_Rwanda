# Stage 00 — Basecalling, demultiplexing and trimming

No scripts in this folder. These steps were performed through the MinKNOW
graphical interface rather than the command line, so they are documented here
with the exact settings used, so that the stage can be reproduced.

## Library preparation and sequencing

| Setting | Value |
| --- | --- |
| Library prep | Native Barcoding Kit v14 — SQK-NBD114.96 and SQK-NBD114.24 |
| Supporting reagents | NEBNext Ultra II; NEBNext Quick Ligation Module (NEB, E6056) |
| Device | MinION Mk1B |
| Flow cell | R10.4.1 (FLO-MIN114) |
| Run control software | MinKNOW 25.09.16 |
| Basecaller | Dorado (via MinKNOW) |

Two kit sizes were used depending on the number of samples pooled per run, as
sample collection accumulated over time. A subset of samples was sequenced in
replicate to assess reproducibility.

## Pass 1 — during the run

Basecalling and demultiplexing were performed live during sequencing using the
high-accuracy (HAC) model.

## Pass 2 — after the run

All runs were re-basecalled with the super-accurate (SUP) model, with adapters
and barcodes trimmed in the same operation:

| Parameter | Value |
| --- | --- |
| Model | super-accurate (SUP) |
| Adapter and barcode trimming | enabled |
| Minimum barcode quality score | 60 |
| Minimum read quality score | 9 |
| Minimum read length filter | none applied |

The SUP-re-basecalled, trimmed reads are the input to stage 01.

## Reproducing this stage on the command line

The published data were produced through the MinKNOW interface with the
settings tabulated above. For orientation, the equivalent Dorado commands are:

    dorado basecaller sup <pod5_dir> \
        --kit-name SQK-NBD114.96 \
        --trim adapters \
        --min-qscore 9 \
        > calls.bam

    dorado demux calls.bam \
        --output-dir <demux_dir> \
        --kit-name SQK-NBD114.96 \
        --emit-fastq

Substitute SQK-NBD114.24 for runs prepared with the 24-barcode kit.
