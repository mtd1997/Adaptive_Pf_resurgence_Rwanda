# Stage 02 — Alignment and coverage

## Scripted

`align_minimap2.sh` aligns each sample's reads to the P. falciparum 3D7
reference (version 3, PlasmoDB release 68) with minimap2 2.30-r1287, retains
mapped reads only, sorts and indexes with samtools 1.22.1, and summarises
per-amplicon depth with bedtools 2.31.

## Targeted amplicons

Seven loci, 25.1 kb in total. Coordinates are in `data/metadata/amplicons.bed`;
the gene names in column 4 become the columns of `coverage_matrix.tsv`, in
BED order.

| Gene | Chromosome | Size (bp) |
| --- | --- | --- |
| MSP2 | Pf3D7_02_v3 | 3,342 |
| CSP  | Pf3D7_03_v3 | 2,034 |
| DHFR | Pf3D7_04_v3 | 2,168 |
| MDR1 | Pf3D7_05_v3 | 7,140 |
| CRT  | Pf3D7_07_v3 | 3,956 |
| DHPS | Pf3D7_08_v3 | 3,161 |
| K13  | Pf3D7_13_v3 | 3,277 |

## Not scripted — performed through graphical interfaces

Two quality-control steps in this stage were carried out interactively and
therefore have no script. They are recorded here for completeness.

**Qualimap 2.2.2** — BAM quality control. Sorted BAM files from this stage were
assessed for mapping quality distribution, coverage uniformity across the
targeted regions, insert statistics and duplication. The equivalent command
would be:

    qualimap bamqc -bam <sample>.bam -gff ../../data/metadata/amplicons.bed \
        -outdir <outdir> --java-mem-size=4G

**IGV 2.16.0** — visual inspection of the amplicons. Alignments were loaded
against the 3D7 reference to confirm amplicon boundaries, inspect coverage
drop-off at amplicon edges, and visually verify variant calls at the
drug-resistance codons. IGV is interactive by nature and is not scriptable.
