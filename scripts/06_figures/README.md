# Stage 06 — Manuscript figures

This stage turns the processed variant data into the figures and tables that
appear in the paper and its supplement.

Scripts are named by **published figure number**, not by the order they were
written. Several figures were renumbered during revision, so an older script
name is not a reliable guide to which figure it produced.

---

## What is in this stage

| Script | Produces |
|---|---|
| `figure05_selection_scan_statistics.R` | all selection statistics (no graphics) |
| `figure05_selection_scan_panels.R` | the three Figure 5 panels |

Further figure scripts are added here as they are brought in from the
analysis working directory.

---

## Which figures came from a script, and which did not

Being explicit about this matters more than making the repository look
complete. Three of the published figures were finished by hand, and a reader
who assumed otherwise would waste time trying to regenerate them.

| Figure | How it was produced |
|---|---|
| Figure 1 — sampling maps | script |
| Figure 2 — resistance marker trends | script |
| Figure 3 — *kelch13* by district | script |
| Figure 4 — PCoA of multilocus variation | script |
| **Figure 5 — selection scan** | **panels from a script, composed by hand in Adobe Illustrator** |
| **Figure 6 — Sanger chromatograms** | **prepared by hand in Adobe Illustrator; no script** |
| **Supplementary Figure S2 — *pfmdr1* Sanger** | **prepared by hand; supplied by the senior author, no script** |

For the three manual figures, what the repository provides is the underlying
data and, where one exists, the script that generated the components. The
final visual arrangement is not reproducible from code, and no script here
claims to reproduce it.

---

## Figure 5 in detail

### Running it

Two steps, in order, from the top of the repository:

    Rscript scripts/06_figures/figure05_selection_scan_statistics.R
    Rscript scripts/06_figures/figure05_selection_scan_panels.R

The first computes everything and writes
`results/tables/figure05/selection_statistics_all_variants.tsv`.
The second reads that table and draws it, computing nothing.

Keeping the split means the statistics can be checked without rendering
anything, and the figure can be restyled without any risk of the numbers
moving.

### Before either step

The statistics script needs an annotation table. Build it with:

    python3 scripts/05_downstream_analysis/build_selection_scan_annotation.py \
        --vcf     data/processed/selection_scan_annotated.vcf.gz \
        --bed     data/metadata/amplicons.bed \
        --markers data/metadata/resistance_markers.tsv \
        --out     results/tables/annotated_variants_snpeff_final.tsv

That script replaces a chain of three annotation scripts from the original
analysis that were not kept. Its header explains what was reconstructed and
how to verify it with `--compare-to`. **Run that verification before
regenerating the figure**: the resistance flags decide which points are drawn
as diamonds, so a disagreement there would change the figure.

### How the published figure was assembled

The panels could not be arranged satisfactorily in ggplot2. Two problems
resisted every attempt:

1. ggplot2's automatic legends would not produce coloured gene symbols, two
   distinct point shapes, and a visibly dashed threshold key at the same
   time.
2. Adding a legend to one panel changed that panel's plotting width, so the
   three panels no longer lined up.

Several variants of the figure script were written to work around this. Two
of the resulting SVGs were passed to a co-author, who composed the final
figure in Adobe Illustrator:

- one had the legend correct,
- the other had the panel proportions correct.

The surviving script, published here, is the one that solved both problems in
ggplot2 — a hand-built legend drawn as its own plot, placed in a fixed
right-hand column that is present in all three rows. Section 10 of that
script explains the approach.

The script behind the second SVG was edited in place during that process and
not kept. It differed only in font sizes; it computed nothing that this
script does not compute. Since the final figure was composed manually in any
case, no script reproduces the published Figure 5, and its absence changes
nothing about what is reproducible.

### What the figure draws, and what it does not claim

The dashed lines are **empirical 95th percentiles** of the observed
distribution of each statistic across the analysed variants. A point above
the line is in the top 5% of what was observed. This is a ranking, not a
significance test — there is no null model and no p-value. The figure legend
in the manuscript should describe it in those terms.

`cumulative_pos` in the plotting table is a **drawing coordinate**, not a
genomic one. The seven amplicons lie on different chromosomes, so each gene
is given its own block on the x axis with a fixed gap between blocks.
Positions within a gene keep their true spacing. The real coordinates are in
the `chrom` and `pos` columns.

CRT is absent from the panels because no CRT variant survived filtering into
the selection-scan VCF. `figure05_gene_accounting.tsv` records that as an
explicit zero rather than leaving it to be inferred from a gap.

---

## Source Data

Nature Communications requires a Source Data file with the values behind each
figure panel. Each script in this stage writes a table containing every value
it draws. For Figure 5 those are:

| File | Covers |
|---|---|
| `results/figures/figure05/figure05_plotting_table.tsv` | all three panels |
| `results/figures/figure05/figure05_labelled_variants.tsv` | the labelled points only |
| `results/figures/figure05/figure05_gene_accounting.tsv` | per-gene counts |

These are assembled into the submitted Source Data workbook, one sheet per
panel, named by published panel letter.

For Figures 6 and S2 the Source Data is the underlying measurements and
sequence traces, not a plotting table, since no plotting step took place.

---

## A note on reproducibility claims

Manual figure assembly is ordinary practice and reviewers accept it. What
they do not accept is a repository that appears to regenerate a figure it
cannot. Every script here states at the top whether it reproduces a published
figure or only its components, and the table above is the single place to
check.
