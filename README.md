Introduction
RMAnno is an R package specifically designed for the annotation of RNA-level modification sites. It addresses the inaccuracies, information loss, and integration difficulties that arise when existing tools—mostly derived from DNA-level annotation methods—are applied to RNA modification data. RMAnno provides researchers with an accurate, flexible, and reproducible workflow for annotating RNA modifications.

Background and Motivation
RNA modification analysis typically involves two phases: upstream (quality control, alignment) and downstream (functional annotation, pathway analysis) (Figure 4). Among these, genomic annotation is a fundamental step in downstream analysis, directly affecting the reliability of subsequent enrichment analyses such as GO and KEGG. However, widely used annotation tools suffer from several limitations:

Misaligned tool design: Tools like ChIPseeker and RNAmod were originally developed for DNA-level annotation; applying them directly to RNA data can lead to erroneous assignments or inaccurate annotations (2–4).

Limited annotation flexibility: Tools such as exomePeak focus primarily on peak calling and basic gene annotation, offering little flexibility for modifying or supplementing annotation details (5).

Inconsistent results: Different tools produce markedly different annotation outputs, posing challenges for data integration and result validation.

Moreover, based on our previous work with database workflows, we observed that RNA‑specific complexities—such as overlapping genes arising from alternative splicing, multiple transcript isoforms, and strand‑specific information—are often overlooked by existing tools, further compromising annotation accuracy.

Key Features
RNA‑level specialization: Designed specifically for RNA modification sites, correctly handling overlapping genes and multi‑transcript annotation.

Strand‑aware annotation: Incorporates strand‑specific information during annotation, avoiding the misassignments common with DNA‑based tools applied to RNA data.

Accurate intergenic region handling: Minimizes bias by properly treating intergenic regions, thereby improving the reliability of downstream functional analyses (e.g., GO and KEGG enrichment).

Extensible molecular interaction annotation: Supports integration of potential molecular interactions (e.g., protein–RNA interactions), enabling users to explore possible regulatory mechanisms associated with modification sites.

Cross‑tool compatibility: Easily integrates outputs from various upstream analysis tools, standardizing annotation results for easier comparison and validation.

Installation
r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yourusername/RMAnno")
Quick Usage Example
r
library(RMAnno)

# Suppose you have a BED file of modification sites
sites <- read.table("mod_sites.bed")

# Perform annotation
result <- RMAnno(sites, genome = "hg38")

# View results
head(result)
For more detailed instructions, please refer to the package vignette:

r
vignette("RMAnno")
Citation
If you use RMAnno in your research, please cite our article (citation details will be added here).

License
GPL-3.0
