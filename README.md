# RMAnno 🧬

![License](https://img.shields.io/github/license/yourusername/RMAnno)
![R](https://img.shields.io/badge/made%20with-R-blue)
![GitHub release](https://img.shields.io/github/v/release/yourusername/RMAnno)

**RMAnno** is an R package specifically designed for the **annotation of RNA-level modification sites**. It addresses the inaccuracies, information loss, and integration difficulties that arise when existing tools—mostly derived from DNA-level annotation methods—are applied to RNA modification data. RMAnno provides researchers with an accurate, flexible, and reproducible workflow for annotating RNA modifications.

---

## 📖 Table of Contents

- [Background](#background)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Citation](#citation)
- [License](#license)

---

## 🔬 Background

RNA modification analysis typically involves two phases: upstream (quality control, alignment) and downstream (functional annotation, pathway analysis). **Genomic annotation** is a fundamental step in downstream analysis, directly affecting the reliability of subsequent enrichment analyses such as GO and KEGG. However, widely used annotation tools suffer from several limitations:

| Limitation | Description |
|------------|-------------|
| **Misaligned tool design** | Tools like ChIPseeker and RNAmod were developed for DNA-level annotation; applying them directly to RNA data can lead to erroneous assignments. |
| **Limited annotation flexibility** | Tools such as exomePeak focus primarily on peak calling, offering little flexibility for modifying annotation details. |
| **Inconsistent results** | Different tools produce markedly different annotation outputs, posing challenges for data integration and result validation. |

Moreover, RNA-specific complexities—such as overlapping genes arising from alternative splicing, multiple transcript isoforms, and strand-specific information—are often overlooked by existing tools, further compromising annotation accuracy.

---

## ✨ Key Features

- ✅ **RNA-level specialization** – correctly handles **overlapping genes** and **multi-transcript** annotation.
- ✅ **Strand-aware annotation** – incorporates strand-specific information, avoiding misassignments common with DNA-based tools.
- ✅ **Accurate intergenic region handling** – minimizes bias and improves downstream functional analyses (GO/KEGG).
- ✅ **Extensible molecular interaction annotation** – supports integration of potential protein–RNA interactions.
- ✅ **Cross-tool compatibility** – easily integrates outputs from various upstream analysis tools.

---

## 📦 Installation

### From GitHub (recommended)
```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("JiongmingMa/RMAnno")
```

### From source (if you have the tar.gz)
```r
install.packages("/path/to/RMAnno.tar.gz", repos = NULL, type = "source")
```

---

## 🚀 Quick Start

Here's a minimal example to annotate RNA modification sites from a BED file:

```r
library(GenomicRanges)
library(RMAnno)

# Read your modification sites (here just example)
gr1 <- GRanges(seqnames = "1", ranges = IRanges(start = c(51579900, 51579900,51180000, 51180000), width = 1), strand = c("+","-","+","-"))

# Extract the information from GTF file (I use human GRCh38.114 GTF file here, from: https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.chr.gtf.gz )
gtf_file <- ExtractGTFFile("path/to/your/GTF_file/Homo_sapiens.GRCh38.114.chr.gtf.gz")

# Annotate using human genome hg38
gr1 <- RMAnno::rmAnno(gr1, gtf_file)

# Check the first few rows
gr1
```

**Expected output** (example):

```plaintext
GRanges object with 4 ranges and 5 metadata columns:
      seqnames    ranges strand |      Region                Gene_id        Gene_name               Bio_type          transcript_id
         <Rle> <IRanges>  <Rle> | <character>            <character>      <character>            <character>            <character>
  [1]        1  51579900      + |      Intron ENSG00000227070/ENSG.. EPS15-AS1/OSBPL9  lncRNA/protein_coding ENST00000371714/ENST..
  [2]        1  51579900      - |      Intron ENSG00000293506/ENSG..        NA/CALR4P lncRNA/transcribed_u.. ENST00000431943/ENST..
  [3]        1  51180000      + |  Intergenic                   <NA>             <NA>                   <NA>                   <NA>
  [4]        1  51180000      - |  Intergenic                   <NA>             <NA>                   <NA>                   <NA>
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

---

## 📚 Detailed Usage

For a comprehensive guide, including:

- How to customize annotation databases
- Working with GRanges objects
- Batch annotation of multiple samples
- Integration with downstream enrichment tools

Please see the package vignette:

```r
browseVignettes("RMAnno")
```

or visit the [online documentation](https://JiongmingMa.github.io/RMAnno/).

---


## 📄 License

This package is licensed under the **GPL-3.0** license. See the [LICENSE](LICENSE) file for details.

---

