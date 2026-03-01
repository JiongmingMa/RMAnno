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
devtools::install_github("yourusername/RMAnno")
```

### From source (if you have the tar.gz)
```r
install.packages("/path/to/RMAnno.tar.gz", repos = NULL, type = "source")
```

---

## 🚀 Quick Start

Here's a minimal example to annotate RNA modification sites from a BED file:

```r
library(RMAnno)

# Read your modification sites (BED format)
sites <- read.table("mod_sites.bed", header = FALSE)

# Annotate using human genome hg38
result <- RMAnno(sites, genome = "hg38")

# Check the first few rows
head(result)
```

**Expected output** (example):

| chr | start   | end     | strand | gene_id | transcript_id | region   |
|-----|---------|---------|--------|---------|---------------|----------|
| chr1| 12345   | 12346   | +      | ENSG001 | ENST001       | CDS      |
| chr1| 67890   | 67891   | -      | ENSG002 | ENST002       | 3'UTR    |

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

or visit the [online documentation](https://yourusername.github.io/RMAnno/).

---

## 📝 Citation

If you use RMAnno in your research, please cite:

> Your Name, et al. (2026). RMAnno: a dedicated tool for RNA modification annotation. *Bioinformatics*, XX(X), XXX-XXX. DOI: 10.1093/bioinformatics/xxxxx

---

## 📄 License

This package is licensed under the **GPL-3.0** license. See the [LICENSE](LICENSE) file for details.

---

## 🙌 Contributing

Contributions are welcome! If you encounter a bug or have a feature request, please open an [issue](https://github.com/yourusername/RMAnno/issues). Pull requests are also appreciated.

---

*Made with ❤️ by the RMAnno team*
