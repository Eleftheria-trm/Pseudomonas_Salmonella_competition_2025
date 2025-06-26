# Pseudomonas-Salmonella Competition Analysis

This repository contains R scripts and data for bacterial competition experiments between Pseudomonas and Salmonella, supporting the research findings published in the companion paper (see citation)

## Directory Contents

### :file_folder: [Figure 2](./Figure%202/README.md)

Contains scripts for visualizing bacterial colony counting data.

### :file_folder: [deseq](./deseq/README.md)

Contains RNA-seq differential expression analysis.

## Usage

Run R scripts individually:

```bash
Rscript "Figure 2/Percentage script_with individual points.R"
Rscript "Figure 2/cfu.seedling_with individual points_script.R"
Rscript deseq/ReportForPaperFinal.R
```

## Requirements

R packages: `readxl`, `dplyr`, `tidyr`, `ggplot2`, `DESeq2`, `data.table`, `emmeans`, `lme4`, `ggtext`, `patchwork`, `openxlsx`

Note that a `./utils/install_dependencies.R` script can help install these.

## Citation

Vimont, N., Bastkowski, S., Savva, G. M., Bloomfield, S. J., Mather, A. E., Webber, M. A., & Trampari, E. (2025).  
*Environmental context as a key driver of Pseudomonas' biocontrol activity against Salmonella*.  
**bioRxiv**. [https://doi.org/10.1101/2025.06.23.661019](https://doi.org/10.1101/2025.06.23.661019)  
[Full text PDF](https://www.biorxiv.org/content/early/2025/06/23/2025.06.23.661019.full.pdf)  
[View at bioRxiv](https://www.biorxiv.org/content/early/2025/06/23/2025.06.23.661019)  


