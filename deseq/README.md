# DESeq2 Differential Expression Analysis

This directory contains DESeq2 analysis of RNA-seq data comparing bacterial gene expression with and without Salmonella treatment.

## Files

- `ReportForPaperFinal.R` - Main DESeq2 analysis script for RNA-seq differential expression
- `In_vitro_all_reads.csv` - *In vitro* RNA-seq count data 
- `Plant raw reads.csv` - *In planta* RNA-seq count data
- `ReportForPaperFinal.html` - Generated HTML report with results

## Running the Analysis

```bash
# From project root
Rscript deseq/ReportForPaperFinal.R

# Or generate HTML report
R -e "rmarkdown::render('deseq/ReportForPaperFinal.R')"
```

## Analysis Overview

- Compares Control vs Treatment conditions in both environments
- Uses paired experimental design with interaction effects
- Identifies environment-specific gene expression responses
- Generates Excel reports with differential expression results