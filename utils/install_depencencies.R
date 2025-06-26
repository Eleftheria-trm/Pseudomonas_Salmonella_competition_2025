# Check and install required R packages for Pseudomonas-Salmonella competition analysis

# Required packages list
required_packages <- c(
  'readxl',      # Reading Excel files
  'dplyr',       # Data manipulation
  'tidyr',       # Data tidying
  'ggplot2',     # Plotting and visualization
  'DESeq2',      # Differential expression analysis
  'data.table',  # High-performance data operations
  'emmeans',     # Statistical modeling
  'lme4',        # Mixed-effects models
  'ggtext',      # Enhanced text in plots
  'patchwork',   # Combining plots
  'openxlsx'     # Excel output generation
)

# Check which packages are installed
installed_packages <- installed.packages()[,'Package']
missing_packages <- required_packages[!(required_packages %in% installed_packages)]

# Report status
cat("=== R Package Status Check ===\n\n")

if(length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse=', '), '\n\n')
} else {
  cat("✓ All required packages are installed!\n\n")
}

cat("Package status:\n")
for(pkg in required_packages) {
  if(pkg %in% installed_packages) {
    cat('✓', pkg, '- INSTALLED\n')
  } else {
    cat('✗', pkg, '- MISSING\n')
  }
}

# Install missing packages
if(length(missing_packages) > 0) {
  cat("\n=== Installing Missing Packages ===\n")
  
  # Set CRAN mirror to avoid mirror selection issues
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
  
  # Check if BiocManager is needed for DESeq2
  if('DESeq2' %in% missing_packages && !require(BiocManager, quietly = TRUE)) {
    cat("Installing BiocManager for Bioconductor packages...\n")
    install.packages("BiocManager")
  }
  
  for(pkg in missing_packages) {
    cat("Installing", pkg, "...\n")
    
    tryCatch({
      if(pkg == 'DESeq2') {
        BiocManager::install(pkg, quiet = TRUE)
      } else {
        install.packages(pkg, quiet = TRUE)
      }
      
      # Verify installation
      if(require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓", pkg, "installed successfully\n")
      } else {
        cat("✗ Failed to install", pkg, "\n")
      }
    }, error = function(e) {
      cat("✗ Error installing", pkg, ":", e$message, "\n")
    })
  }
  
  cat("\n=== Installation Complete ===\n")
} else {
  cat("\nNo packages need to be installed.\n")
}

cat("\nYou can now run the analysis scripts:\n")
cat("- Rscript 'Figure 2/Percentage script_with individual points.R'\n")
cat("- Rscript 'Figure 2/cfu.seedling_with individual points_script.R'\n") 
cat("- Rscript deseq/ReportForPaperFinal.R\n")