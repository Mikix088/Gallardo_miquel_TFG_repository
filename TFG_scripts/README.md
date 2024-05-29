# Brief scripts introduction

## GWAS Data Quality Control Script (00)

This script is designed to perform a series of quality control (QC) steps on genotype data from Genome-Wide Association Studies (GWAS). The script ensures the integrity and accuracy of the dataset by addressing common issues such as duplicate variants, incorrect SNP identifiers, and filtering based on minor allele frequency (MAF), Hardy-Weinberg equilibrium (HWE), and genotyping call rates.

The QC process involves:
- Removing duplicate SNPs based on their chromosome and base pair positions.
- Updating SNP identifiers using a reference file to ensure consistency.
- Filtering SNPs and samples based on predefined quality thresholds to ensure robust and reliable downstream analyses.

By following these steps, the script helps in maintaining high-quality genetic data, crucial for accurate GWAS results and subsequent biological interpretations.

## Multi-Data Integration and Preprocessing Script (01)

This script serves as a comprehensive solution for integrating and preprocessing genetic and phenotypic data for further analysis. It leverages R libraries such as readxl and data.table to efficiently handle data manipulation tasks.

Key functionalities include:

- Integration of multiple datasets into a unified format, ensuring consistency across different data sources.
- Quality assessment and cleaning procedures to identify and handle missing values, outliers, and inconsistencies.
- Generation of covariate and phenotype files required for downstream genetic association analyses.
- Exporting the cleaned and processed data in various formats for easy access and compatibility with statistical analysis tools.

By automating these preprocessing tasks, the script streamlines the data preparation pipeline, enabling researchers to focus on interpreting results and deriving insights from complex genetic and phenotypic datasets.

## GWAS Data Preparation Script (02)

This script is designed to prepare genotype data for Genome-Wide Association Studies (GWAS) after it's imputation in TopMed server, by performing various preprocessing steps. The script utilizes a combination of R and Bash commands to handle data loading, filtering, and formatting tasks efficiently.

Key functionalities of the script include:

1. **Data Loading and Quality Control**: The script loads genotype information from imputed data files and conducts quality control checks to exclude SNPs with low imputation quality (Rsq <= 0.6). It also generates lists of SNPs passing QC thresholds.

2. **Binary File Creation**: Using PLINK, the script creates binary files (.bed, .bim, .fam) for each chromosome after excluding SNPs with low imputation quality.

3. **Path List Generation**: It generates a list of file paths for the merged binary files to facilitate further processing steps.

4. **BIM File Transformation**: The script transforms the V2 column in BIM files to match the required format for downstream analysis.

5. **Merging Chromosomes**: It merges binary files for all chromosomes into a single dataset using PLINK.

6. **Data Filtering**: Finally, the script applies additional filtering criteria to remove samples with high genotyping missingness (geno >= 0.05) and low minor allele frequency (maf <= 0.01).

By automating these steps, the script streamlines the data preparation process, ensuring that the GWAS dataset is clean, formatted, and ready for subsequent analyses.

## GWAS & PCA Analysis Script (03)

This script processes genotype association data from GWAS studies, conducting Principal Component Analysis (PCA) to uncover population structure. It utilizes R packages such as `data.table`, `dplyr`, `SNPRelate`, and others for data manipulation and analysis.

1. **Data Loading**: Genotype association data is loaded and preprocessed to remove file extensions and ensure data integrity.

2. **Subset Generation**: Subsets containing specific information (e.g., collaterals' IDs) are created and saved for further analysis.

3. **PLINK Operations**: PLINK is used to extract genotype data based on IDs and perform LD-based SNP pruning.

4. **PCA Calculation**: PCA is conducted using the `snpgdsPCA` function, converting binary files to GDS format and calculating principal components.

5. **Population Structure Analysis**: PCA results are visualized to identify population clusters, allowing for the investigation of genetic associations with different phenotypic outcomes.

## GWAS Analysis Pipeline (04)

This repository contains a comprehensive pipeline for Genome-Wide Association Study (GWAS) analysis focused on stroke-related traits. It includes data preparation, linear model analysis, Genome-wide Complex Trait Analysis (GCTA), and quality control.
Using PLINK and R, the pipeline handles data preprocessing, applies linear models accounting for covariates, and conducts association tests. Quality control involves Hardy-Weinberg Equilibrium (HWE) tests. 
Results interpretation includes visualization and annotation of significant associations. This pipeline facilitates the systematic investigation of stroke genetics, aiding in understanding disease mechanisms and identifying potential therapeutic targets.

## Statistical Summary (05)

This R script analyzes data from the 4 variables of interest for the dichotomic ones include different categories (all cases, good cases, and bad cases). It loads data from an Excel file, subsets it based on specific participant IDs, and performs various calculations and analyses. The script includes:

- Loading the data from an Excel file using `read_excel`.
- Subsetting the data based on participant IDs.
- Processing and analyzing data for all cases, good cases, and bad cases separately.
- Calculating statistics such as percentage of females, median and quantiles of age, smoking status, aspects_b, hypertension, diabetes mellitus, alcohol consumption, NIHSS basal scores, mRS alta, mRS previo, NIHSS alta, cardio ischemic disease, previous atrial fibrillation, cerebral ischemic preconditioning, symptomatic hemorrhagic transformation, and TOAST classification.
- Generating p-values for various factors using chi-square tests and Kruskal-Wallis tests.

This script provides insights into different clinical factors and their association with the variables of interest outcomes, offering valuable information for further analysis and decision-making.

## TWAS (06)

### Loading Libraries
- The script starts by loading necessary libraries including `data.table`, `colochelpR`, `reshape2`, `dplyr`, and `SNPlocs.Hsapiens.dbSNP155.GRCh38`.

### Data Processing and Analysis
- The script iterates over different phenotypes (`neuros`, `collaterals`, `rank_transformation`, `brainomix`) and performs the following steps:
  - Loads summary statistics data from GWAS files.
  - Processes the data, extracting relevant columns like `rsid`, `CHR`, `BP`, `A1`, `A2`, `BETA`, and `SE`.
  - Filters out SNPs without rs IDs and performs allele flipping.
  - Matches summary stats with reference data and performs quality control.
  - Loads weights and iterates over each gene to compute TWAS Z-scores and p-values.
  - Optionally performs COLOC analysis to assess colocalization between TWAS and GWAS signals.
  - Outputs results to files and computes permutation tests if specified.

### Reading Results
- After TWAS and COLOC analysis, the script reads the results for each phenotype separately:
- For `brainomix`, `neuros`, `collaterals`, and `rank_transformation`, it reads TWAS results from corresponding files, computes FDR-adjusted p-values, and stores the data in separate data frames.

This script performs genome-wide association studies (GWAS), transcriptome-wide association studies (TWAS), and colocalization analyses to identify potential genetic associations with various phenotypes, providing valuable insights for further investigation.

## PWAS (07)

### Loading Packages
Loads necessary R packages including data.table, R.utils, reshape2, tidyr, plyr, bigsnpr, and coloc.

### Creating Recurrent Objects
 Defines recurrent objects such as the trait name, file paths for data and analysis, and cohort information.

### PWAS Analysis
 Creates folders for each cohort.
 Conducts PWAS analysis for each chromosome, generating PWAS result files.

### Loading Results
 Reads PWAS results for each cohort and phenotype.
 Combines PWAS results into a single data frame.
 Filters significant features based on adjusted p-value threshold.
 Identifies unique genes with the best p-value in each group.

### Additional Data Processing
 Matches Ensembl IDs with gene symbols using annotable package.
 Saves the PWAS results and significant features to text files.

 Provides insights into polygenic associations with the specified trait across different cohorts.

## Single-cell (08)

This script conducts a detailed analysis of RNA expression in endothelial and vascular smooth muscle cells to uncover molecular mechanisms underlying vascular dysfunction. 
It performs differential expression analysis on integrated single-cell RNA sequencing data to identify significant genes associated with exon mutations, sex differences, and cognitive impairment.
The analysis involves data integration, cell type identification, and statistical evaluation of gene expression differences. Significant findings are validated across discovery and replication datasets, 
with results formatted for further functional enrichment and interaction network analysis using the STRING database, offering insights into the biological pathways and processes involved in vascular health and disease.

## Locus Zoom Analysis Script (E1)

This script utilizes the `locuszoomr` package in R to perform locus zoom analysis on GWAS data. It also incorporates LDlinkR for linkage disequilibrium (LD) analysis. 

1. **Package Loading**: Necessary R packages like `data.table`, `locuszoomr`, `EnsDb.Hsapiens.v75`, and `LDlinkR` are loaded.

2. **Data Loading and Preprocessing**: GWAS data is loaded and filtered to include only alleles "A", "T", "C", and "G". Loci with p-values below 5e-8 are extracted.

3. **Locus Zoom Analysis**:
   - **UKBB Analysis**: Locus zoom analysis is performed on the GWAS data using the `locus` function, focusing on specific index SNPs.
   - **Collaterals Analysis**: Another locus zoom analysis is conducted for collaterals GWAS data, again focusing on specific index SNPs.

4. **LD Analysis**: LD analysis is performed using LDlinkR to explore the linkage disequilibrium patterns around the index SNPs.

5. **Visualization**: Locus zoom plots are generated to visualize the association signals and LD patterns around the index SNPs for both UKBB and collaterals datasets.

## GWAS Data Processing Script (E2)

This script preprocesses genotype association data from multiple GWAS studies, cleaning and standardizing the format for downstream analyses. It utilizes R packages like `data.table`, `dplyr`, and `SNPlocs.Hsapiens.dbSNP155.GRCh37` for efficient data manipulation and conversion.

Key steps performed by the script include:

1. **Data Loading**: Genotype association data from different GWAS studies are loaded into separate data frames (`gwas1`, `gwas2`, `gwas3`, `gwas4`) using the `fread` function from the `data.table` package.

2. **SNP Identifier Processing**: The script extracts numeric SNP identifiers (`no_rs`) from the SNP column and removes any rows with missing identifiers.

3. **Genomic Location Conversion**: SNP identifiers are converted to genomic locations (chromosome and position) in the GRCh37 assembly using the `convert_rs_to_loc` function from the `SNPlocs.Hsapiens.dbSNP155.GRCh37` package. Chromosome and position information is extracted and added to the data frames (`df1`, `df2`, `df3`, `df4`).

4. **Position Data Formatting**: Genomic positions are converted from float to integer format for consistency and compatibility with downstream analyses. A helper function `sci` identifies scientific notation positions and ensures proper formatting.

5. **Data Writing**: Processed data frames are written to compressed files (`geno_assoc_brainomix_clean_grch37.gwas.gz`, `geno_assoc_neuros_clean_grch37.gwas.gz`, `geno_assoc_collaterals_clean_grch37.gwas.gz`, `geno_assoc_rank_transformation_clean_grch37.gwas.gz`) using the `fwrite` function from the `data.table` package.

By standardizing genotype association data and converting SNP identifiers to genomic locations, this script prepares the data for further analysis, facilitating genetic variant discovery and association studies.

