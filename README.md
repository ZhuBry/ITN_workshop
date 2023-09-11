# README
# Acute Myeloid Leukemia Heatmap

This README provides an overview and instructions for running the Acute Myeloid Leukemia Heatmap analysis. The analysis is adapted from the [refine.bio-examples notebook](https://alexslemonade.github.io/refinebio-examples/03-rnaseq/clustering_rnaseq_01_heatmap.html) and uses the [acute myeloid leukemia sample dataset](https://www.refine.bio/experiments/SRP070849) from Shih et al., 2017. The purpose of the analysis is to cluster and visualize gene expression data for different types of AML under controlled treatment conditions.

## Setup

Before running the analysis, set up the necessary folders for storing data, plots, and results. The following code snippet creates the required folders if they don't already exist:

```R
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Create the plots folder if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Create the results folder if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}
```

## Install Libraries

To perform clustering and create a heatmap, we will use the `pheatmap` package by Slowikowski et al., 2017. If the package is not already installed, you can install it using the following code:

```R
if (!("pheatmap" %in% installed.packages())) {
  install.packages("pheatmap", update = FALSE)
}
```

After installing the package, load the `pheatmap` and `magrittr` libraries:

```R
library(pheatmap)
library(magrittr)
```

## Import and Set Up Data

The next step is to import and set up the data for analysis. The code snippet below reads in the metadata and gene expression data from TSV files and prepares them for clustering:

```R
# Read in metadata TSV file
metadata <- readr::read_tsv("data/metadata_SRP070849.tsv")

# Read in data TSV file and set gene IDs as row names
expression_df <- readr::read_tsv("data/SRP070849/SRP070849.tsv") %>%
  tibble::column_to_rownames("Gene")

# Ensure the metadata and data are in the same sample order
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if the sample order is the same
all.equal(colnames(expression_df), metadata$refinebio_accession_code)
```

## Clustering Heatmap

The code below performs clustering and creates a heatmap using the `pheatmap` function:

```R
# Choose genes of interest based on variance
variances <- apply(expression_df, 1, var)
upper_var <- quantile(variances, 0.75)
df_by_var <- data.frame(expression_df) %>%
  dplyr::filter(variances > upper_var)

# Save the filtered genes to a TSV file
readr::write_tsv(df_by_var, "results/top_90_var_genes.tsv")

# Prepare metadata for annotation
annotation_df <- metadata %>%
  dplyr::mutate(
    mutation = dplyr::case_when(
      startsWith(refinebio_title, "TET2") ~ "TET2",
      startsWith(refinebio_title, "IDH2") ~ "IDH2",
      startsWith(refinebio_title, "WT") ~ "WT",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(
    refinebio_accession_code,
    mutation,
    refinebio_treatment
  ) %>%
  tibble::column_to_rownames("refinebio_accession_code")

# Create annotated heatmap
heatmap_annotated <- pheatmap(df_by_var, cluster_rows = TRUE, cluster_cols = TRUE,
                              annotation_col = annotation_df)
```

## Additional Resources

For more information and examples on creating heatmaps using the `pheatmap` package, you can refer to the following resources:

- [Making a heatmap in R with the pheatmap package](https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/) [^0^]
- [Heatmap with pheatmap package in R](https://www.reneshbedre.com/blog/heatmap-with-pheatmap-package-r.html) [^2^]
- [pheatmap function from the pheatmap package](https://biocorecrg.github.io/CRG_RIntroduction_2021/pheatmap-function-from-the-pheatmap-package.html) [^6^]
- [How to Use pheatmap() in R to Create Heatmaps](https://www.statology.org/pheatmap-r