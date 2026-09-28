# KEGG Enrichment Analysis for Significantly Rhythmic Genes
# Author: Script for analyzing rhythmic genes from regression analysis
# Date: 2025-10-15

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# Set working directory

# Create output directory if it doesn't exist
if (!dir.exists("data/enrichment")) {
  dir.create("data/enrichment", recursive = TRUE)
}

# List all sig_genes_regression_periodic_12h files
gene_files <- list.files("data/Merged/regression",
                         pattern = "^sig_genes_regression_periodic_12h.*\\.txt$",
                         full.names = TRUE)

cat("Found", length(gene_files), "gene files to process:\n")
print(basename(gene_files))

# Initialize a list to store all results
all_results <- list()

# Process each file
for (file_path in gene_files) {
  cat("\n========================================\n")
  cat("Processing:", basename(file_path), "\n")
  cat("========================================\n")

  # Extract group name from filename
  group_name <- basename(file_path)

  # Read the gene file
  gene_data <- read_tsv(file_path, show_col_types = FALSE)

  cat("Number of genes:", nrow(gene_data), "\n")

  # Extract gene symbols and ENSG IDs
  gene_symbols <- gene_data$Gene
  ensg_ids <- gene_data$ENSG

  # Convert gene symbols to ENTREZ IDs
  # First try with gene symbols
  gene_entrez <- bitr(gene_symbols,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

  # If some genes didn't convert, try with ENSG IDs
  missing_symbols <- setdiff(gene_symbols, gene_entrez$SYMBOL)
  if (length(missing_symbols) > 0) {
    cat("Trying to convert", length(missing_symbols), "missing genes using ENSG IDs...\n")
    missing_ensg <- ensg_ids[gene_symbols %in% missing_symbols]
    # Remove version numbers from ENSG IDs
    missing_ensg_clean <- gsub("\\..*", "", missing_ensg)

    gene_entrez_ensg <- bitr(missing_ensg_clean,
                             fromType = "ENSEMBL",
                             toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)

    # Combine the results
    if (nrow(gene_entrez_ensg) > 0) {
      # Just combine the ENTREZ IDs
      gene_entrez <- unique(c(gene_entrez$ENTREZID, gene_entrez_ensg$ENTREZID))
    } else {
      gene_entrez <- gene_entrez$ENTREZID
    }
  } else {
    gene_entrez <- gene_entrez$ENTREZID
  }

  cat("Successfully converted", length(unique(gene_entrez)), "genes to ENTREZ IDs\n")

  # Perform KEGG enrichment analysis
  if (length(unique(gene_entrez)) > 0) {
    kegg_result <- enrichKEGG(gene = unique(gene_entrez),
                              organism = 'hsa',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)

    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      cat("Found", nrow(kegg_result@result), "significant KEGG pathways\n")

      # Extract results and calculate enrichment ratio
      kegg_df <- as.data.frame(kegg_result)

      # Parse GeneRatio and BgRatio to calculate enrichment ratio
      kegg_df$GeneRatio_num <- sapply(strsplit(kegg_df$GeneRatio, "/"), function(x) as.numeric(x[1]))
      kegg_df$GeneRatio_denom <- sapply(strsplit(kegg_df$GeneRatio, "/"), function(x) as.numeric(x[2]))
      kegg_df$BgRatio_num <- sapply(strsplit(kegg_df$BgRatio, "/"), function(x) as.numeric(x[1]))
      kegg_df$BgRatio_denom <- sapply(strsplit(kegg_df$BgRatio, "/"), function(x) as.numeric(x[2]))

      # Calculate enrichment ratio (similar to E_ratio in GO results)
      kegg_df$E_ratio <- (kegg_df$GeneRatio_num / kegg_df$GeneRatio_denom) /
                         (kegg_df$BgRatio_num / kegg_df$BgRatio_denom)

      # Calculate -log10(p-value) for comparison with GO results
      kegg_df$neg_log10_pvalue <- -log10(kegg_df$pvalue)

      # Create output format similar to GO results
      output_df <- data.frame(
        Group = group_name,
        ID = kegg_df$ID,
        GO_Term = kegg_df$Description,  # Using Description as pathway name
        Term = "KEGG pathway",  # Indicating this is a KEGG pathway
        m = kegg_df$GeneRatio_num,  # Number of genes in overlap
        M = kegg_df$GeneRatio_denom,  # Total genes in query
        n = kegg_df$BgRatio_num,  # Total genes in pathway
        N = kegg_df$BgRatio_denom,  # Total genes in background
        E_ratio = kegg_df$E_ratio,
        P_value = kegg_df$neg_log10_pvalue,
        stringsAsFactors = FALSE
      )

      # Sort by P_value (descending, as it's -log10)
      output_df <- output_df[order(output_df$P_value, decreasing = TRUE), ]

      # Store results
      all_results[[group_name]] <- output_df

      # Print top 5 pathways
      cat("\nTop 5 enriched KEGG pathways:\n")
      print(head(kegg_df[, c("ID", "Description", "pvalue", "qvalue")], 5))
    } else {
      cat("No significant KEGG pathways found for this group\n")
    }
  } else {
    cat("No genes could be converted to ENTREZ IDs\n")
  }
}

# Combine all results
if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)

  # Save combined results
  output_file <- "data/enrichment/full_term_KEGG_Regression_CHIRAL_12h.txt"
  write_tsv(combined_results, output_file)

  cat("\n========================================\n")
  cat("ANALYSIS COMPLETE!\n")
  cat("========================================\n")
  cat("Total groups processed:", length(all_results), "\n")
  cat("Total KEGG pathways found:", nrow(combined_results), "\n")
  cat("Results saved to:", output_file, "\n")

  # Print summary statistics
  cat("\nSummary by group:\n")
  summary_table <- combined_results %>%
    group_by(Group) %>%
    summarise(
      n_pathways = n(),
      mean_pvalue = mean(P_value),
      max_pvalue = max(P_value)
    ) %>%
    arrange(desc(n_pathways))

  print(summary_table)

} else {
  cat("\nNo enrichment results found for any group\n")
}
