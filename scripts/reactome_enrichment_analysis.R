# Reactome Enrichment Analysis for Significantly Rhythmic Genes
# Author: Script for analyzing rhythmic genes from regression analysis
# Date: 2025-10-15

library(dplyr)
library(readr)

# Set working directory

# Create output directory if it doesn't exist
if (!dir.exists("data/enrichment")) {
  dir.create("data/enrichment", recursive = TRUE)
}

# Read Reactome annotation file
cat("Loading Reactome annotation...\n")
reactome_anno <- read_tsv("data/Ref/Reactome_anno_HS.txt",
                          col_names = c("ENSG", "Pathway_ID", "URL", "Pathway_Name", "Evidence", "Species"),
                          show_col_types = FALSE)

cat("Total Reactome annotations:", nrow(reactome_anno), "\n")
cat("Unique genes:", length(unique(reactome_anno$ENSG)), "\n")
cat("Unique pathways:", length(unique(reactome_anno$Pathway_ID)), "\n\n")

# Create a pathway-to-genes mapping
pathway_genes <- reactome_anno %>%
  select(Pathway_ID, Pathway_Name, ENSG) %>%
  distinct() %>%
  group_by(Pathway_ID, Pathway_Name) %>%
  summarise(genes = list(ENSG), .groups = "drop") %>%
  mutate(pathway_size = sapply(genes, length))

# Get background gene set (all genes in Reactome)
background_genes <- unique(reactome_anno$ENSG)
N_background <- length(background_genes)

cat("Background genes in Reactome:", N_background, "\n")
cat("Pathways with gene counts:", nrow(pathway_genes), "\n\n")

# List all sig_genes_regression_periodic_12h files
gene_files <- list.files("data/Merged/regression",
                         pattern = "^sig_genes_regression_periodic_12h.*\\.txt$",
                         full.names = TRUE)

cat("Found", length(gene_files), "gene files to process:\n")
print(basename(gene_files))

# Initialize a list to store all results
all_results <- list()

# Fisher's exact test function for enrichment
fisher_enrichment <- function(query_genes, pathway_genes_list, all_pathway_genes, N_background) {
  # Number of query genes in pathway
  m <- length(intersect(query_genes, pathway_genes_list))

  # Total number of query genes
  M <- length(query_genes)

  # Total number of genes in pathway
  n <- length(pathway_genes_list)

  # Total background genes
  N <- N_background

  # Only test if there's overlap
  if (m == 0) {
    return(NULL)
  }

  # Create contingency table
  # m: genes in both query and pathway
  # M-m: genes in query but not in pathway
  # n-m: genes in pathway but not in query
  # N-M-n+m: genes in neither

  contingency_matrix <- matrix(c(m, M - m, n - m, N - M - n + m), nrow = 2)

  # Fisher's exact test
  fisher_result <- fisher.test(contingency_matrix, alternative = "greater")

  # Calculate enrichment ratio
  expected_ratio <- n / N
  observed_ratio <- m / M
  enrichment_ratio <- observed_ratio / expected_ratio

  return(list(
    m = m,
    M = M,
    n = n,
    N = N,
    p_value = fisher_result$p.value,
    enrichment_ratio = enrichment_ratio
  ))
}

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

  # Extract ENSG IDs (remove version numbers)
  ensg_ids <- gsub("\\..*", "", gene_data$ENSG)

  # Filter to genes that are in Reactome background
  query_genes <- intersect(ensg_ids, background_genes)

  cat("Number of genes in Reactome background:", length(query_genes), "\n")

  if (length(query_genes) == 0) {
    cat("No genes found in Reactome background, skipping...\n")
    next
  }

  # Test each pathway
  pathway_results <- list()

  for (i in 1:nrow(pathway_genes)) {
    pathway_id <- pathway_genes$Pathway_ID[i]
    pathway_name <- pathway_genes$Pathway_Name[i]
    pathway_gene_list <- pathway_genes$genes[[i]]

    result <- fisher_enrichment(query_genes, pathway_gene_list, background_genes, N_background)

    if (!is.null(result)) {
      pathway_results[[pathway_id]] <- c(
        Pathway_ID = pathway_id,
        Pathway_Name = pathway_name,
        result
      )
    }
  }

  if (length(pathway_results) == 0) {
    cat("No pathways with overlapping genes found\n")
    next
  }

  # Convert to data frame
  results_df <- bind_rows(pathway_results)

  # Convert character columns to appropriate types
  results_df$m <- as.numeric(results_df$m)
  results_df$M <- as.numeric(results_df$M)
  results_df$n <- as.numeric(results_df$n)
  results_df$N <- as.numeric(results_df$N)
  results_df$p_value <- as.numeric(results_df$p_value)
  results_df$enrichment_ratio <- as.numeric(results_df$enrichment_ratio)

  # Calculate adjusted p-values (Benjamini-Hochberg)
  results_df$q_value <- p.adjust(results_df$p_value, method = "BH")

  # Filter for significant results (p < 0.05)
  significant_results <- results_df %>%
    filter(p_value < 0.05) %>%
    arrange(p_value)

  cat("Found", nrow(significant_results), "significant Reactome pathways (p < 0.05)\n")

  if (nrow(significant_results) > 0) {
    # Calculate -log10(p-value) for comparison with GO/KEGG results
    significant_results$neg_log10_pvalue <- -log10(significant_results$p_value)

    # Create output format similar to GO/KEGG results
    output_df <- data.frame(
      Group = group_name,
      ID = significant_results$Pathway_ID,
      GO_Term = significant_results$Pathway_Name,
      Term = "Reactome pathway",
      m = significant_results$m,
      M = significant_results$M,
      n = significant_results$n,
      N = significant_results$N,
      E_ratio = significant_results$enrichment_ratio,
      P_value = significant_results$neg_log10_pvalue,
      stringsAsFactors = FALSE
    )

    # Sort by P_value (descending, as it's -log10)
    output_df <- output_df[order(output_df$P_value, decreasing = TRUE), ]

    # Store results
    all_results[[group_name]] <- output_df

    # Print top 5 pathways
    cat("\nTop 5 enriched Reactome pathways:\n")
    print(head(significant_results[, c("Pathway_ID", "Pathway_Name", "m", "n", "p_value", "q_value")], 5))
  }
}

# Combine all results
if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)

  # Save combined results
  output_file <- "data/enrichment/full_term_Reactome_Regression_CHIRAL_12h.txt"
  write_tsv(combined_results, output_file)

  cat("\n========================================\n")
  cat("ANALYSIS COMPLETE!\n")
  cat("========================================\n")
  cat("Total groups processed:", length(all_results), "\n")
  cat("Total Reactome pathways found:", nrow(combined_results), "\n")
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
