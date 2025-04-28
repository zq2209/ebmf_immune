library(Matrix)
# library(MatrixExtra)
library(flashier)
library(fastTopics)
library(dplyr)
library(readr)


# Taken from https://github.com/stephenslab/pathways/blob/master/inst/code/read_gene_set_data.R
read_gene_info <- function (file) {

  # Read the data into a data frame.
  out <- suppressMessages(read_delim(file,delim = "\t",col_names = TRUE))
  class(out) <- "data.frame"
  dbXrefs    <- out$dbXrefs
  out        <- out[c("GeneID","Symbol","Synonyms","chromosome")]

  # Set any entries with a single hyphen to NA, and convert the
  # "chromosome" column to a factor.
  out$chromosome[out$chromosome == "-"] <- NA
  out$Synonyms[out$Synonyms == "-"]     <- NA
  dbXrefs[dbXrefs == "-"]               <- NA
  out <- transform(out,chromosome = factor(chromosome))

  # Extract the Ensembl ids. Note that a small number of genes map to
  # more than one Ensembl id; in those cases, we retain the first
  # Ensembl id only.
  dbXrefs <- strsplit(dbXrefs,"|",fixed = TRUE)
  out$Ensembl <- sapply(dbXrefs,function (x) {
    i <- which(substr(x,1,8) == "Ensembl:")
    if (length(i) > 0)
      return(substr(x[i[1]],9,nchar(x[i[1]])))
    else
      return(as.character(NA))
  })

  # For human genes, extract the HGNC (HUGO Gene Nomenclature
  # Committee) ids.
  out$HGNC <- sapply(dbXrefs,function (x) {
    i <- which(substr(x,1,10) == "HGNC:HGNC:")
    if (length(i) > 0)
      return(substr(x[i[1]],6,nchar(x[i[1]])))
    else
      return(as.character(NA))
  })

  # Return the processed gene data.
  return(out)
}

homo_sapien_geno_info <- read_gene_info('../data/Homo_sapiens.gene_info.gz')
data <- readRDS('../data/human_brain_COP_cells.rds')
counts <- t(data$RNA$data)

# Keep genes that are in the gene info file and remove those without nonzero counts

reduced_counts <-
  counts[, colnames(counts) %in% homo_sapien_geno_info$Ensembl]
cols_to_keep <- colSums(reduced_counts != 0, na.rm = TRUE) > 0
reduced_counts <- reduced_counts[, cols_to_keep]


# compute shifted log counts
n  <- nrow(reduced_counts)
x  <- rpois(1e7, 1/n)
s1 <- sd(log(x + 1))
a <- 1
size_factors <- rowSums(reduced_counts)
size_factors <- size_factors / mean(size_factors)
shifted_log_counts <- log1p(reduced_counts / (a * size_factors))
# shifted_log_counts <- mapSparse(reduced_counts / (a * size_factors),
#                                 fn = log1p)


flashier_fit <- flash(shifted_log_counts,
                      ebnm_fn = ebnm_point_exponential,
                      var_type = 2,
                      greedy_Kmax = 50,
                      S = s1,
                      backfit = TRUE)

flashier_fit_semi <- flash(shifted_log_counts,
                           ebnm_fn = c(ebnm_point_exponential, ebnm_point_laplace),
                           var_type = 2,
                           greedy_Kmax = 50,
                           S = s1,
                           backfit = TRUE)


fasttopics_fit_50 <- fit_topic_model(reduced_counts, k = 50)
fasttopics_fit_10 <- fit_topic_model(reduced_counts, k = 10)

save(flashier_fit, flashier_fit_semi, fasttopics_fit_50, fasttopics_fit_10,
     file = '../data/human_brain_COP_cells_fit.RData')
