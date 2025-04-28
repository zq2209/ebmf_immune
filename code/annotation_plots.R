# This function creates an "effect plot". The effects_matrix input
# should be a matrix in which rows are features and columns are
# dimensions (e.g., factors in a matrix factorization). Effects
# smaller than zero_value in magnitude are not included in the plot.
effect_plot <- function (effects_matrix, font_size = 9,
                         zero_value = 0.01) {
  features <- rownames(effects_matrix)
  pdat <- data.frame(feature_name = features,
                     stringsAsFactors = FALSE)
  pdat <- cbind(pdat,effects_matrix)
  pdat <- melt(pdat,id.vars = "feature_name",variable.name = "dim",
               value.name = "value")
  pdat <- transform(pdat,
                    effect_size = abs(value),
                    effect_sign = factor(sign(value),c(-1,0,1)),
                    feature_name = factor(feature_name,rev(features)))
  pdat$effect_size[pdat$effect_size < zero_value] <- NA
  return(ggplot(pdat,aes(x = dim,y = feature_name,size = effect_size,
                         fill = effect_sign)) +
         geom_point(color = "white",shape = 21,na.rm = TRUE) +
         scale_size(range = c(1,6),
                    breaks = range(pdat$effect_size,na.rm = TRUE),
                    labels = round(range(pdat$effect_size,na.rm = TRUE),
                                   digits = 2)) +
         scale_fill_manual(values = c("navy","gray","orangered"),
                           drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(size = 3))) +
         labs(x = "dimension",y = "feature name",
              fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size))
}

# This function selects the top "driving" genes for the selected
# dimensions ("dims"). The effects_matrix input should be a matrix in
# which rows are genes and columns are dimensions. It is assumed the
# names of the rows give the names or ids of the genes. If
# select_down_effects = TRUE, the genes with the largest negative
# effects are also included. Note that a gene is never selected more
# than once.
select_driving_genes <- function (effect_matrix, dims, n,
                                  select_down_effects = TRUE) {
  genes <- NULL
  for (i in dims) {
    genes <- c(genes,head(order(effect_matrix[,i],decreasing = TRUE),n))
    if (select_down_effects)
      genes <- c(genes,head(order(effect_matrix[,i],decreasing = FALSE),n))
  }
  genes <- unique(genes)
  return(rownames(effect_matrix)[genes])
}

# This creates an effect plot for a matrix factorization model, and
# adds a layer of automation by handpicking the genes (features) to
# show in the effect plot.  This may be overriden with the "genes"
# argument.
driving_genes_heatmap <-
  function (effects_matrix, dims = 1:ncol(effects_matrix), n = 3,
            genes = select_driving_genes(effects_matrix,dims = dims,n = n))
  effect_plot(F[genes,]) +
    labs(x = "factor",y = "gene",size = "size",fill = "sign")

# Compute the "least extreme" (l.e.) effect differences. The
# l.e. effect difference for dimension k is defined as the smallest
# difference between the effect of dimension k and the effect of any
# other dimension compared (as specified by the "compare_dims"
# argument). Note that this only works for a non-negative effects
# matrix, e.g., a non-negative matrix factorization.
compute_le_diff <- function (effects_matrix,
                             compare_dims = 1:ncol(effects_matrix)) {
  n <- nrow(effects_matrix)
  k <- ncol(effects_matrix)
  le_effects <- matrix(0,n,k)
  rownames(le_effects) <- rownames(effects_matrix)
  colnames(le_effects) <- colnames(effects_matrix)
  for (i in 1:k) {
    le_effects[,i] <-
      effects_matrix[,i] -
        apply(effects_matrix[,setdiff(compare_dims,i)],1,max)
  }
  return(le_effects)
}

# This creates a "distinctive genes plot"; this is a plot in which the
# effect estimate is shown on the x axis and the "least extreme"
# difference between the estimated effects is shown on the y axis. The
# idea is that these scatterplots should better highlight the
# "interesting" genes for a given dimension/factor. The "label_gene"
# input argument is a function that returns TRUE when the gene should
# be labeled in the plot; the default is that it always returns FALSE
# (so that no genes are labeled in the plot).
distinctive_genes_scatterplot <-
  function (effects_matrix, k, compare_dims = 1:ncol(effects_matrix),
            font_size = 12, label_size = 2.25, max_overlaps = Inf,
            label_gene = function (x, y) rep(FALSE,length(x))) {
  le_diff <- compute_le_diff(effects_matrix,compare_dims)
  genes   <- rownames(effects_matrix)
  pdat    <- data.frame(gene    = genes,
                        effect  = effects_matrix[,k],
                        le_diff = le_diff[,k])
  i <- which(!label_gene(pdat$effect,pdat$le_diff))
  pdat[i,"gene"] <- NA
  return(ggplot(pdat,aes(x = effect,y = le_diff,label = gene)) +
         geom_point(color = "dodgerblue") +
         geom_hline(yintercept = 0,color = "magenta",linetype = "dotted",
                    linewidth = 0.5) +
         geom_text_repel(color = "black",size = label_size,
                         fontface = "italic",segment.color = "black",
                         segment.size = 0.25,min.segment.length = 0,
                         max.overlaps = max_overlaps,na.rm = TRUE) +
         labs(x = "effect estimate",y = "least extreme difference") +
         theme_cowplot(font_size = font_size))
}
