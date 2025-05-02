plot_factors <- function(x, row_names, col_names, row_order, col_order,
                         row_title = "Trait", col_title = "Factor"){
  if(is.null(x)) return(NULL)
  if(missing(row_names)) row_names <-  seq(nrow(x))
  if(missing(col_names)) col_names <-  seq(ncol(x))
  if(missing(row_order)) row_order <- seq(nrow(x))
  if(missing(col_order)) col_order <- seq(ncol(x))
  rownames(x) <- row_names
  colnames(x) <- col_names
  x <- x[row_order, col_order]
  meltx <- melt(x)
  names(meltx)[1:2] <- c("R", "C")
  meltx <- meltx %>%
    mutate( R = factor(R, levels = row_names[row_order]),
            C = factor(C, levels = col_names[col_order]))
  ggplot(data = meltx, aes(x=C, y=R, fill=value)) +
    ylab(row_title) +
    xlab(col_title) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle=90))
}


#' Prepare Annotation Data for Heatmap Visualization
#'
#' Creates annotation data structures needed for pheatmap with biological system grouping
#'
#' @param F A factor matrix (rows = traits, columns = factors)
#' @param affected_system A character vector specifying biological systems for each trait
#' @param color_palette A color palette function (default: viridis option "D")
#' @param factor_prefix Prefix for factor column names (default: "k")
#'
#' @return A list containing:
#'   - F_ordered: The factor matrix ordered by biological system
#'   - annotation_df: Data frame for pheatmap annotation
#'   - annotation_colors: Color mapping for pheatmap
#'   - group_colors: Named vector of group colors
#'
#' @examples
#' # Prepare annotation data
#' prep <- prepare_annotation_data(F_matrix, affected_systems)
#' # Use with pheatmap:
#' pheatmap(prep$F_ordered, annotation_row = prep$annotation_df,
#'          annotation_colors = prep$annotation_colors)
plot_factor_annotation_data <- function(F, affected_system,
                                    color_palette = function(n) viridis::viridis(n, option = "D"),
                                    factor_prefix = "k") {

  # Input validation
  if (nrow(F) != length(affected_system)) {
    stop("Length of affected_system must match number of rows in F")
  }

  if (!is.matrix(F) && !is.data.frame(F)) {
    stop("F must be a matrix or data.frame")
  }

  # Create annotation data frame
  annotation_df <- data.frame(Groups = affected_system)
  rownames(annotation_df) <- rownames(F)

  # Add factor names if they don't exist
  if (is.null(colnames(F))) {
    colnames(F) <- paste0(factor_prefix, 1:ncol(F))
  }

  # Order by biological system
  trait_order <- order(annotation_df$Groups)
  F_ordered <- F[trait_order, , drop = FALSE]
  annotation_df_ordered <- annotation_df[trait_order, , drop = FALSE]

  # Create color mapping
  group_levels <- unique(annotation_df_ordered$Groups)
  group_colors <- setNames(color_palette(length(group_levels)), group_levels)
  annotation_colors <- list(Groups = group_colors)

  # Return all components needed for plotting
  return(list(
    F_ordered = F_ordered,
    annotation_df = annotation_df_ordered,
    annotation_colors = annotation_colors,
    group_colors = group_colors
  ))
}
