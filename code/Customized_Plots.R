# DimPlotSagnik: Creates a customized UMAP or PCA plot for high-dimensional data visualization, 
# specifically for single-cell data. It allows for various plot customizations, such as point color, 
# size, cell grouping, and shape by metadata. Additional options include faceting by metadata columns, 
# shuffling plot order, and adding labels to clusters. The function also supports combining multiple 
# plots to allow grouped visualizations.

DimPlotSagnik <- function(
    data,
    dims = c(1, 2),         # Dimensions to plot (e.g., UMAP1 and UMAP2)
    cells = NULL,           # Subset of cells to plot
    cols = NULL,            # Colors for different groups
    pt.size = NULL,         # Point size for plotting
    reduction = NULL,       # Dimensional reduction technique (e.g., "umap" or "pca")
    group.by = NULL,        # Metadata column to color by
    split.by = NULL,        # Metadata column to facet by
    shape.by = NULL,        # Metadata column to shape by
    order = NULL,           # Custom order for plotting
    shuffle = FALSE,        # Option to shuffle data points to randomize plotting order
    seed = 1,               # Random seed for shuffling
    label = FALSE,          # Option to add labels to clusters
    label.size = 4,         # Size of cluster labels
    label.color = 'black',  # Color of the labels
    label.box = FALSE,      # Option to put labels in a box
    repel = FALSE,          # Use repelling labels to avoid overlap
    cells.highlight = NULL, # Subset of cells to highlight
    cols.highlight = '#DE2D26', # Color for highlighted cells
    sizes.highlight = 1,    # Size of highlighted cells
    na.value = 'grey50',    # Color for NA values
    ncol = NULL,            # Number of columns for faceting
    combine = TRUE,         # Combine plots into one
    raster = NULL,          # Option to use raster graphics
    raster.dpi = c(512, 512) # Resolution for raster graphics
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  colnames(data) <- paste0("UMAP", dims)
  data <- as.data.frame(x = data)
  dims <- paste0("UMAP", dims)
  data <- cbind(data, group.by)
  orig.groups <- group.by
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      plot <- if (is.null(x = orig.groups)) {
        plot + labs(title = NULL)
      } else {
        plot + labs(title = NULL)
      }
    }
  )
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

# DimPlotContinuous: A customized plotting function to visualize high-dimensional data, typically UMAP or PCA,
# with an option for continuous color scales. If the `group.by` parameter is a numeric vector, the function 
# uses a continuous color gradient, where higher values are represented by deeper colors. This function is 
# designed for flexible plotting with additional options for faceting, highlighting specific cells, labeling 
# clusters, and shuffling data points for randomized plotting order.

DimPlotContinuous <- function(
    data,
    dims = c(1, 2),         # Dimensions to plot (e.g., UMAP1 and UMAP2)
    cells = NULL,           # Subset of cells to plot
    cols = NULL,            # Colors for different groups
    pt.size = NULL,         # Point size for plotting
    reduction = NULL,       # Dimensional reduction technique (e.g., "umap" or "pca")
    group.by = NULL,        # Metadata column to color by
    split.by = NULL,        # Metadata column to facet by
    shape.by = NULL,        # Metadata column to shape by
    order = NULL,           # Custom order for plotting
    shuffle = FALSE,        # Option to shuffle data points to randomize plotting order
    seed = 1,               # Random seed for shuffling
    label = FALSE,          # Option to add labels to clusters
    label.size = 4,         # Size of cluster labels
    label.color = 'black',  # Color of the labels
    label.box = FALSE,      # Option to put labels in a box
    repel = FALSE,          # Use repelling labels to avoid overlap
    cells.highlight = NULL, # Subset of cells to highlight
    cols.highlight = '#DE2D26', # Color for highlighted cells
    sizes.highlight = 1,    # Size of highlighted cells
    na.value = 'grey50',    # Color for NA values
    ncol = NULL,            # Number of columns for faceting
    combine = TRUE,         # Combine plots into one
    raster = NULL,          # Option to use raster graphics
    raster.dpi = c(512, 512) # Resolution for raster graphics
) {
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  colnames(data) <- paste0("UMAP", dims)
  data <- as.data.frame(x = data)
  dims <- paste0("UMAP", dims)
  data <- cbind(data, group.by)
  orig.groups <- group.by
  group.by <- colnames(x = data)[3:ncol(x = data)]
  
  # Detect if `group.by` contains numeric data for continuous coloring
  is_continuous <- is.numeric(data[[group.by]])
  
  # Convert group.by to factor if it’s not numeric (categorical coloring)
  if (!is_continuous) {
    data[[group.by]] <- factor(data[[group.by]])
  }
  
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  
  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      plot <- SingleDimPlot(
        data = data[, c(dims, x, split.by, shape.by)],
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )
      
      # Apply continuous color scale if `group.by` is numeric
      if (is_continuous) {
        plot <- plot + scale_color_viridis_c(option = "D")  # Use viridis color scale
      }
      
      if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = x,
          repel = repel,
          size = label.size,
          split.by = split.by,
          box = label.box,
          color = label.color
        )
      }
      
      if (!is.null(x = split.by)) {
        plot <- plot + FacetTheme() +
          facet_wrap(
            facets = vars(!!sym(x = split.by)),
            ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
              length(x = unique(x = data[, split.by]))
            } else {
              ncol
            }
          )
      }
      
      plot <- plot + labs(title = NULL)
    }
  )
  
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  
  return(plots)
}


# create_elbow_plot: Generates an elbow plot using the top variable features in a Seurat object 
# based on standardized variance, helping to identify the most informative genes for downstream 
# analyses. The elbow plot shows the top `k` genes by variance, allowing users to identify an 
# optimal number of genes based on the point of inflection in the variance distribution.


create_elbow_plot <- function(object, assay = "RNA", k = 5000) {
  # Get the high variable feature information
  hvf.info <- HVFInfo(object = object, assay = assay)
  
  # Sort the features by standardized variance in descending order
  sorted_standardized_variances <- hvf.info[order(hvf.info$variance.standardized, decreasing = TRUE), ]
  
  # Convert row names to a column named "Gene"
  sorted_standardized_variances <- sorted_standardized_variances %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(Gene_Index = row_number())
  
  # Select the top k genes based on variance
  top_interest <- head(sorted_standardized_variances, k)
  
  # Create the elbow plot
  ggplot(top_interest, aes(x = Gene_Index, y = variance.standardized)) +
    geom_line() +
    labs(title = paste("Elbow Plot of Top", k, "Genes by Standardized Variance"),
         x = paste("Top", k, "Genes (ordered by variance)"),
         y = "Standardized Variance") +
    theme_minimal()
}

