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