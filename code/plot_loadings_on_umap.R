#' Plot loadings on UMAP for a given factor
#'
#' @param umap a matrix of UMAP embeddings
#' @param loading a vector of loadings
#' @param factor_num an integer or string indicating the factor number
#' @param size size of the points
#'
#' @return gg object
#' @export
#'
#' @examples
#' plt <- plot_loadings_on_umap(Embeddings(data$UMAP), ldf(flashier_fit_semi, type="m"), 1)
#' print(plt)
plot_loadings_on_umap <- function(umap, loading, factor_num, size = 3) {
  gplot <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = size, aes(alpha = loading, color = loading), show.legend = FALSE) +
    labs(title = paste("factor ", factor_num) , x = "UMAP 1", y = "UMAP 2") +
    theme_minimal() +
    scale_color_gradient2(low="#f0f0f0", mid="#bdbdbd", high= "black") ##636363")
  return(gplot)
}
