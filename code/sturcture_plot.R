
embed_method <- function(fit) {
  if (nrow(fit$L) < 20) {
    return(rnorm(nrow(fit$L)))
  } else {
    d <- dim(fit$L)
    message(sprintf("Running tsne on %s x %s matrix.",d[1],d[2]))
    return(drop(suppressMessages(tsne_from_topics(fit,dims = 1))))
  }
}


structure_plot_ggplot_call <- function (dat, colors, ticks = NULL,
                                        font.size = 9, linewidth = 3) {
  ggplot(dat,aes_string(x = "sample",y = "prop",fill = "topic")) +
    geom_col(linewidth = linewidth,width = 0.9) +
    scale_fill_manual(values = colors) +
    labs(x = "",y = "membership") +
    theme_cowplot(font.size) +
    theme(axis.line   = element_blank(),
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1))
}


compile_structure_plot_data <- function (L, topics) {
  n <- nrow(L)
  k <- length(topics)
  if (is.null(rownames(L))) {
    rownames(L) <- traits  # Fallback: use numeric indices
  }
  dat <- data.frame(sample = rownames(L),
                    topic  = rep(topics,each = n),
                    prop   = c(L[,topics]))
  dat$topic <- factor(dat$topic,topics)
  return(dat)
}




fit <- fl_snmf_ldf$F
topics <- trait_factors
grouping <- traits
gap <- 3
#perplexity = 70,
n <- 500
loadings_order <- "embed"

# Check if 'F' is a matrix or model object
if (is.matrix(fit)) {
  L <- fit
  F <- matrix(1, nrow(L), ncol(L))
  fit <- list(F = F, L = L)
  class(fit) <- "poisson_nmf_fit"
} else {
  if (!(inherits(F, "poisson_nmf_fit") | inherits(F, "multinom_topic_model_fit")))
    stop("Invalid input class.")
  if (inherits(F, "poisson_nmf_fit"))
    fit <- poisson2multinom(F)  # Assume this function exists
}

n0 <- nrow(fit$L)  # Number of samples
k  <- ncol(fit$L)   # Number of topics

if (is.null(colnames(fit$L)))
  colnames(fit$L) <- paste0("k", 1:k)

if (is.null(topics))
  topics <- order(colMeans(fit$L))
if (!is.character(topics))
  topics <- colnames(fit$L[,topics,drop = FALSE])
if (!(length(topics) > 1 & all(is.element(topics,colnames(fit$L)))))
  stop("Input argument \"topics\" should be a subset of at least two ",
       "topics (columns of fit$L) specified by their names or column ",
       "indices")

if (is.null(grouping))
  grouping <- factor(rep(1,n0))
if (!is.factor(grouping))
  grouping <- as.factor(grouping)
if (length(grouping) != n0)
  stop("Input argument \"grouping\" should be a factor with one entry ",
       "for each row of fit$L")

colors9 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
             "#ffff33","#a65628","#f781bf","#999999")


if (k < 10) {
  colors <- colors9
} else if (k < 22) {
    colors <- kelly()[-1]
} else {
    colors <- glasbey()[-1]
}

if (length(colors) < k)
  stop("There must be at least as many colours as topics")

names(colors) <- colnames(fit$L)
colors <- colors[topics]


if (all(loadings_order == "embed")) {

  # If necessary, randomly subsample the rows of L.
  if (n < n0) {
    rows <- sample(n0,n)
    fit <- select_loadings(fit,rows)
    grouping <- grouping[rows,drop = FALSE]
  }

  # The ordering of the rows is not provided, so determine an
  # ordering by computing a 1-d embedding of L.
  if (nlevels(grouping) == 1) {
    y <- embed_method(fit)
    loadings_order <- order(y)
  } else {
    loadings_order <- NULL
    for (group in levels(grouping)) {
      i <- which(grouping == group)
      if (length(i) > 0)
        y <- embed_method(select_loadings(fit,i))
      loadings_order <- c(loadings_order,i[order(y)])
    }
  }
} else {
  if (!is.null(n))
    warning("Input argument \"n\" is ignored when \"loadings_order\" is ",
            "not \"embed\"")
  if (is.character(loadings_order))
    loadings_order <- match(loadings_order,rownames(fit$L))
}

fit$L <- fit$L[loadings_order,, drop = FALSE]
grouping <- grouping[loadings_order,drop = TRUE]

dat <- compile_structure_plot_data(fit$L,topics)

p1 <- structure_plot_ggplot_call(dat, colors)



fit <- fl_mf_ldf$F
topics <- trait_factors
grouping <- NULL
gap <- 3
#perplexity = 70,
n <- 500
loadings_order <- "embed"

# Check if 'F' is a matrix or model object
if (is.matrix(fit)) {
  L <- fit
  F <- matrix(1, nrow(L), ncol(L))
  fit <- list(F = F, L = L)
  class(fit) <- "poisson_nmf_fit"
} else {
  if (!(inherits(F, "poisson_nmf_fit") | inherits(F, "multinom_topic_model_fit")))
    stop("Invalid input class.")
  if (inherits(F, "poisson_nmf_fit"))
    fit <- poisson2multinom(F)  # Assume this function exists
}

n0 <- nrow(fit$L)  # Number of samples
k  <- ncol(fit$L)   # Number of topics

if (is.null(colnames(fit$L)))
  colnames(fit$L) <- paste0("k", 1:k)

if (is.null(topics))
  topics <- order(colMeans(fit$L))
if (!is.character(topics))
  topics <- colnames(fit$L[,topics,drop = FALSE])
if (!(length(topics) > 1 & all(is.element(topics,colnames(fit$L)))))
  stop("Input argument \"topics\" should be a subset of at least two ",
       "topics (columns of fit$L) specified by their names or column ",
       "indices")

if (is.null(grouping))
  grouping <- factor(rep(1,n0))
if (!is.factor(grouping))
  grouping <- as.factor(grouping)
if (length(grouping) != n0)
  stop("Input argument \"grouping\" should be a factor with one entry ",
       "for each row of fit$L")

colors9 <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
             "#ffff33","#a65628","#f781bf","#999999")


if (k < 10) {
  colors <- colors9
} else if (k < 22) {
  colors <- kelly()[-1]
} else {
  colors <- glasbey()[-1]
}

if (length(colors) < k)
  stop("There must be at least as many colours as topics")

names(colors) <- colnames(fit$L)
colors <- colors[topics]


if (all(loadings_order == "embed")) {

  # If necessary, randomly subsample the rows of L.
  if (n < n0) {
    rows <- sample(n0,n)
    fit <- select_loadings(fit,rows)
    grouping <- grouping[rows,drop = FALSE]
  }

  # The ordering of the rows is not provided, so determine an
  # ordering by computing a 1-d embedding of L.
  if (nlevels(grouping) == 1) {
    y <- embed_method(fit)
    loadings_order <- order(y)
  } else {
    loadings_order <- NULL
    for (group in levels(grouping)) {
      i <- which(grouping == group)
      if (length(i) > 0)
        y <- embed_method(select_loadings(fit,i))
      loadings_order <- c(loadings_order,i[order(y)])
    }
  }
} else {
  if (!is.null(n))
    warning("Input argument \"n\" is ignored when \"loadings_order\" is ",
            "not \"embed\"")
  if (is.character(loadings_order))
    loadings_order <- match(loadings_order,rownames(fit$L))
}

fit$L <- fit$L[loadings_order,, drop = FALSE]
grouping <- grouping[loadings_order,drop = TRUE]

dat <- compile_structure_plot_data(fit$L,topics)

p2 <- structure_plot_ggplot_call(dat, colors)







