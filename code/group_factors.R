ANOVA_one_factor <- function(l_vec, group_vec, stats = "R2"){
  if(stats %in% c("F value", "Pr(>F)")){
    linear_mod <- lm(l_vec ~ group_vec)
    return(anova(linear_mod)[[stats]][1])
  }
  else if(stats == "R2"){
    linear_mod <- lm(l_vec ~ group_vec)
    return(summary(linear_mod)$r.squared)
  }
  else if(stats == "abs_diff"){
    linear_mod <- lm(l_vec ~ -1 + group_vec)
    coef_vec <- coef(linear_mod)
    return(abs(diff(range(coef_vec))/(mean(coef_vec)-min(coef_vec))))
  }
  else{
    stop("stats should be one of 'F value', 'Pr(>F)', 'R2', 'abs_diff'")
  }
}

ANOVA_factors <- function(L, group_vec, stats = "R2"){
  k <- ncol(L)
  stats_values <- numeric(k)
  group_vec <- as.factor(group_vec)

  for (i in 1:k){
    stats_values[i] <- ANOVA_one_factor(L[,i], group_vec, stats = stats)
  }

  ordered_df <- data.frame(
    rank = 1:k,
    stats = stats_values[order(stats_values, decreasing = TRUE)],
    factor = order(stats_values, decreasing = TRUE)
  )

  return(ordered_df)
}

structure_plot_group <- function(L, group_vec, cutoff = NULL, stats = "R2", group_name = NULL){
  colnames(L) <- paste0("k",1:ncol(L))
  ordered_factor <- ANOVA_factors(L = L, group_vec = group_vec, stats = stats)
  if (is.null(cutoff)){
    # by default, we only keep the top 20% factors
    cutoff <- quantile(ordered_factor$stats, 0.8)
    ordered_factor <- ordered_factor[ordered_factor$stats >= cutoff,]
  } else{
    ordered_factor <- ordered_factor[ordered_factor$stats >= cutoff,]
  }
  L_selected <- L[,ordered_factor$factor]

  p <- fastTopics::structure_plot(L_selected, grouping = group_vec,
                                  gap = 3,perplexity = 70,n = 500)

  if(is.null(group_name)){
    p <- p +
      labs(y = "membership",title = "group-specific factors",
           fill = "factor",color = "factor")
  }else{
    p <- p +
      labs(y = "membership",title = paste0(group_name,"-specific factors"),
           fill = "factor",color = "factor")
  }
  return(p)
}
