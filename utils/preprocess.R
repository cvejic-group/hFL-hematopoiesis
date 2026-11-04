# notin
# from: https://stackoverflow.com/questions/38351820/negation-of-in-in-r
# https://stackoverflow.com/questions/5831794/opposite-of-in-exclude-rows-with-values-specified-in-a-vector
`%notin%` <- Negate(`%in%`)
`%!in%` <- Negate(`%in%`)
`%nin%` <- Negate(`%in%`)
`%ni%` <- Negate(`%in%`)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# determine the number of PCs
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# 1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2. The point where the percent change in variation between the consecutive PCs is less than 0.1%.
# 3. choose the smaller one
determine_pc_num <- function(seurat_obj) {
  # Determine percent of variation associated with each PC
  pct <- seurat_obj[["pca"]]@stdev^2 / sum(seurat_obj[["pca"]]@stdev^2) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # choose the smaller one
  return(min(c(co1, co2)))
}

# https://pegasus.readthedocs.io/en/stable/_modules/pegasus/plotting/plot_library.html#elbowplot
##################################
#                                #
#        Selection number        #
#           components           #
#                                #
##################################

## ---------------------------------------------##
## ------------------FUNCTIONS------------------##
## ---------------------------------------------##

#' our_largest_variance_from_random_matrix
#' @description
#' Calculates the largest variance of a random matrix based on Johnstone 2001
#' (https://projecteuclid.org/journals/annals-of-statistics/volume-29/issue-2/On-the-distribution-of-the-largest-eigenvalue-in-principal/10.1214/aos/1009210544.full)
#' and Shekhar et al. 2022(https://elifesciences.org/articles/73809)
#' @param ncells number of cells
#' @param nfeatures number of features
#' @param pval pvalue cutoff (either 0.01 or 0.05)
#' @return variance treshold
largest_variance_from_random_matrix <- function(ncells = NULL, nfeatures = NULL, pval = .05) {
  # Check function inputs
  if (!pval %in% c(0.01, 0.05, "0.01", "0.05")) {
    warning("pval was not 0.01 or 0.05, was by default assigned to 0.01")
    pval <- "0.01"
  }
  # Quantiles from the Tracy-Widom distribution of order 1.
  quantiles <- c("0.01" = 2.023335, "0.05" = 0.9792895)

  # Calculate treshold depending on pval
  val1 <- (ncells - 1)**0.5
  val2 <- nfeatures**0.5
  mu <- (val1 + val2)**2
  sigma <- (val1 + val2) * (1.0 / val1 + 1.0 / val2)**(1.0 / 3.0)
  res <- (quantiles[as.character(pval)] * sigma + mu) / (ncells - 1)

  return(res)
}

#' ncomps_signi_pegasus
#' @description
#' Generate ElbowPlot and suggest number of components to select based on random matrix theory as in Pegasus
#' @param seu seurat object
#' @param reduction type of reduction to assess, default "pca"
#' @param pvalue pvalue cutoff to define variance cutoff (0.01 or 0.05), default 0.05, see largest_variance_from_random_matrix()
#' @return ncomps
ncomps_signi_pegasus <- function(seu = NULL, reduction = "pca", pvalue = 0.05) {
  # Calculate variance treshold
  tresh <- largest_variance_from_random_matrix(
    ncells = ncol(seu),
    nfeatures = length(VariableFeatures(seu)),
    pval = pvalue
  )
  # Variance of each comp
  var <- (seu@reductions[[reduction]]@stdev)^2 ## eigen values

  # Components with variance above treshold
  ncomps <- sum(var > tresh)

  # cat(ncomps, "components are recommended over the", length(var), "tested", "\n")

  # ElbowPlot
  # df <- data.frame(var = var, cumsum = cumsum(var),
  #                 rank = 1:length(var))
  # print(ggplot(df, aes(cumsum, var, label = rank, color = rank > ncomps)) +
  #        geom_text(size = 3) + theme_bw() + scale_color_manual(values = c("#e7298a", "#666666")))

  return(ncomps)
}
