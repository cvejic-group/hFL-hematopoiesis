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


# rm doublets, using scDblFinder
remove_doublet <- function(srat_obj = NULL) {
  srat_obj <- as.Seurat(scDblFinder(as.SingleCellExperiment(srat_obj)))
  doublet_barcodes <- names(srat_obj$scDblFinder.class)[srat_obj$scDblFinder.class == "doublet"]
  srat_filt <- subset(srat_obj, subset = scDblFinder.class == "singlet")
  # reasons
  reasons <- data.frame(
    barcode = colnames(srat_obj)
  ) |>
    mutate(doublet = case_when(
      barcode %in% doublet_barcodes ~ TRUE,
      TRUE ~ FALSE
    ))
  srat_filt@misc$reasons <- reasons
  return(srat_filt)
}


# filter low-quality cells
# with static threshold
fast_qc <- function(srat_obj = NULL) {
  reasons <- data.frame(
    barcode = colnames(srat_obj),
    low_gene = srat_obj$nFeature_RNA < 250,
    high_gene = srat_obj$nFeature_RNA > 8000,
    low_umi = srat_obj$nCount_RNA < 750,
    high_umi = srat_obj$nCount_RNA > 110000,
    high_mt = srat_obj$percent.mt > 20
  )
  reasons$discard <- reasons$low_gene | reasons$high_gene |
    reasons$low_umi | reasons$high_umi |
    reasons$high_mt
  srat_filt <- srat_obj[, !reasons$discard]
  srat_filt@misc$reasons <- reasons
  return(srat_filt)
}

# filter low-quality cells
# with dynamic and static threshold
lowQuality_filter <- function(seurat_obj = NULL,
                              remove_cell = TRUE,
                              gene_fixed_low = 250,
                              gene_fixed_high = 10000,
                              umi_fixed_low = 500,
                              umi_fixed_high = 100000,
                              mt_fixed_high = NULL,
                              ...) {
  # filter low quality
  filt_df <- scuttle::quickPerCellQC(seurat_obj@meta.data,
    sum.field = "nCount_RNA",
    detected.field = "nFeature_RNA",
    sub.fields = c(
      "percent.mt",
      "percent.ribo"
    ),
    ...
  ) |>
    as.data.frame() |>
    mutate(barcode = colnames(seurat_obj)) |>
    mutate(high_n_features = FALSE) %>% # add one column for gene_fixed_high filter
    mutate(high_lib_size = FALSE) %>% # add one column for umi_fixed_high filter
    dplyr::select(
      barcode,
      low_lib_size, high_lib_size,
      low_n_features, high_n_features,
      high_percent.mt, high_percent.ribo
    ) # change column order

  # further filter by fixed threshold
  if (!is.null(gene_fixed_low)) {
    filt_df$low_n_features <- filt_df$low_n_features | (seurat_obj$nFeature_RNA < gene_fixed_low)
  }
  if (!is.null(gene_fixed_high)) {
    filt_df$high_n_features <- seurat_obj$nFeature_RNA > gene_fixed_high
  }
  if (!is.null(umi_fixed_low)) {
    filt_df$low_lib_size <- filt_df$low_lib_size | (seurat_obj$nCount_RNA < umi_fixed_low)
  }
  if (!is.null(umi_fixed_high)) {
    filt_df$high_lib_size <- seurat_obj$nCount_RNA > umi_fixed_high
  }
  if (!is.null(mt_fixed_high)) {
    filt_df$high_percent.mt <- filt_df$high_percent.mt | (seurat_obj$percent.mt > mt_fixed_high)
  }
  # discard
  reasons <- seurat_obj@misc$reasons |>
    left_join(filt_df, by = "barcode") |>
    mutate(discard = doublet | low_lib_size | high_lib_size |
      low_n_features | high_n_features |
      high_percent.mt | high_percent.ribo) |>
    mutate(discard = replace_na(discard, FALSE))
  # filter?
  if (remove_cell) {
    srat_filt <- seurat_obj[, !reasons[!reasons$doublet, ]$discard]
    srat_filt@misc$reasons <- reasons
    return(srat_filt)
  } else {
    seurat_obj@misc$reasons <- reasons
    return(seurat_obj)
  }
}

# fast clustering using default parameters
fast_cluster <- function(seurat_obj = NULL,
                         vars.to.regress = NULL,
                         res = .8) {
  .seurat_obj <- NormalizeData(seurat_obj, verbose = F)
  .seurat_obj <- FindVariableFeatures(.seurat_obj, verbose = F)
  .seurat_obj <- ScaleData(.seurat_obj,
    vars.to.regress = vars.to.regress,
    verbose = F
  )
  # remove existing
  if ("pca" %in% names(.seurat_obj@reductions)) {
    .seurat_obj[["pca"]] <- NULL
  }
  if ("PCA" %in% names(.seurat_obj@reductions)) {
    .seurat_obj[["PCA"]] <- NULL
  }
  .seurat_obj <- RunPCA(.seurat_obj,
    npcs = 100, verbose = F
  )
  pc_num <- ncomps_signi_pegasus(seu = .seurat_obj)
  .seurat_obj@misc$ncomps <- pc_num
  .seurat_obj <- FindNeighbors(.seurat_obj, dims = 1:pc_num, verbose = F)
  .seurat_obj <- FindClusters(.seurat_obj, resolution = res, verbose = F)
  # remove existing
  if ("umap" %in% names(.seurat_obj@reductions)) {
    .seurat_obj[["umap"]] <- NULL
  }
  if ("UMAP" %in% names(.seurat_obj@reductions)) {
    .seurat_obj[["UMAP"]] <- NULL
  }
  .seurat_obj <- RunUMAP(.seurat_obj, dims = 1:pc_num, verbose = F)
  return(.seurat_obj)
}

# modified from RIDDLER
# Turn bed reads from all cells into window x cell matrix
window_matrix <- function(bed_file = NULL, w_file = NULL, out_file = NULL, n_frag = 1000) {
  library(data.table)
  print(paste0("Reading ", bed_file))
  print(paste0("Write output to ", out_file))
  chr_list <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY")

  # Load fragments, get barcodes with at least 1000 reads
  frags <- fread(bed_file, data.table = FALSE)
  print("Fragments loaded")
  b_tab <- table(frags[, 4])
  barcodes <- names(b_tab[b_tab >= n_frag])
  frags <- frags[frags[, 4] %in% barcodes, ]
  print(paste0("Barcodes filtered: ", length(barcodes)))
  print(paste0("Fragments to read: ", nrow(frags)))
  # Read in window file
  window <- read.table(w_file, header = FALSE, stringsAsFactors = FALSE, )
  n <- nrow(window)
  window <- cbind(window, 1:n)

  hash_window <- list()
  for (i in 1:length(chr_list)) {
    hash_window[[i]] <- window[window[, 1] == chr_list[i], ]
  }
  names(hash_window) <- chr_list

  # Create matrix of total reads per window
  counts <- matrix(0, nrow = n, ncol = length(barcodes))
  colnames(counts) <- barcodes

  t <- proc.time()[3]
  for (chr in chr_list) {
    frags_c <- frags[frags[1] == chr, ]
    chr_wnd <- hash_window[[chr]]
    # quickly find intersecting windows
    s_ids <- findInterval(frags_c[, 2], chr_wnd[, 2])
    e_ids <- findInterval(frags_c[, 3], chr_wnd[, 3]) + 1
    c_ids <- which(s_ids != e_ids)
    if (length(c_ids) > 0) {
      # find best overlapping window
      offset <- apply(abs(frags_c[c_ids, c(2, 3)] - chr_wnd[s_ids[c_ids], 3]), 1, which.max) - 1
      s_ids[c_ids] <- s_ids[c_ids] + offset
    }
    # Adjust chr ids to global ids
    w_ids <- chr_wnd[s_ids, 5]
    # Add to counts
    tab_frag <- table(w_ids, frags_c[, 4])
    counts[as.numeric(rownames(tab_frag)), colnames(tab_frag)] <- tab_frag

    print(paste(chr, proc.time()[3] - t))
    t <- proc.time()[3]
  }

  # save sums
  write.table(cbind(window[, 1:4], counts), file = out_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# filter out features that not in seurat_obj
filter_features <- function(seurat_obj = NULL, features = NULL) {
  tmp <- unique(features)
  features <- c()
  for (i in tmp) {
    if (i %in% rownames(seurat_obj)) {
      features <- c(features, i)
    }
  }
  return(features)
}

# harmony by groups
# modified from: https://github.com/kaizhang/SnapATAC2/blob/main/snapatac2-python/python/snapatac2/preprocessing/_harmony.py
harmony_by_group <- function(
    srat_obj = NULL,
    batch_keys = NULL,
    use_rep = "pca",
    use_dims = NULL,
    groupby = NULL,
    key_added = NULL,
    plot_convergence = TRUE,
    ...) {
  if (!is.null(groupby)) {
    # take out pca
    mat <- srat[[use_rep]]@cell.embeddings
    if (!is.null(use_dims)) {
      mat <- mat[, 1:use_dims]
    }
    if (is.character(groupby)) {
      groupby <- srat_obj@meta.data[[groupby]]
    }
    groups <- unique(groupby)
    for (group in groups) {
      group_idx <- which(groupby == group)
      meta_data <- srat_obj@meta.data[group_idx, ]
      mat[group_idx, ] <- RunHarmony(mat[group_idx, ], meta_data, batch_keys,
        plot_convergence = plot_convergence, ...
      )
      srat_obj[[key_added]] <- CreateDimReducObject(embeddings = mat, assay = DefaultAssay(srat_obj))
    }
  } else {
    srat_obj <- RunHarmony(srat_obj,
      group.by.vars = batch_keys,
      dims.use = 1:use_dims,
      plot_convergence = plot_convergence,
      reduction.save = key_added,
      ...
    )
  }
  return(srat_obj)
}
