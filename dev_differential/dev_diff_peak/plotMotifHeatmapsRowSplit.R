# Standalone monaLisa heatmap function with row_split support.
# Requires an installed monaLisa package and its normal dependencies.
# source("plotMotifHeatmapsRowSplit.R"), then use plotMotifHeatmapsRowSplit(...).
# Group vectors/data-frame rows and hclust leaves use the input motif order.
# Numeric cutree groups matching a tree cut preserve the full dendrogram.
# Other groups preserve the supplied tree separately within each slice.
# Use factor levels for custom slice order. A scalar split requires clustering.
# Returns the usual heatmap list; doPlot = FALSE suppresses drawing.

plotMotifHeatmapsRowSplit <- local({
  validateSplit <- function(row_split, n) {
    if (is.null(row_split)) return(invisible(NULL))
    if (is.numeric(row_split) && length(row_split) == 1L &&
        is.null(dim(row_split))) {
      if (!is.finite(row_split) || row_split != floor(row_split) ||
          row_split < 1 || row_split > n) {
        stop("'row_split' must be an integer between 1 and the number of motifs.")
      }
    } else if (is.data.frame(row_split)) {
      if (nrow(row_split) != n || ncol(row_split) == 0L ||
          !all(vapply(row_split, function(z) is.atomic(z) && is.null(dim(z)),
                      logical(1)))) {
        stop("'row_split' must have one row per motif and atomic grouping columns.")
      }
    } else if (!is.atomic(row_split) || !is.null(dim(row_split)) ||
               length(row_split) != n) {
      stop("'row_split' must be a grouping vector with one entry per motif, a data frame, or a single number of slices.")
    }
    invisible(NULL)
  }

  prepareSplit <- function(row_split, cluster, n) {
    validateSplit(row_split, n)
    if (is(cluster, "hclust") && length(cluster$order) != n) {
      stop("'cluster' must contain exactly one leaf per motif, in input row order.")
    }
    result <- list(row_split = row_split, cluster_rows = cluster)
    if (is.null(row_split)) return(result)
    if (is.numeric(row_split) && length(row_split) == 1L) {
      if (identical(cluster, FALSE)) {
        stop("A numeric 'row_split' requires 'cluster = TRUE' or an hclust object.")
      }
      return(result)
    }
    if (identical(cluster, FALSE)) return(result)

    # A numeric cutree vector can use ComplexHeatmap's native tree cutting.
    if (is.numeric(row_split) && !anyNA(row_split)) {
      k <- length(unique(row_split))
      tree_groups <- stats::cutree(cluster, k = k)
      if (identical(match(row_split, row_split), match(tree_groups, tree_groups))) {
        result$row_split <- k
        return(result)
      }
    }

    # For other groupings, retain the original tree within each slice.
    tree <- stats::as.dendrogram(cluster)
    prune <- function(node, selected) {
      if (stats::is.leaf(node)) {
        original <- as.integer(node)
        index <- match(original, selected)
        if (is.na(index)) return(NULL)
        node[] <- index
        attr(node, "label") <- as.character(original)
        return(node)
      }
      children <- Filter(Negate(is.null), lapply(node, prune, selected = selected))
      if (length(children) == 0L) return(NULL)
      if (length(children) == 1L) return(children[[1]])
      node[] <- children
      members <- vapply(children, function(z) attr(z, "members"), numeric(1))
      midpoints <- vapply(children, function(z) {
        if (stats::is.leaf(z)) 0 else attr(z, "midpoint")
      }, numeric(1))
      attr(node, "members") <- sum(members)
      offsets <- c(0, cumsum(members)[-length(members)])
      attr(node, "midpoint") <- mean(offsets + midpoints)
      node
    }
    result$cluster_rows <- function(mat) prune(tree, as.integer(rownames(mat)))
    result
  }

  function(x,
           which.plots = c("negLog10P", "pearsonResid",
                           "negLog10Padj", "log2enr"),
           width = 4,
           col.enr = c("#053061", "#2166AC", "#4393C3",
                       "#92C5DE", "#D1E5F0", "#F7F7F7",
                       "#FDDBC7", "#F4A582", "#D6604D",
                       "#B2182B", "#67001F"),
           col.sig = c("#F0F0F0", "#D9D9D9", "#BDBDBD",
                       "#969696", "#737373", "#525252",
                       "#252525", "#000000"),
           col.gc = c("#F7FCF5", "#E5F5E0", "#C7E9C0",
                      "#A1D99B", "#74C476", "#41AB5D",
                      "#238B45", "#006D2C", "#00441B"),
           maxEnr = NULL,
           maxSig = NULL,
           highlight = NULL,
           cluster = FALSE,
           show_dendrogram = FALSE,
           show_motif_GC = FALSE,
           show_seqlogo = FALSE,
           show_bin_legend = FALSE,
           width.seqlogo = 1.5,
           use_raster = FALSE,
           na_col = "white",
           doPlot = TRUE,
           row_split = NULL,
           ...) {
    stopifnot(exprs = {
      is(x, "SummarizedExperiment")
      all(which.plots %in% assayNames(x))
      "bins" %in% names(metadata(x))
      (!show_motif_GC || "motif.percentGC" %in% colnames(rowData(x)))
    })
    b <- metadata(x)$bins
    .assertScalar(x = width, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = show_dendrogram, type = "logical")
    .assertScalar(x = show_motif_GC, type = "logical")
    .assertScalar(x = show_seqlogo, type = "logical")
    .assertScalar(x = show_bin_legend, type = "logical")
    .assertScalar(x = width.seqlogo, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = use_raster, type = "logical")
    .assertScalar(x = na_col, type = "character")
    .assertScalar(x = doPlot, type = "logical")
    stopifnot(exprs = {
      ncol(x) == nlevels(b)
      all(which.plots %in% c("negLog10P", "negLog10Padj",
                             "pearsonResid", "log2enr"))
      is.null(highlight) || (is.logical(highlight) &&
                               length(highlight) == nrow(x))
    })
    validateSplit(row_split, nrow(x))
    bincols <- attr(getColsByBin(b), "cols")
    if (identical(cluster, TRUE)) {
      clAssayName <- "pearsonResid"
      clAssay <- assay(x, clAssayName)
      allNA <- rowSums(is.na(clAssay)) == ncol(clAssay)
      if (any(allNA)) {
        warning("removing motifs without finite values in '",
                clAssayName, "': ",
                paste(rownames(clAssay)[allNA], collapse = ", "))
        x <- x[!allNA, ]
        clAssay <- clAssay[!allNA, , drop = FALSE]
        if (!is.null(highlight)) highlight <- highlight[!allNA]
        if (is.data.frame(row_split)) {
          row_split <- row_split[!allNA, , drop = FALSE]
        } else if (!is.null(row_split) && length(row_split) > 1L) {
          row_split <- row_split[!allNA]
        }
      }
      clres <- hclust(dist(clAssay))
    } else if (identical(cluster, FALSE)) {
      clres <- FALSE
    } else if (is(cluster, "hclust")) {
      clres <- cluster
    } else {
      cli_abort("{.arg cluster} must be either {.code TRUE}, {.code FALSE} or an {.cls hclust}-object.")
    }
    rowLayout <- prepareSplit(row_split, clres, nrow(x))
    hmBin <- HeatmapAnnotation(
      df = data.frame(bin = factor(colnames(x),
                                   levels = colnames(x))),
      name = "bin",
      col = list(bin = bincols),
      show_annotation_name = FALSE,
      which = "column", width = unit(width,"inch"),
      annotation_height = unit(width / 16, "inch"),
      show_legend = show_bin_legend)
    tmp <- matrix(if (!is.null(highlight)) {
      as.character(highlight)
    } else {
      rep(NA, nrow(x))
    },
    ncol = 1, dimnames = list(unname(rowData(x)$motif.name), NULL))
    labelArgs <- list(...)
    if (is.function(rowLayout$cluster_rows)) {
      # Internal keys handle duplicated motif display names.
      rownames(tmp) <- as.character(seq_len(nrow(x)))
      if (!"row_labels" %in% names(labelArgs))
        labelArgs$row_labels <- unname(rowData(x)$motif.name)
      if (!"cluster_row_slices" %in% names(labelArgs))
        labelArgs$cluster_row_slices <- FALSE
      if (!"row_dend_reorder" %in% names(labelArgs))
        labelArgs$row_dend_reorder <- FALSE
    }
    hmSeqlogo <- NULL
    if (show_seqlogo) {
      pfms <- rowData(x)$motif.pfm
      maxwidth <- max(vapply(TFBSTools::Matrix(pfms), ncol, 0L))
      grobL <- lapply(pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")
      hmSeqlogo <- HeatmapAnnotation(
        logo = annoSeqlogo(grobL = grobL, which = "row",
                           space = unit(0.5, "mm"),
                           width = unit(width.seqlogo, "inch")),
        show_legend = FALSE, show_annotation_name = FALSE, which = "row")
    }
    hmMotifs <- do.call(Heatmap, c(list(
      matrix = tmp, name = "names",
      width = unit(if (!is.null(highlight)) .2 else 0, "inch"),
      na_col = NA, col = c("TRUE" = "green3", "FALSE" = "white"),
      cluster_rows = rowLayout$cluster_rows, row_split = rowLayout$row_split,
      show_row_dend = show_dendrogram,
      cluster_columns = FALSE, show_row_names = TRUE,
      row_names_side = "left", show_column_names = FALSE,
      show_heatmap_legend = FALSE, left_annotation = hmSeqlogo
    ), labelArgs))

    assayNameMap1 <- c(negLog10P = "P value",
                       negLog10Padj = "adj. P value",
                       pearsonResid = "Pearson residual",
                       log2enr = "log2 enrichment")
    assayNameMap2 <- c(negLog10P = "P value (-log10)",
                       negLog10Padj = "adj. P value (-log10)",
                       pearsonResid = "Pearson residual (o-e)/sqrt(e)",
                       log2enr = "enrichment (log2)")
    L <- list(labels = hmMotifs)
    if (show_motif_GC) {
      tmp <- as.matrix(rowData(x)[, "motif.percentGC", drop = FALSE])
      hmPercentGC <- Heatmap(
        matrix = tmp, name = "Percent G+C",
        width = unit(0.2, "inch"), na_col = NA,
        col = colorRamp2(breaks = c(0, seq(20, 80, length.out = 254), 100),
                         colors = colorRampPalette(col.gc)(256)),
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(color_bar = "continuous"),
        use_raster = use_raster,
        ...
      )
      L <- c(L, list("percentGC" = hmPercentGC))
    }
    ret <- c(L, lapply(which.plots, function(w) {
      dat <- assay(x, w)
      if ((w == "pearsonResid") | (w == "log2enr")) {
        rng <- c(-1, 1) * if (is.null(maxEnr)) {
          quantile(abs(dat), .995, na.rm = TRUE)
        } else {
          maxEnr
        }
        cols <- col.enr
      } else {
        rng <- c(0,
                 if (is.null(maxSig)) {
                   quantile(dat, .995, na.rm = TRUE)
                 } else {
                   maxSig
                 })
        cols <- col.sig
      }
      Heatmap(
        matrix = dat,
        name = assayNameMap1[w],
        width = unit(width,"inch"),
        column_title = assayNameMap2[w],
        col = colorRamp2(breaks = seq(rng[1], rng[2], length.out = 256),
                         colors = colorRampPalette(cols)(256)),
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        # column_names_side = "bottom",
        # column_names_max_height = unit(1.5,"inch"),
        top_annotation = hmBin, show_heatmap_legend = TRUE,
        heatmap_legend_param = list(color_bar = "continuous"),
        use_raster = use_raster,
        na_col = na_col,
        ...
      )
    }))
    names(ret)[seq(length(ret) - length(which.plots) + 1L, length(ret))] <-
      which.plots
    if (doPlot) {
      show(Reduce(ComplexHeatmap::add_heatmap, ret))
    }
    invisible(ret)
  }
}, envir = new.env(parent = asNamespace("monaLisa")))
