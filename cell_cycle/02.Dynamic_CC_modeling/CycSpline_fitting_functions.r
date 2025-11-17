library(mgcv)
library(Matrix)
library(progressr)

# 1) helper function for input preparation ---------------------------
.prepare_inputs <- function(counts, pseudotime, pcw, libsize, eps = 1e-8){
  stopifnot(length(counts) == length(pseudotime),
            length(counts) == length(pcw),
            length(counts) == length(libsize))
  
  data.frame(
    counts      = as.numeric(counts),
    pseudotime  = as.numeric(pseudotime),
    PCW         = factor(pcw),
    log_libsize = log(pmax(as.numeric(libsize), 1e-8))
  )
}

# 2) Single gene fitting ---------------------------
fit_cycspline <- function(counts_df, k=40,knots=6, method = "ML"){
  if (!is.data.frame(counts_df)) {
    stop("Input counts_df must be a data frame.")
  }  
  if (!all(c("counts", "pseudotime", "PCW", "log_libsize") %in% colnames(counts_df))) {
    stop("Input counts_df must contain the following columns: counts, pseudotime, PCW, log_libsize.")
  }
  
  f_full <- counts ~ s(pseudotime, bs = "cc", k = knots) + PCW + offset(log_libsize)
  f_red  <- counts ~ PCW + offset(log_libsize)  # reduced model (no periodic terms)
  
  fit_full <- gam(f_full, family = nb(), data = counts_df,gamma = 1.4, method = method)
  
  ## use same theta in reduced model to make the LRT stable
  theta_hat <- fit_full$family$getTheta(TRUE)
  fit_red   <- gam(f_red, family = nb(theta = theta_hat), data = counts_df, method = method)
  
  ## joint LRT for sin & cos via anova on nested models
  a <- anova(fit_red, fit_full, test = "Chisq")
  p_value <- a$`Pr(>Chi)`[2]  # H0: beta_sin = beta_cos = 0
  
  
  ## make grid prediction
  grid <- data.frame(
    pseudotime = seq(0, 2*pi, length.out = k),
    PCW = 18,
    log_libsize = median(counts_df$log_libsize)
  )
  yhat <- predict(fit_full, newdata = grid, type = "response")
  
  peak_phase <- grid$pseudotime[ which.max(yhat) ]
  eta <- predict(fit_full, newdata = grid, type = "link")
  upper <- mean( eta[ yhat >= quantile(yhat, 0.9) ] )
  lower <- mean( eta[ yhat <= quantile(yhat, 0.1) ] )
  avg_log2fc <- (upper - lower) / log(2)
  
  
  ## return results
  list(yhat = yhat, 
       p_value = p_value, 
       avg_log2fc = avg_log2fc, 
       peak_phase = peak_phase)
}

# 3) sample a subset of genes (rows) ---------------------------
sample_genes <- function(counts_matrix, n = 500, seed = 123) {
  set.seed(seed)
  g <- nrow(counts_matrix)
  sample.int(g, n)
}


# 4) Assess cyclic spline k ---------------------------
assess_cyclic_k <- function(counts_matrix, pseudotime, pcw, libsize,
                            k_grid = c(6, 8, 12, 16),
                            n_genes = 600,
                            max_fail = 0.10,
                            edf_cap_frac = 0.90,
                            method = "ML",
                            seed = 1L) {
  
  stopifnot(ncol(counts_matrix) == length(pseudotime),
            length(pseudotime) == length(pcw),
            length(pcw) == length(libsize))
  
  set.seed(seed)
  gi <- sample_genes(counts_matrix, n_genes)
  pt_range <- range(as.numeric(pseudotime), finite = TRUE)
  
  mkdf <- function(i) .prepare_inputs(
    counts     = as.numeric(counts_matrix[i, ]),
    pseudotime = as.numeric(pseudotime),
    pcw        = pcw,
    libsize    = libsize
  )
  
  diag <- NULL
  
  options(progressr.enable = TRUE, progressr.clear = FALSE)
  handlers(handler_txtprogressbar())
  
  progressr::with_progress(enable = TRUE,{
    p <- progressr::progressor(steps = length(k_grid) * length(gi),
                               message = "Assessing cyclic k")
    
    diag <- lapply(k_grid, function(k_smooth) {
      fails <- nearcap <- 0L; nfit <- 0L
      for (ix in seq_along(gi)) {
        i  <- gi[ix]
        df <- mkdf(i)
        fit <- try(mgcv::gam(
          counts ~ s(pseudotime, bs = "cc", k = k_smooth) + PCW + offset(log_libsize),
          family = mgcv::nb(), data = df, method = method,
          knots  = list(pseudotime = pt_range)
        ), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          nfit <- nfit + 1L
          edf1 <- as.numeric(summary(fit)$s.table[1, "edf"])
          if (!is.na(edf1) && edf1 > edf_cap_frac * (k_smooth - 1)) nearcap <- nearcap + 1L
          kc <- try(mgcv::k.check(fit, subsample = min(5000, nrow(df)), n.rep = 200), silent = TRUE)
          if (!inherits(kc, "try-error")) {
            kidx <- suppressWarnings(as.numeric(kc[1, "k-index"]))
            pval <- suppressWarnings(as.numeric(kc[1, "p-value"]))
            if (!is.na(pval) && !is.na(kidx) && pval < 0.05 && kidx < 1) fails <- fails + 1L
          }
        }
        p(sprintf("k=%d | gene %d/%d", k_smooth, ix, length(gi)))  # <-- PROGRESS TICK
      }
      data.frame(
        k_smooth      = k_smooth,
        n_fit         = nfit,
        fail_rate     = if (nfit) fails / nfit else NA_real_,
        near_cap_rate = if (nfit) nearcap / nfit else NA_real_
      )
    })
  })
  
  diag <- do.call(rbind, diag)
  ok <- with(diag, !is.na(fail_rate) & fail_rate <= max_fail & near_cap_rate <= max_fail)
  chosen <- if (any(ok)) min(diag$k_smooth[ok]) else max(k_grid)
  list(chosen_k = chosen, diagnostics = diag)
}




# 5) wrap up function for all genes ---------------------------
fit_all <- function(counts_matrix,
                    pseudotime,
                    pcw,
                    libsize,
                    knots=10,
                    k = 40,
                    workers = max(1, parallel::detectCores() - 1),
                    method = "ML",
                    block_size = 200,          # genes per task
                    globals_max_gb = 8) {      # cap for big objects
  
  stopifnot(ncol(counts_matrix) == length(pseudotime),
            length(pseudotime) == length(pcw),
            length(pcw) == length(libsize))
  
  # 1) Allow big globals (needed on Windows / multisession)
  options(future.globals.maxSize = globals_max_gb * 1024^3)
  
  # 2) Pick a plan that avoids copies when possible
  os <- tolower(Sys.info()[["sysname"]])
  if (os == "windows") {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::multicore, workers = workers)  # fast: shared memory
  }
  
  # 3) Make row slicing fast (once)
  if (!inherits(counts_matrix, "dgRMatrix")) {
    counts_r <- methods::as(counts_matrix, "RsparseMatrix")
  } else {
    counts_r <- counts_matrix
  }
  
  G <- nrow(counts_r)
  gene_ids <- rownames(counts_r)
  
  # Preallocate outputs
  yhats <- matrix(NA_real_, nrow = G, ncol = k,
                  dimnames = list(gene_ids, NULL))
  gene_fit <- data.frame(
    p_value     = rep(NA_real_, G),
    avg_log2fc = rep(NA_real_, G),
    peak_phase  = rep(NA_real_, G),
    row.names   = gene_ids
  )
  
  # 4) Indices in blocks to amortize overhead
  idx_blocks <- split(seq_len(G), ceiling(seq_len(G) / block_size))
  
  options(
    progressr.enable  = TRUE,   # allow in non-interactive sessions
    progressr.clear   = FALSE,  # don't erase progress output (important for logs)
    progressr.interval = 0,     # flush updates promptly
    progressr.delay    = 0
  )
  
  progressr::handlers(handler_txtprogressbar())
  progressr::handlers(global = TRUE)
  progressr::with_progress({
    p <- progressr::progressor(steps = length(idx_blocks))
    
    # One future per block (hundreds of genes)
    outs <- future.apply::future_lapply(
      idx_blocks,
      function(ii) {
        yh_loc <- matrix(NA_real_, nrow = length(ii), ncol = k)
        gf_loc <- matrix(NA_real_, nrow = length(ii), ncol = 3)
        
        for (jj in seq_along(ii)) {
          i <- ii[jj]
          counts_i <- as.numeric(counts_r[i, ])
          
          counts_df <- .prepare_inputs(
            counts     = counts_i,
            pseudotime = pseudotime,
            pcw        = pcw,
            libsize    = libsize
          )
          
          ans <- try(fit_cycspline(counts_df, k = k, knots = 10, method = method),
                     silent = TRUE)
          
          if (!inherits(ans, "try-error")) {
            yh_loc[jj, ] <- ans$yhat
            gf_loc[jj, ] <- c(ans$p_value, ans$avg_log2fc, ans$peak_phase)
          }
        }
        p()
        list(rows = ii, yhat = yh_loc, fit = gf_loc)
      },
      future.seed = TRUE
    )
  })
  
  # Stitch results back
  for (o in outs) {
    yhats[o$rows, ] <- o$yhat
    gene_fit$p_value[o$rows]     <- o$fit[, 1]
    gene_fit$avg_log2fc[o$rows] <- o$fit[, 2]
    gene_fit$peak_phase[o$rows]  <- o$fit[, 3]
  }
  gene_fit$p_adjust <- p.adjust(gene_fit$p_value, method = "fdr")
  
  list(yhats = yhats, gene_fit = gene_fit)
}

