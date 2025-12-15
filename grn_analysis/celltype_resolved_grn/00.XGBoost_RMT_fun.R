###############################
### XGBoost helper function ###
###############################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

setwd("~/work/")
source("./00.Initialization.R")

# Dependencies
library(xgboost)
library(SHAPforxgboost)
library(pROC)
library(future)
library(furrr)
library(progressr)
library(matrixStats)


train_OVR_models_targetCT <- function(
    feature_matrix_all,
    cell_types,
    model_type = "xgboost",
    nfold = 5,
    nrounds = 500,
    max_depth = 3,
    shap_ratio_cutoff = 2,
    shap_ratio_diff_cutoff = -500,
    seed = 123,
    do_zscore = TRUE,
    do_MinMax = FALSE,
    balance_classes = TRUE,
    n_iter = 30,
    parallel = TRUE,
    n_workers = 5,
    nthread_xgb = 15
) {
  
  # Parallel & progress setup
  plan(if (parallel) multisession else sequential, workers = n_workers)
  handlers("txtprogressbar")
  
  all_results <- list()
  unique_cts <- unique(cell_types)
  
  # Preprocess full matrix
  feature_matrix_all.Z <- feature_matrix_all
  if (do_zscore) feature_matrix_all.Z <- scale(feature_matrix_all)
  if (do_MinMax) feature_matrix_all.Z <- t(apply(feature_matrix_all.Z, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-6)))
  
  # Main loop per cell type
  for (ct in unique_cts) {
    message("Processing cell type: ", ct)
    featZ    <- feature_matrix_all.Z
    feat_raw <- feature_matrix_all
    labels   <- ifelse(cell_types == ct, 1, 0)
    features <- colnames(featZ)
    p <- length(features)
    
    # Initialize accumulators
    real_pos_sum     <- numeric(p)
    real_neg_sum     <- numeric(p)
    null_pos_sum     <- numeric(p)
    null_neg_sum     <- numeric(p)
    ratio_sum        <- numeric(p)
    ratio_neg_sum    <- numeric(p)
    imp_gain_sum     <- numeric(p)
    imp_cover_sum    <- numeric(p)
    imp_freq_sum     <- numeric(p)
    auc_sum          <- 0
    logFC_mat        <- matrix(NA, nrow = p, ncol = n_iter, dimnames = list(features, NULL))
    
    # Iterative subsampling with progress
    with_progress({
      p_fun <- progressr::progressor(along = seq_len(n_iter))
      iter_res <- future_map(seq_len(n_iter), function(i) {
        p_fun()
        set.seed(seed + i)
        pos_idx <- which(labels == 1)
        neg_idx <- which(labels == 0)
        if (balance_classes) {
          neg_sample <- sample(neg_idx, length(pos_idx), replace = FALSE)
          idx <- c(pos_idx, neg_sample)
        } else {
          idx <- seq_along(labels)
        }
        X_sub <- featZ[idx, , drop = FALSE]
        y_sub <- labels[idx]
        raw_sub <- feat_raw[idx, , drop = FALSE]
        
        # CV AUC
        dcv <- xgb.DMatrix(data = X_sub, label = y_sub)
        params <- list(objective = "binary:logistic", eval_metric = "auc",
                       max_depth = max_depth, 
                       eta = 0.1, 
                       nthread = nthread_xgb)
        cv <- xgb.cv(params, dcv, nfold = nfold, nrounds = nrounds,
                     verbose = 0, prediction = TRUE)
        auc_i <- pROC::auc(y_sub, cv$pred)
        
        # Train final model
        model <- xgb.train(params, dcv, nrounds = nrounds, verbose = 0)
        
        # SHAP real
        shap_r <- shap.values(xgb_model = model, X_train = X_sub)$shap_score
        rp <- colMeans(shap_r[y_sub == 1, , drop = FALSE])
        rn <- colMeans(shap_r[y_sub == 0, , drop = FALSE])
        
        # SHAP null
        null_mat <- matrix(rnorm(nrow(X_sub) * p), nrow = nrow(X_sub))
        colnames(null_mat) <- features
        dnull <- xgb.DMatrix(data = null_mat, label = y_sub)
        mnull <- xgb.train(params, dnull, nrounds = nrounds, verbose = 0)
        shap_n <- shap.values(xgb_model = mnull, X_train = null_mat)$shap_score
        np <- colMeans(shap_n[y_sub == 1, , drop = FALSE])
        nn <- colMeans(shap_n[y_sub == 0, , drop = FALSE])
        
        # SHAP ratio per iter
        ratio_i     <- abs(rp) / (abs(np) + 1e-6)
        ratio_neg_i <- abs(rn) / (abs(nn) + 1e-6)
        
        # Feature importance
        imp <- xgb.importance(model = model)
        gain_vec  <- setNames(imp$Gain, imp$Feature)[features];  gain_vec[is.na(gain_vec)]   <- 0
        cover_vec <- setNames(imp$Cover, imp$Feature)[features]; cover_vec[is.na(cover_vec)] <- 0
        freq_vec  <- setNames(imp$Frequency, imp$Feature)[features]; freq_vec[is.na(freq_vec)]  <- 0
        
        # logFC
        pos_expr <- colMeans(raw_sub[y_sub == 1, , drop = FALSE])
        neg_expr <- colMeans(raw_sub[y_sub == 0, , drop = FALSE])
        pos_expr[!is.finite(pos_expr)] <- 0
        neg_expr[!is.finite(neg_expr)] <- 0
        logFC_i <- log2(pos_expr + 1e-6) - log2(neg_expr + 1e-6)
        
        list(rp = rp, rn = rn, np = np, nn = nn,
             ratio = ratio_i, ratio_neg = ratio_neg_i,
             gain = gain_vec, cover = cover_vec, freq = freq_vec,
             auc = auc_i, logFC = logFC_i)
      }, .options = furrr_options(seed = TRUE))
    })
    
    # Aggregate results across iterations
    for (j in seq_len(n_iter)) {
      res <- iter_res[[j]]
      real_pos_sum      <- real_pos_sum      + res$rp
      real_neg_sum      <- real_neg_sum      + res$rn
      null_pos_sum      <- null_pos_sum      + res$np
      null_neg_sum      <- null_neg_sum      + res$nn
      ratio_sum         <- ratio_sum         + res$ratio
      ratio_neg_sum     <- ratio_neg_sum     + res$ratio_neg
      imp_gain_sum      <- imp_gain_sum      + res$gain
      imp_cover_sum     <- imp_cover_sum     + res$cover
      imp_freq_sum      <- imp_freq_sum      + res$freq
      auc_sum           <- auc_sum           + res$auc
      logFC_mat[, j]    <- res$logFC
    }
    
    # Compute means / medians
    shap_real_avg     <- real_pos_sum      / n_iter
    shap_real_neg_avg <- real_neg_sum      / n_iter
    shap_null_avg     <- null_pos_sum      / n_iter
    shap_null_neg_avg <- null_neg_sum      / n_iter
    shap_ratio_avg    <- ratio_sum         / n_iter
    shap_ratio_neg_avg<- ratio_neg_sum     / n_iter
    shap_ratio_diff   <- shap_ratio_avg - shap_ratio_neg_avg
    
    importance_df <- data.frame(
      Feature   = features,
      Gain      = imp_gain_sum  / n_iter,
      Cover     = imp_cover_sum / n_iter,
      Frequency = imp_freq_sum  / n_iter,
      cell_type = ct
    )
    auc_cv_mean <- auc_sum / n_iter
    logFC_med   <- rowMedians(logFC_mat)
    
    # Feature selection
    possible <- features[shap_ratio_avg > shap_ratio_cutoff & 
                           shap_ratio_diff > shap_ratio_diff_cutoff &
                           logFC_med > 0]
    confident<- possible[shap_real_avg[possible] > 0 & 
                           shap_real_neg_avg[possible] <= 0]
    
    # SHAP dataframe
    shap_df <- data.frame(
      Feature            = features,
      SHAP               = shap_real_avg,
      SHAP_neg           = shap_real_neg_avg,
      SHAP_random        = shap_null_avg,
      SHAP_random_neg    = shap_null_neg_avg,
      SHAP_ratio         = shap_ratio_avg,
      SHAP_ratio_neg     = shap_ratio_neg_avg,
      shap_ratio_diff    = shap_ratio_diff,
      logFC              = logFC_med,
      cell_type          = ct,
      CT_eReg_possible   = features %in% possible,
      CT_eReg_confident  = features %in% confident,
      selected           = features %in% confident
    )
    
    all_results[[ct]] <- list(
      AUC_CV_mean       = auc_cv_mean,
      importance        = importance_df,
      shap              = shap_df,
      possible_features = possible,
      confident_features= confident
    )
  }
  return(all_results)
}