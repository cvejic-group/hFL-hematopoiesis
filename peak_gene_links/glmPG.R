
suppressPackageStartupMessages({
  #.libPaths("~/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
  library(tidyverse)
  library(lme4)
})
options(stringsAsFactors = F)

# for loop
lmmPG <- function(object=NULL,
                   inter_pcw = FALSE,
                   correct_batch=TRUE # include Sex, libraryID, donorID, batch in the model. TRUE for HSC, FALSE for others
                   ) {
  # info
  peak_info <- object$links
  num_link <- nrow(peak_info)
  meta_data <- object$mdata
  gm <- object$rna
  pm <- object$atac

  # helper: choose fixed/random effects formula
  f_base <- f_full <- NULL
  if (inter_pcw) {
    if (correct_batch) {
      rhs_base <- "atac + PCWsca + m_mt + Sex + (1|donorID) + (1|libraryID) + (1|Batch)"
      rhs_full <- paste0(rhs_base, " + atac:PCWsca")
    } else {
      rhs_base <- "atac + PCWsca + m_mt"
      rhs_full <- paste0(rhs_base, " + atac:PCWsca")
    }
  } else {
    if (correct_batch) {
      rhs_base <- "m_mt + Sex + (1|donorID) + (1|libraryID) + (1|Batch)"
      rhs_full <- paste0("atac + ", rhs_base)
    } else {
      rhs_base <- "m_mt"
      rhs_full <- paste0("atac + ", rhs_base)
    }
  }
  f_base <- as.formula(paste("exp ~", rhs_base))
  f_full <- as.formula(paste("exp ~", rhs_full))
  # choose fitting function
  fit_fun <- if (correct_batch) lmer else lm

  # loop
  res_lst <- lapply(seq_len(num_link), function(n){
    this_gene <- peak_info[n, 1] # GENE is FIRST COLUMN OF PEAK.INFO
    this_peak <- peak_info[n, 2] # PEAK is SECOND COLUMN OF PEAK.INFO
    cat(paste0(this_gene, ":", this_peak, ", ", n, "/", num_link), "\n")
    # prep
    df2 <- meta_data
    df2$exp <- as.numeric(gm[this_gene, , drop = TRUE])
    df2$atac <- as.numeric(pm[this_peak, , drop = TRUE])

    # fit models
    fit_base <- fit_fun(f_base, data = df2)
    fit_full <- fit_fun(f_full, data = df2)

    # model comparison
    comp <- anova(fit_base, fit_full)
    aic_table <- AIC(fit_base, fit_full)

    # extract model comparison stats
    if (correct_batch) {
      loglik_base <- comp$logLik[1]
      loglik_full <- comp$logLik[2]
      chisq       <- comp$Chisq[2]
      chidf       <- comp$Df[2]
      p_value     <- comp$`Pr(>Chisq)`[2]
    } else {
      loglik_base <- NA
      loglik_full <- NA
      chisq       <- NA
      chidf       <- comp$Df[2]
      p_value     <- comp$`Pr(>F)`[2]
    }

    # summary on full model
    sum_full <- if (correct_batch) {
      fit_full_reml <- update(fit_full, REML = TRUE)
      summary(fit_full_reml)
    } else {
      summary(fit_full)
    }
    coefs_atac <- sum_full$coefficients["atac", ]

    # create output dataframe
    if (inter_pcw) {
      beta_pcw <- sum_full$coefficients["PCWsca", ]
      beta_inter <- sum_full$coefficients["atac:PCWsca", ]
      data.frame(
        row.names = NULL,
        gene = this_gene,
        peak = this_peak,
        beta_atac = coefs_atac["Estimate"],
        se_atac = coefs_atac["Std. Error"],
        t_atac = coefs_atac["t value"],
        beta_pcw = beta_pcw["Estimate"],
        se_pcw = beta_pcw["Std. Error"],
        t_pcw = beta_pcw["t value"],
        beta_inter = beta_inter["Estimate"],
        se_inter = beta_inter["Std. Error"],
        t_inter = beta_inter["t value"],
        aic_base = aic_table$AIC[1],
        aic_full = aic_table$AIC[2],
        loglik_base = loglik_base,
        loglik_full = loglik_full,
        chisq = chisq,
        chidf = chidf,
        p_value = p_value
      )
    } else {
      data.frame(
        row.names = NULL,
        gene = this_gene,
        peak = this_peak,
        beta_atac = coefs_atac["Estimate"],
        se_atac = coefs_atac["Std. Error"],
        t_atac = coefs_atac["t value"],
        aic_base = aic_table$AIC[1],
        aic_full = aic_table$AIC[2],
        loglik_base = loglik_base,
        loglik_full = loglik_full,
        chisq = chisq,
        chidf = chidf,
        p_value = p_value
      )
    }
  })
  res <- do.call(rbind, res_lst)
  return(res)
}

# Test code for lmmPG function
test_lmmPG <- function() {
  mdata <- qs::qread("/work/home/project/20231127_DevM/devm_r432/output/HSC.pg_prep/mdata.qs")
  gm <- qs::qread("/work/home/project/20231127_DevM/devm_r432/output/HSC.pg_prep/gm.qs")
  pm <- qs::qread("/work/home/project/20231127_DevM/devm_r432/output/HSC.pg_prep/pm.qs")
  link.list <- qs::qread("/work/home/project/20231127_DevM/devm_r432/output/HSC.pg_prep/cis.g2p_list.qs")

  # Get the corresponding data.frame from the list:
  links <- as.data.frame(link.list[[1]])[1:10,]
  # object_lst
  obj_lst <- list(
    mdata = mdata,
    rna = gm,
    atac = pm,
    links = links
  )
  res_df1 <- lmmPG(object=obj_lst, inter_pcw=FALSE)
  res_df2 <- lmmPG(object=obj_lst, inter_pcw=TRUE)
}

# Test code for lmmPG function
test_lmmPG2 <- function() {
  mdata <- qs::qread("/work/home/project/20231127_DevM/E2G/lmmPG/GP/mc_mdata.qs")
  gm <- qs::qread("/work/home/project/20231127_DevM/E2G/lmmPG/GP/mc_gm_scale.qs")
  pm <- qs::qread("/work/home/project/20231127_DevM/E2G/lmmPG/GP/mc_pm_scale.qs")
  link.list <- qs::qread("/work/home/project/20231127_DevM/E2G/lmmPG/GP/cis.g2p_list.qs")

  # Get the corresponding data.frame from the list:
  links <- as.data.frame(link.list[[1]])[1:10,]
  # object_lst
  obj_lst <- list(
    mdata = mdata,
    rna = gm,
    atac = pm,
    links = links
  )
  res_df1 <- lmmPG(object=obj_lst, inter_pcw=FALSE, correct_batch = FALSE)
}
