#!/usr/bin/env Rscript

DOCNAME <- "fracOpen"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)

# srat
library(Seurat)
library(ArchR)
# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 12)

# model
library(lme4)
library(glmmTMB)

# Function to detect functional equivalence
is_functionally_equivalent <- function(col1, col2) {
  # Check if each unique value in col1 maps to exactly one unique value in col2
  map1 <- all(tapply(df[[col2]], df[[col1]], function(x) length(unique(x)) == 1))

  # Check the reverse mapping
  map2 <- all(tapply(df[[col1]], df[[col2]], function(x) length(unique(x)) == 1))

  return(map1 && map2)  # Both directions must hold for full equivalence
}

# Function to detect functional equivalence
not_equivalent <- function(col1=NULL, col2=NULL) {
  # Check if each unique value in col1 maps to exactly one unique value in col2
  map1 <- all(tapply(col2, col1, function(x) length(unique(x)) == 1))
  # Check the reverse mapping
  map2 <- all(tapply(col1, col2, function(x) length(unique(x)) == 1))
  return(!(map1 && map2))  # Both directions must hold for full equivalence
}


data_dir <- "archr_48FL.eachCell"
cells <- c("HSC", "MEMP-t", "MastP-t", "MDP", "LMPP", "LP", "Cycling-LP",
           "GP", "Granulocyte",
           "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP", "Mast",
           "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T")

for (cell in cells) {
  print(cell)
  # load
  proj_dir <- here::here('output', data_dir, cell)
  proj <- loadArchRProject(path = proj_dir)
  # get pm
  proj_se <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "PeakMatrix"
  )
  pm <- assays(proj_se)$PeakMatrix
  # fracOpen
  f_pm <- colSums(pm > 0)/nrow(pm)
  # get tm
  proj_se <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    binarize = TRUE
  )
  tm <- assays(proj_se)$TileMatrix
  # fracOpen
  f_tm <- colSums(tm > 0)/nrow(tm)
  # build df
  df <- as.data.frame(proj@cellColData) %>%
    rownames_to_column("barcode") %>%
    mutate(logFrag = log(nFrags)) %>%
    dplyr::select(barcode, nFrags, logFrag, libraryID, donorID, PCW, Sex, Batch) %>%
    left_join(data.frame(barcode = names(f_pm), fracOpenPm = f_pm), by = "barcode") %>%
    left_join(data.frame(barcode = names(f_tm), fracOpenTm = f_tm), by = "barcode") %>%
    mutate(libraryID = factor(libraryID),
           donorID = factor(donorID),
           Sex = factor(Sex),
           Batch = factor(Batch))
  # Apply pairwise
  cat_vars <- names(df[, c("libraryID", "donorID", "Sex", "Batch")])
  comb <- combn(cat_vars, 2)
  # Identify functionally equivalent pairs
  apply(comb, 2, function(pair) {
    if (is_functionally_equivalent(pair[1], pair[2])) {
      cat(sprintf("Columns '%s' and '%s' are functionally equivalent (one-to-one mapping).\n",
                  pair[1], pair[2]))
    }
  })
  # residuls (PM)
  f <- "fracOpenPm ~ logFrag + (1|libraryID) + (1|donorID) + (1|Batch) + (1|Sex)"
  f <- as.formula(f)
  model_lmm <- lmer(f,  data = df, control = lmerControl(optimizer = "bobyqa"))
  model_beta <- glmmTMB(f, family = beta_family(), data = df)
  df$fracOpenPm_lmer <- residuals(model_lmm)
  df$fracOpenPm_beta <- residuals(model_beta)
  # residuls (TM)
  f1 <- "fracOpenTm ~ logFrag + (1|libraryID) + (1|donorID) + (1|Batch) + (1|Sex)"
  f1 <- as.formula(f1)
  model_lmm1 <- lmer(f1,  data = df, control = lmerControl(optimizer = "bobyqa"))
  model_beta1 <- glmmTMB(f1, family = beta_family(), data = df)
  df$fracOpenTm_lmer <- residuals(model_lmm1)
  df$fracOpenTm_beta <- residuals(model_beta1)
  # save
  saveRDS(df, file = here::here("output", DOCNAME, paste0(cell, ".fracOpen.rds")))
}

# load some cellmeta
cell_order <- c("HSC", "GP", "Granulocyte",
                "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
                "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
                "NK", "ILCP", "T",
                "Hepatocyte", "Endothelia")
df_meta = read.csv('/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_cellmeta.v01.csv') |>
  dplyr::select(barcode = X, anno_wnn_v51, Batch, S.Score = rna.S.Score, G2M.Score = rna.G2M.Score, Phase = rna.Phase) %>%
  mutate(barcode = str_replace(barcode, "_", "#"),
         anno_wnn_v51 = factor(anno_wnn_v51, levels = cell_order),
         Batch = factor(Batch)) %>%
  droplevels()


for (pcw in 5:18) {
  print(pcw)
  # load
  proj_dir <- file.path("/work/DevM_analysis/data/ArchR_perPCW", paste0("archr_48FL.nonBinarized.PCW", pcw))
  proj <- loadArchRProject(path = proj_dir)
  # get tm
  proj_se <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "TileMatrix",
    binarize = TRUE
  )
  tm <- assays(proj_se)$TileMatrix
  # fracOpen
  f_tm <- colSums(tm > 0)/nrow(tm)
  # build df
  df <- as.data.frame(proj@cellColData) %>%
    rownames_to_column("barcode") %>%
    mutate(logFrag = log(nFrags)) %>%
    left_join(df_meta, by = "barcode") %>%
    dplyr::select(barcode, nFrags, logFrag, libraryID, donorID, PCW, Sex, Batch, anno_wnn_v51, S.Score, G2M.Score, Phase) %>%
    left_join(data.frame(barcode = names(f_tm), fracOpenTm = f_tm), by = "barcode") %>%
    mutate(libraryID = factor(libraryID),
           donorID = factor(donorID),
           Sex = factor(Sex)) %>%
    droplevels()
  # covariates
  covariates <- c("libraryID")
  # donorID
  if (length(unique(df$donorID)) != 1) {
    if (not_equivalent(df$libraryID, df$donorID)) {
      covariates <- c(covariates, "donorID")
    }
  }
  # Sex
  if (length(unique(df$Sex)) != 1 &&
      not_equivalent(df$libraryID, df$Sex) &&
      (!("donorID" %in% covariates) || not_equivalent(df$donorID, df$Sex))) {
    covariates <- c(covariates, "Sex")
  }
  # Batch
  if (length(unique(df$Batch)) != 1 &&
      not_equivalent(df$libraryID, df$Batch) &&
      (!("donorID" %in% covariates) || not_equivalent(df$donorID, df$Batch)) &&
      (!("Sex" %in% covariates) || not_equivalent(df$Sex, df$Batch))) {
    covariates <- c(covariates, "Batch")
  }
  # residuls (TM)
  f1 <- paste0("fracOpenTm ~ logFrag + ", paste0("(1|", covariates, ")", collapse = "+"))
  f1 <- as.formula(f1)
  model_lmm1 <- lmer(f1,  data = df, control = lmerControl(optimizer = "bobyqa"))
  model_beta1 <- glmmTMB(f1, family = beta_family(), data = df)
  df$fracOpenTm_lmer <- residuals(model_lmm1)
  df$fracOpenTm_beta <- residuals(model_beta1)
  # save
  saveRDS(df, file = here::here("output", DOCNAME, paste0("PCW", pcw, ".fracOpen.rds")))
}




sessionInfo()
Sys.Date()
