#!/usr/bin/env Rscript

# cd ~/project/20231127_DevM/devm_r432
# nohup Rscript code/05.get_rna_pb_by_donor.R

DOCNAME <- "RNA.pb_by_donor"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
# srat
library(Seurat)

# load
srat <- readRDS("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_rna_clustered.v02.rds")
srat

# pb by donor
srat_pb <- AggregateExpression(srat, group.by = "donorID", return.seurat = TRUE)
srat_pb

# save
saveRDS(srat_pb, file = here::here("output", DOCNAME, "srat_pb_by_donor.rds"))

sessionInfo()
Sys.Date()

