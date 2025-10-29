#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

# load script
#.libPaths("~/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
source("~/project/20231127_DevM/E2G/lmmPG/lmmPG.R")

# INPUTS
option_list = list(
  make_option(c("-n", "--node"), type="integer", default=NULL,
              help="JOB ARRAY number"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="object list"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output folder"),
  make_option(c("--time"), action="store_true", default=FALSE,
              help="Interaction with time"),
  make_option(c("--nocorrect"), action="store_true", default=FALSE,
              help="Correct covariates (Sex, library, donor, batch)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#print(opt)

# Load:
obj_lst = qs::qread(opt$input)

# Get the corresponding data.frame from the list:
obj_lst$links <- as.data.frame(obj_lst$link_lst[[opt$node]])

# Run
if (opt$nocorrect) {
  res_df <- lmmPG(object=obj_lst, inter_pcw=opt$time, correct_batch=FALSE)
} else {
  res_df <- lmmPG(object=obj_lst, inter_pcw=opt$time, correct_batch=TRUE)
}

# out per block
filename <- paste0(opt$output, "/lmmPG_", opt$node, ".tsv")
write.table(res_df, file = filename, row.names = F, col.names = T, quote = F)

