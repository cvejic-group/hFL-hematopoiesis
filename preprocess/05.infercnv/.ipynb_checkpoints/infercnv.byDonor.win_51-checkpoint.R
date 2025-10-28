#!/usr/bin/env Rscript

# usage
# nohup Rscript analysis/infercnv.byDonor.win_51.R donor > analysis/infercnv.byDonor.win_51.log &

# nested expression evaluated
options(expressions=10000)
# traceback
options(error = function() traceback(2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(infercnv))

# args
args <- commandArgs(trailingOnly=TRUE)

# test args
if (length(args) == 0) {
  stop("At least one argument must be supplied (patient name).", call.=FALSE)
} else {
  patient_name <- args[1]
}

message(patient_name)

# test
#patient_name = "FL1"
#patient_name = "FL2"

# out dir
DOCNAME <- "infercnv.byDonor.win_51"
out_dir <- here::here('output', DOCNAME, patient_name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# gene order file
genes_pos <- here::here("data", "genes_pos.txt")

# the ref cells
srat_ref <- readRDS(here::here("data", "mergedHDBRCtrl.4_infercnv.rds"))
ref_group_names <- unique(srat_ref$libraryID)

# case samples
seu_dir <- "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj"
case_files <- grep("rna", grep(patient_name, list.files(seu_dir), value = TRUE), value = TRUE)
srat_case_lst <- lapply(case_files, function(f){
  sample_name <- str_split(f, '_', simplify = T)[1]
  x <- readRDS(file.path(seu_dir, f))
  x[["RNA"]] <- as(object = x[["RNA"]], Class = "Assay")
  x <- RenameCells(x, new.names = paste0(sample_name, "_", colnames(x)))
  x
})
srat_case <- purrr::reduce(srat_case_lst, merge)
srat_case

# remove the sample from ref_group_names if it is a control
`%ni%` <- Negate(`%in%`)
sample_names <- str_split(case_files, "_", simplify = T)[,1]
ref_group_names <- ref_group_names[ref_group_names %ni% c(sample_names)]

# group ref by donor
ref_group_names <- str_split(ref_group_names, "-", simplify = T)[,1] |> unique()

# remove cells from srat_ref if it is a control
srat_ref <- subset(srat_ref, subset = libraryID %in% c(sample_names), invert = TRUE)

# merge case/ctrl
srat <- merge(srat_case, srat_ref)

# exp mat
expMat <- GetAssayData(srat, assay = 'RNA', layer = 'counts')

# cell group info
groupFile <- file.path(out_dir, 'groupFile.txt')
tmp_ref <- srat_ref@meta.data |>
  rownames_to_column("barcode") |>
  select(barcode, Sample = donorID)
# pooled libs
if (patient_name %in% c("FLpool1", "FLpool2", "FLpool3", "FLpool4", "FLpool5", "FLpool6")) {
  tmp_case <- srat_case@meta.data |>
    rownames_to_column("barcode") |>
    select(barcode, Sample = donorID)
} else {
  tmp_case <- srat_case@meta.data |>
    rownames_to_column("barcode") |>
    select(barcode, Sample = libraryID)
}
groupInfo <- rbind(tmp_ref, tmp_case)
write.table(groupInfo, file = groupFile, sep = '\t',
            quote = F, col.names = F, row.names = F)

# prepare infercnv obj
z <- CreateInfercnvObject(raw_counts_matrix = expMat,
                          annotations_file = groupFile,
                          gene_order_file = genes_pos,
                          ref_group_names = ref_group_names,
                          chr_exclude = c("chrY", "chrM"),
                          delim = "\t")
# run
Cstack_info()
z <- infercnv::run(z,
                   cutoff = 0.1,
                   window_length = 51,
                   cluster_by_groups = TRUE,
                   denoise = TRUE,
                   num_threads = 16,
                   output_format = 'pdf',
                   out_dir = out_dir)

