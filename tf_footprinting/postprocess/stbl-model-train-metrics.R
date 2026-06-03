library(tidyverse)
library(jsonlite)
library(magrittr)
library(writexl)

OUT_XLSX="analysis/tables/chrombpnet-model-training-metrics.xlsx"

## Load celltype metrics
FILES_JSON_PATTERN_CT_CBPN="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/*/fold_*/evaluation/*._chrombpnet_metrics.json"
FILES_JSON_PATTERN_CT_BIAS="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/*/fold_*/evaluation/bias_metrics.json"
COLS_CELLTYPES=c(
  "HSC"="#E41A1C",
  "GP"="#E0FFFF",  "Granulocyte"="#B3CDE3",
  "MEMP-t"="#E6AB02", "MEMP"="#FF7F00", "MEP"="#CD661D", "MEMP-Mast-Ery"="#FDCDAC",
  "MEMP-Ery"="#E9967A", "Early-Ery"="#CD5555", "Late-Ery"="#8B0000",
  "MEMP-MK"="#663C1F", "MK"="#40E0D0",
  "MastP-t"="#1E90FF", "MastP"="#1F78B4", "Mast"="#253494",
  "MDP"="#E6F5C9", "Monocyte"="#005A32", "Kupffer"="#00EE00",
  "cDC1"="#ADFF2F", "cDC2"="#B3DE69", "pDC"="#4DAF4A", "ASDC"="#CDC673",
  "LMPP"="#FFF2AE", "LP"="#FFD92F", "Cycling-LP"="#FFFF33",
  "PreProB"="#FFF0F5", "ProB-1"="#FFB5C5", "ProB-2"="#E78AC3",
  "Large-PreB"="#CD1076", "Small-PreB"="#FF3E96", "IM-B"="#FF00FF",
  "NK"="#A020F0", "ILCP"="#49006A", "T"="#984EA3",
  "Hepatocyte"="#666666", "Endothelia"="#333333"
)

df_metrics_ct_cbpn <- tibble(file=Sys.glob(FILES_JSON_PATTERN_CT_CBPN)) %>%
  mutate(sample=str_match(file, "([^/]+)/fold")[,2]) %>%
  filter(sample %in% names(COLS_CELLTYPES)) %>%
  mutate(sample=fct(sample, levels=names(COLS_CELLTYPES))) %>%
  mutate(fold=str_match(file, "fold_([0-9]+)/")[,2]) %>%
  mutate(raw_metrics=map(file, fromJSON)) %>%
  mutate(spearmanr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$spearmanr),
         pearsonr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$pearsonr),
         mse=sapply(raw_metrics, function (l) l$counts_metrics$peaks$mse),
         median_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_jsd),
         median_norm_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_norm_jsd)) %>%
  mutate(raw_metrics=NULL, file=NULL)

df_metrics_ct_bias <- tibble(file=Sys.glob(FILES_JSON_PATTERN_CT_BIAS)) %>%
  mutate(sample=str_match(file, "([^/]+)/fold")[,2]) %>%
  filter(sample %in% names(COLS_CELLTYPES)) %>%
  mutate(sample=fct(sample, levels=names(COLS_CELLTYPES))) %>%
  mutate(fold=str_match(file, "fold_([0-9]+)/")[,2]) %>%
  mutate(raw_metrics=map(file, fromJSON)) %>%
  mutate(bias_spearmanr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$spearmanr),
         bias_pearsonr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$pearsonr),
         bias_mse=sapply(raw_metrics, function (l) l$counts_metrics$peaks$mse),
         bias_median_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_jsd),
         bias_median_norm_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_norm_jsd)) %>%
  mutate(raw_metrics=NULL, file=NULL)

df_metrics_ct <- full_join(df_metrics_ct_cbpn, df_metrics_ct_bias, by=c("sample", "fold")) %>%
  arrange(sample, fold)

## Load HSC metrics
FILES_JSON_PATTERN_HSC_CBPN="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/HSC_PCW*/fold_*/evaluation/*._chrombpnet_metrics.json"
FILES_JSON_PATTERN_HSC_BIAS="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/HSC_PCW*/fold_*/evaluation/bias_metrics.json"
COLS_PCW=c(
  "HSC_PCW5"="#FFBF00",
  "HSC_PCW6"="#ff986d",
  "HSC_PCW7"="#ff8782",
  "HSC_PCW8"="#ff7b9b",
  "HSC_PCW9"="#f176b4",
  "HSC_PCW10"="#d97fce",
  "HSC_PCW11"="#b989e2",
  "HSC_PCW12"="#9094ed",
  "HSC_PCW13"="#60a4f2",
  "HSC_PCW14"="#32b1eb",
  "HSC_PCW15"="#29badb",
  "HSC_PCW16"="#4bc0c8",
  "HSC_PCW17"="#8dd7db",
  "HSC_PCW18"="#c0dedf")

df_metrics_hsc_cbpn <- tibble(file=Sys.glob(FILES_JSON_PATTERN_HSC_CBPN)) %>%
  mutate(sample=str_match(file, "([^/]+)/fold")[,2]) %>%
  filter(sample %in% names(COLS_PCW)) %>%
  mutate(sample=fct(sample, levels=names(COLS_PCW))) %>%
  mutate(fold=str_match(file, "fold_([0-9]+)/")[,2]) %>%
  mutate(raw_metrics=map(file, fromJSON)) %>%
  mutate(spearmanr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$spearmanr),
         pearsonr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$pearsonr),
         mse=sapply(raw_metrics, function (l) l$counts_metrics$peaks$mse),
         median_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_jsd),
         median_norm_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_norm_jsd)) %>%
  mutate(raw_metrics=NULL, file=NULL)

df_metrics_hsc_bias <- tibble(file=Sys.glob(FILES_JSON_PATTERN_HSC_BIAS)) %>%
  mutate(sample=str_match(file, "([^/]+)/fold")[,2]) %>%
  filter(sample %in% names(COLS_PCW)) %>%
  mutate(sample=fct(sample, levels=names(COLS_PCW))) %>%
  mutate(fold=str_match(file, "fold_([0-9]+)/")[,2]) %>%
  mutate(raw_metrics=map(file, fromJSON)) %>%
  mutate(bias_spearmanr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$spearmanr),
         bias_pearsonr=sapply(raw_metrics, function (l) l$counts_metrics$peaks$pearsonr),
         bias_mse=sapply(raw_metrics, function (l) l$counts_metrics$peaks$mse),
         bias_median_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_jsd),
         bias_median_norm_jsd=sapply(raw_metrics, function (l) l$profile_metrics$peaks$median_norm_jsd)) %>%
  mutate(raw_metrics=NULL, file=NULL)

df_metrics_hsc <- full_join(df_metrics_hsc_cbpn, df_metrics_hsc_bias, by=c("sample", "fold")) %>%
  arrange(sample, fold)

## append
df_metrics_all <- bind_rows(df_metrics_ct, df_metrics_hsc)

write_xlsx(df_metrics_all, path=OUT_XLSX)
