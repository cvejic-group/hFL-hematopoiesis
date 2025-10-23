library(tidyverse)

sam_info <- read_csv("~/project/20231127_DevM/sam_info.csv")
samples <- sam_info$Sample

# collect outputs
df_doublet_lst <- lapply(samples, function(i){
  message(i)
  data_dir = "~/project/20231127_DevM/scDblFinder"
  fname <- file.path(data_dir, paste0(i, ".csv"))
  read_csv(fname) |>
    mutate(barcode = paste0(i, "_", barcode))
})
df_doublet <- do.call(rbind, df_doublet_lst)
write_csv(df_doublet, file = "scDblFinder_out_mergerd.csv")

# sample level summary
df_lst = lapply(1:length(samples), function(i) {
  data_dir = "~/project/20231127_DevM/scDblFinder"
  df = df_doublet_lst[[i]] %>%
    mutate(predicted_doublet = case_when(
      scDblFinder.class == "singlet" ~ FALSE,
      TRUE ~ TRUE
    ))
  n_cell = dim(df)[1]
  data.frame(
    Sample = samples[i],
    n_cell = n_cell,
    scDblFinder_doublet = sum(df$predicted_doublet),
    percent.doublet = sum(df$predicted_doublet)/n_cell * 100
  )
})

do.call(rbind, df_lst) %>%
  write_csv("scDblFinder_doublet_rate.csv")


