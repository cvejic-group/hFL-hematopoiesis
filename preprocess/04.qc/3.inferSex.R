###################################
#                                 #
#         02.2 Infer sex          #
#                                 #
###################################


## Project: Developmental multiome
## Date: 08.2024
## Author: Pernille
## Script name: 02.2.Infer-sex.R

## Script description: Infer the sex of each sample based on sex-specific gene expression

.libPaths("/work/home/Software/R")

library(dplyr)
library(tidyr)
library(ggplot2)

source("Developmental_multiome/Analysis/Utils/GeneralUtils.R")
source("/work/home/Developmental_multiome/Analysis/Utils/ListMarkers.R")

##---------------------------------------------##
##---------------0. Loading data---------------##
##---------------------------------------------##

# from https://www.sciencedirect.com/science/article/pii/S009286742300908X
M_F_markers <- c("XIST", # Female
                 "KDM5D",  "ZFY", "PRKY", "USP9Y", "TTTY14", "UTY", "RPS4Y1", "DDX3Y", "EIF1AY") # Male
                #"TTTY15", "PSMA6P1", "GYG2P1", "TXLNGY", very low in all samples

# seurat path
my_data <- list.files(path = "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj",
                      full.names = TRUE, recursive = F)
my_data <- my_data[grep("rna",my_data)]
my_names <- sapply(strsplit(sapply(strsplit(my_data, "/"), '[',8), "_"), '[',1)
names(my_data) <- my_names

#For CB:
#my_data <- my_data[grep("CB", my_data)]

##---------------------------------------------##
##-------1. Extract sex genes expression-------##
##---------------------------------------------##

full_exp_mat <- setNames(data.frame(matrix(ncol = length(M_F_markers) + 2,
                                           nrow = 0)),
                         c("donorID", "sampleID", M_F_markers))
for (d in names(my_data)){

  cat(d, "starts", "\n")
  rna <- readRDS(my_data[[d]])

  exp_mat <- t(as.matrix(rna[["RNA"]]$data[rownames(rna) %in% M_F_markers,]))
  add <- matrix(0, nrow = ncol(rna), ncol = sum(!M_F_markers %in% rownames(rna)))
  rownames(add) <- rownames(exp_mat); colnames(add) <- M_F_markers[!M_F_markers %in% rownames(rna)]
  exp_mat <- cbind(exp_mat, add)
  exp_mat <- bind_cols(rna@meta.data[, c("donorID", "sampleID"), drop = F], as.data.frame(exp_mat))
  rownames(exp_mat) <- paste0(exp_mat$donorID, "_", rownames(exp_mat))
  full_exp_mat <- rbind(full_exp_mat, exp_mat[, colnames(full_exp_mat)])
}

##---------------------------------------------##
##-----------------2. Dot plot-----------------##
##---------------------------------------------##
#Same done for sampleID
data <- pivot_longer(full_exp_mat %>% select(!sampleID),
                     -donorID, names_to="Gene", values_to="Expression")
data_summary <- data %>%
  group_by(donorID, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)
data_summary$Gene <- factor(data_summary$Gene, levels = M_F_markers)

dot_plot <- ggplot(data_summary, aes(x=Gene, y=donorID)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% detected", range = c(0,6)) +
  scale_fill_gradientn(colours = rev(viridisLite::plasma(100)),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Average\nexpression") +
  ylab("Cluster") + xlab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=14))
dot_plot

ggsave(plot = dot_plot, filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Dotplot_sexGenes_perdonorID.pdf",
       height = 12, width = 6)
#ggsave(plot = dot_plot, filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/CB_Dotplot_sexGenes_perdonorID.pdf",
#       height = 7, width = 6)

ggplot(full_exp_mat %>% filter(donorID == "FL6PCW2d0") %>%
         #rowwise() %>%
         mutate(meanY = rowMeans(select(., M_F_markers[-1]))),
       aes(XIST, meanY)) + geom_point()

##---------------------------------------------##
##-----------------2. Save sex-----------------##
##---------------------------------------------##

Sex <- as.data.frame(Sex) %>% tibble::rownames_to_column("sampleID")
write.table(Sex, "/work/DevM_analysis/utils/sex_info.DevM.txt", sep = "\t")

sam_info <- read.table("/work/DevM_analysis/utils/sam_info.DevM.20.08.24.txt", sep = "\t", header = T)
sam_info <- left_join(sam_info, Sex)

write.table(sam_info, "/work/DevM_analysis/utils/sam_info.DevM.20.08.24.txt", sep = "\t")

