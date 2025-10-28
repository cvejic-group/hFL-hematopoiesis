###################################
#                                 #
#            Moscot plot          #
#             per PCW             #
#                                 #
###################################

## Project: Developmental multiome
## Date: 02.2025
## Author: Pernille
## Script name: 3. Summary_barplot

## Script description: Summary bar plots of moscot results

.libPaths("/work/home/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)

setwd("/work/DevM_analysis/03.trajectory/Moscot_Temporal/")

##---------------------------------------------##
##-----------------Prepare env-----------------##
##---------------------------------------------##

mashaGgplot2Theme <- list(
  theme_classic(base_size = 14) + 
    theme(text = element_text(size = 14)) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "white", size = 0.5,
                                          linetype = 2))
)
aes_turn_x_text <-  theme(axis.text.x = element_text(angle = 90, hjust = 1))

annot2col <- c("HSC" = "#E41A1C",
               "GP" = "#E0FFFF",
               "Granulocyte" = "#B3CDE3",
               "MEMP-t" = "#E6AB02",
               "MEMP" = "#FF7F00",
               "MEP" = "#CD661D",
               "MEMP-Mast-Ery" = "#FDCDAC",
               "MEMP-Ery" = "#E9967A",
               "Early-Ery" = "#CD5555",
               "Late-Ery" = "#8B0000",
               "MEMP-MK" = "#663C1F",
               "MK" = "#40E0D0",
               "MastP-t" = "#1E90FF",
               "MastP" = "#1F78B4",
               "Mast" = "#253494",
               "MDP" = "#E6F5C9",
               "Monocyte" = "#005A32",
               "Kupffer" = "#00EE00",
               "cDC1" = "#B3DE69",
               "cDC2" = "#ADFF2F",
               "pDC" = "#4DAF4A",
               "ASDC" = "#CDC673",
               "LMPP" = "#FFF2AE",
               "LP" = "#FFD92F",
               "Cycling-LP" = "#FFFF33",
               "PreProB" = "#FFF0F5",
               "ProB-1" = "#FFB5C5",
               "ProB-2" = "#E78AC3",
               "Large-PreB" = "#CD1076",
               "Small-PreB" = "#FF3E96",
               "IM-B" = "#FF00FF",
               "NK" = "#A020F0",
               "ILCP" = "#49006A",
               "T" = "#984EA3")

mapping2col <- c("#FFBF00","#ff986d", "#ff8782","#ff7b9b","#f176b4","#d97fce","#b989e2",
                 "#9094ed","#60a4f2","#32b1eb","#29badb","#4bc0c8","#8dd7db","#c0dedf")
mapping_lvl <- c(); for(w in 5:17){mapping_lvl <- c(mapping_lvl, paste0(w, "_to_", w+1))}
names(mapping2col) <- mapping_lvl

##---------------------------------------------##
##-------------------Function------------------##
##---------------------------------------------##
plot <- function(selection, data_sum = dataSum, data = moscot_res, 
                 min_value = 0, leg.pos = "none"){
  
  data_sum <- data_sum %>% filter(source %in% selection)
  data_sum <- data_sum %>% filter(M > min_value)
  data <- left_join(data_sum[, c("source", "target")], data)
  
  #dataSum$M_plus_S <- dataSum$M + dataSum$S
  #dataSum$M_min_S <- dataSum$M - dataSum$S
  
  p <- ggplot(data_sum, mapping = aes(x = target,  y = M, fill = target))  +  
    geom_col(position = "dodge", alpha = 0.9) +
    geom_errorbar(mapping = aes(ymin = pmax(M - S, 0), ymax = pmin(M + S, 1)), 
                  width = .2,
                  position = position_dodge(.9)) +
    geom_jitter(data, mapping = aes(x = target, y = value, col = mapping),
                size = 0.5, width = .3) + 
    facet_grid(~source, scales = "free_x", space = "free_x", switch = "x") + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank()) + 
    scale_fill_manual(values = annot2col) +
    scale_color_manual(values = mapping2col) +
    mashaGgplot2Theme + aes_turn_x_text +
    ylim(c(0,1)) + 
    ylab("Descendency probability") +
    xlab("") + theme(legend.position = leg.pos)
  
  return(p)
}

##---------------------------------------------##
##--------------------Main---------------------##
##---------------------------------------------##

for(direction in c("forward", "backward")){
  
  ## 1. Load data
  moscot_res <- list()
  for(w in 5:17){
    moscot_res[[paste0(w, "_to_", w+1)]] <- read.csv(paste0("data/Moscot_allFL_", w,
                                                            "_to_", w+1, "_", direction, ".csv"), 
                                                     header = T)
  }
  
  ## 2. Melt 
  for(w in names(moscot_res)){
    moscot_res[[w]] <- reshape2::melt(moscot_res[[w]], id.vars = "X")
    if(direction == "forward"){
      colnames(moscot_res[[w]]) <- c("source", "target", "value")
    }
    if(direction == "backward"){
      colnames(moscot_res[[w]]) <- c("target", "source", "value")
    }
    moscot_res[[w]]$mapping <- w
  }
  # Check if correctly assigned dependeing on the direction 
  #(forward rowsum == 1, backward Colsum == 1)
  testit::assert(round(sum((moscot_res[[1]] %>% filter(source == "HSC"))["value"]), 
                       0.0000000000001) == 1)
  
  ## 3. Unlist and prepare data
  moscot_res <- data.table::rbindlist(moscot_res)
  moscot_res$target <- as.character(moscot_res$target)
  moscot_res$target <- gsub("\\.", "-", moscot_res$target)
  moscot_res$source <- as.character(moscot_res$source)
  moscot_res$source <- gsub("\\.", "-", moscot_res$source)
  moscot_res$mapping <- factor(moscot_res$mapping,
                               levels = mapping_lvl)
  
  ## 4. Calculate mean and std
  moscot_res$source <- as.character(moscot_res$source)
  moscot_res$target <- as.character(moscot_res$target)
  
  dataSum <- moscot_res[, .(M = mean(value, na.rm = T), 
                            S = sd(value, na.rm = T)),
                        by = .(source, target)]
  ## 5. Factor
  dataSum$source <- factor(as.character(dataSum$source),
                           levels = names(annot2col))#, "unassigned"))
  dataSum$target <- factor(as.character(dataSum$target),
                           levels = names(annot2col))
  
  moscot_res$source <- factor(moscot_res$source, levels = names(annot2col))
  moscot_res$target <- factor(moscot_res$target, levels = names(annot2col))
  
  
  ## 6. Plots
  p1 <- plot(c("HSC", "GP", "MEMP-t", "MEMP", "MEP", "MEMP-Mast", "MEMP-Mast-Ery", 
               "MEMP-Ery", "Early-Ery", "Late-Ery", "MEMP-MK"),
             min_value = 0.02) 
  p2 <- plot(c("Granulocyte", "MK", "MastP-t", "MastP", "Mast", "MDP", "Monocyte", 
               "Kupffer", "cDC1", "cDC2", "pDC", "ASDC"), 
             min_value = 0.02)
  p3 <- plot(c("LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", 
               "Large-PreB",  "Small-PreB", "IM-B", "NK", "T", "ILCP"),
             min_value = 0.02, leg.pos = "bottom")
  leg <- ggpubr::get_legend(p3)
  p3 <- p3 + theme(legend.position = "none")
  ggsave(plot = gridExtra::grid.arrange(p1,p2,p3, ncol = 1),
         filename = paste0("plot/Moscot_summary-barplot_FL_blood_", direction, "_on3rows.pdf"), 
         width = 15, height = 15)
  
  p1 <- plot(c("HSC", "GP", "Granulocyte", "MEMP-t", "MEMP", "MEP"), min_value = 0.02)
  p2 <- plot(c("MEMP-Mast", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery", "MEMP-MK", "MK"),min_value = 0.02)
  p3 <- plot(c("MastP-t", "MastP", "Mast", "Monocyte", "Kupffer"),min_value = 0.02)
  p4 <- plot(c("MDP", "cDC1", "cDC2", "pDC", "ASDC"), min_value = 0.02)
  p5 <- plot(c("LMPP", "LP", "Cycling-LP", "PreProB"), min_value = 0.02)
  p6 <- plot(c("ProB-1", "ProB-2", "Large-PreB",  "Small-PreB", "IM-B", "NK", "T", "ILCP"), min_value = 0.02)
  
  ggsave(plot = gridExtra::grid.arrange(p1,p2,p3, p4,p5,p6, ncol = 1),
         filename = paste0("plot/Moscot_summary-barplot_FL_blood_", direction, "_on6rows.pdf"), 
         width = 12.2, height = 30)
  
  ggsave(plot = ggpubr::as_ggplot(leg), 
         filename = paste0("plot/Blood/LEGEND_Summary_barplot_", direction, ".pdf"))
  
}
