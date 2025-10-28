
#' my_stacked_barplot
#' @description 
#' @param data a data.frame (rownames will be in x) or melted data.frame
#' @param melt if F data are already melted with colnames sample (x), variable (y), value
#' @param factor_levels_variables levels from top to bottom
#' @param col from top to bottom
#' @param min_value minimal value to print on the barplot
#' @param shift_text shift text along y axis to fit nicely the plot
my_stacked_barplot <- function(data, melt = F,
                               factor_levels_variables = NULL, col = NULL,
                               factor_level_samples = NULL,
                               min_value = 0, v.just = 0.5,
                               size_text = 3, col_text = "white"){
  if(melt){
    data$sample <- rownames(data)
    data <- reshape2::melt(data, id.vars = "sample")
  }
  
  if(!is.null(factor_levels_variables)){
    data$variable <- factor(data$variable, 
                            levels = factor_levels_variables)
  }else{
    data$variable <- as.factor(data$variable)
  }
  data <- data[order(data$sample),]
  
  if(!is.null(factor_level_samples)){
    data$sample <- factor(data$sample, levels = factor_level_samples)
  }
  
  
  p <- ggplot(data, aes(x = sample, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + mashaGgplot2Theme 
  
  data$to_write <- factor(data$value > min_value, levels = c("TRUE", "FALSE"))
  p <- p +
    geom_text(data, 
              mapping = aes(sample, value, label = value, col = to_write),
              position = position_stack(vjust = v.just), col = col_text,
              fontface = "bold", size = size_text)
  
  if(min_value != 0){
    p <- p +
      scale_color_manual(values = c(col_text, "#00000000"))
  }
  
  if(!is.null(col)){
    p <- p + scale_fill_manual(values = col)
  }
  
  return(p)
}

# Plot look and feel
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
