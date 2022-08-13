area_plot  <- function(data = NULL ,
                       x_labs = expression(paste(beta[A])),
                       y_labs = "% infected",
                       legend_labs = "Infection\n status \n       AB",
                       tags = NULL,
                       ind_var =  "betaA",
                       def_param_plot = NULL,
                       colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                       def_line_colour = "grey",
                       line_type = 1,
                       legend_pos = "NULL",
                       defaults = c(sigmaA == 0.05 & nuA == 0.005),
                       xbreaks = NULL,
                       ybreaks = NULL,
                       titles = NULL,
                       other_param_plot = NULL,
                       other_line_colour = NULL
                       ){
  model_det <- data[["simulation_details"]]
  all_results <- data[["simulations"]]

  k <- (model_det$endo_no_per_sp+1) ^ (model_det$endo_species)
  colours <- colorRampPalette(colouring)((k))

  colours <- c("black", colours)


   defaults <- dplyr::enquo(defaults)

   xbreaks <- dplyr::enquo(xbreaks)

  col <- Reorder(res = all_results, mod_det = model_det)




  col_labs <- gsub("\\N", "", col)
  if(model_det$endo_species == 1){
    col_labs <- gsub('.{1}$', '', col_labs)
  }
 all_res <-  all_results %>% dplyr::filter(!!defaults)


datasum <- all_res %>% pivot_longer(all_of(col)) %>% group_by(.data[[ind_var]]) %>% summarise(sum = sum(value))


 data2 <- cbind(all_res, sum = datasum$sum)


data2 %>% pivot_longer(all_of(col)) %>% mutate(name = forcats::fct_relevel(name, all_of(col)))%>% ggplot(aes(x = .data[[ind_var]], y = value/sum, fill = name))+
    geom_area(size = 0.5, color = "black")+
    scale_fill_manual(values = colours, labels = col_labs) +
    scale_y_continuous(expand = c(0, 0), breaks = ybreaks)+
    scale_x_continuous(expand = c(0, 0),
                       breaks = xbreaks) +
    geom_vline(xintercept = def_param_plot,
               linetype = line_type,
               colour = def_line_colour, size = 1)+
  geom_vline(xintercept = other_param_plot,
             linetype = line_type,
             colour = other_line_colour, size = 1)+
  theme_bw()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.title.y = element_text(angle = 90, size = 15),
        axis.title.x = element_text(size = 20),
        plot.title = element_text(color="Black", size=12, face="bold"),
                                  plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"),
        legend.position = legend_pos,
        legend.title = element_text(size=15),
        plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0.85, 0.05), legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 15))+
  ylab(label = y_labs) +
  xlab(label = x_labs)+
  labs(fill = legend_labs, tag = tags)+
  ggtitle(titles)
}

#combine plots



mult_graphs_no_dem <- ggpubr::ggarrange(plots = plot_comb,
                                        labels = c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3, common.legend = TRUE, legend = "right")
