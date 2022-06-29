area_plot  <- function(data = NULL ,
                       x_labs = expression(paste(beta[A])),
                       y_labs = "% infected",
                       legend_labs = "Infection\n status \n       AB",
                       ind_var =  "betaA",
                       def_param_value = NULL,
                       colouring = c("dodgerblue4", "darkolivegreen3", "chartreuse2")){
  model_det <- data[["simulation_details"]]
  parameters <- data[["param_combos"]]
  all_results <- data[["simulations"]]

  k <- (model_det$endo_species + 1) ^ (model_det$endo_no_per_sp) - 1
  colours <- colorRampPalette(colouring)((k))



  col <- Reorder(res = all_results, mod_det = model_det)
  col_labs <- gsub("\\N", "", col)
  data1 <-
    all_results %>% pivot_longer(all_of(col)) %>% group_by(ind_var) %>% summarise(sum = sum(value))

  data2 <- cbind(data1, sum = data1$sum)

  data3 <- data2[, c(col)]

  data2 <- cbind(data2[, 1:12], data3, sum = data1$sum)

  #for use with single endosymbionts

  data_test <- data2 %>% gather(Group, value, all_of(col))

  data4 <- pivot_longer(data2, all_of(col), names_to = "Group")

  #levels(data4$value) <- coln_dub

  ind_var <- enquo(ind_var)

  #double endosymbionts

  data2 %>% pivot_longer(all_of(coln_dub), names_to = "Group") %>% factor(data2$"Group", levels = col) %>% ggplot(aes(x = !!ind_var, y = value /
                                                                                                                        sum, fill = Group)) +
    ggplot(data4, aes(x = betaA, y = value / sum, fill = Group)) +
    geom_area(size = 0.5, colour = "black") +
    theme_classic() +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 15)) +
    theme(axis.text = element_text(size = 12)) +
    theme(axis.title = element_text(size = 15)) +
    theme(axis.title.y = element_text(angle = 90, size = 15)) +
    ylab(label = y_labs) +
    xlab(label = x_labs) +
    scale_fill_manual(values = c("black", colours), labels = labs) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = seq(0.0001, 0.001, 0.0001)) +
    theme(legend.position = "right") +
    theme(element_line(size = 1)) +
    geom_vline(xintercept = def_param_value,
               linetype = 1,
               colour = "white") +
    labs(fill = legend_labs)
}
