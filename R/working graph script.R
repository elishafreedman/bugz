#figures script

 #single endosymbionts

betaA_single_endo_no_dem <- area_plot(data = single_endo_no_dems_run_beta,
                               x_labs = expression(paste(beta)),
                               y_labs = NULL,
                               legend_labs = "# infected",
                               ind_var =  "betaA",
                               def_param_plot = 0.0005,
                               colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                               def_line_colour = "yellow",
                               line_type = 1,
                               legend_pos = "NULL",
                               defaults = c(sigmaA == 0 & nuA == 0.005),
                               xbreaks = seq(0.0001,0.001, 0.0001),
                               titles = NULL)


betaA_single_endo <- area_plot(data = single_endo_run ,
                                      x_labs = expression(paste(beta)),
                                      y_labs = NULL,
                                      legend_labs = "# infected",
                                      ind_var =  "betaA",
                                      def_param_plot = 0.0005,
                                      colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                                      def_line_colour = "yellow",
                                      line_type = 1,
                                      legend_pos = "NULL",
                                      defaults = c(sigmaA == 0 & nuA == 0.005),
                                      xbreaks = seq(0.0001,0.001, 0.0001),
                               title = NULL)


sigmaA_single_endo <- area_plot(data = single_endo_run ,
                                x_labs = expression(paste(sigma)),
                                y_labs = NULL,
                                legend_labs = "# infected",
                                ind_var =  "sigmaA",
                                def_param_plot =  0,
                                colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                                def_line_colour = "yellow",
                                line_type = 1,
                                legend_pos = "NULL",
                                defaults = c(betaA == 0.0005 & nuA == 0.005),
                                xbreaks = seq(0, 1, 0.1))

sigmaA_single_endo_no_dem <- area_plot(data = single_endo_no_dems_run_sigma ,
                                       x_labs = expression(paste(sigma)),
                                       y_labs = NULL,
                                       legend_labs = "# infected",
                                       ind_var =  "sigmaA",
                                       def_param_plot =  0,
                                       colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                                       def_line_colour = "yellow",
                                       line_type = 1,
                                       legend_pos = "NULL",
                                       defaults = c(betaA == 0.0005 & nuA == 0.005),
                                       xbreaks = seq(0, 1, 0.1))


nuA_single_endo <- area_plot(data = single_endo_run,
                             x_labs = expression(paste(nu)),
                             y_labs = NULL,
                             legend_labs = "# infected",
                             ind_var =  "nuA",
                             def_param_plot = 0.005,
                             colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                             def_line_colour = "yellow",
                             line_type = 1,
                             legend_pos = "NULL",
                             defaults = c(sigmaA == 0 & betaA == 0.0005),
                             xbreaks = seq(0.001, 0.01, 0.001))

nuA_single_endo_no_dem <- area_plot(data = single_endo_no_dems_run_nu,
                                    x_labs = expression(paste(nu)),
                                    y_labs = NULL,
                                    legend_labs = "# infected",
                                    ind_var =  "nuA",
                                    def_param_plot = 0.005,
                                    colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                                    def_line_colour = "yellow",
                                    line_type = 1,
                                    legend_pos = "NULL",
                                    defaults = c(sigmaA == 0 & betaA == 0.0005),
                                    xbreaks = seq(0.001, 0.01, 0.001))

#multiple plots

mult_graphs_single_endo <- ggpubr::ggarrange(betaA_single_endo ,
                                             betaA_single_endo_no_dem ,
                                             sigmaA_single_endo ,
                                             sigmaA_single_endo_no_dem ,
                                             nuA_single_endo ,
                                             nuA_single_endo_no_dem,
                                         ncol = 2, nrow = 3,
                                         labels = c("A", "D", "B", "E", "C", "F"),
                                        common.legend = TRUE, legend = "right", vjust = 1)
mult_graphs_single_endo

annotate_figure(mult_graphs_single_endo,
                left = text_grob("porportion infected", rot = 90, vjust = 1, size = 15),
                top = text_grob("With demographics                                                     Without demographics", face = "bold", size = 20))

library(ggpp)
