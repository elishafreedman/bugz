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
                               titles = NULL, tags = NULL)


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
                               title = NULL, tags = NULL)


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
                                xbreaks = seq(0, 1, 0.1), tags = NULL)

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
                                       xbreaks = seq(0, 1, 0.1), tags = NULL)


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
                             xbreaks = seq(0.001, 0.01, 0.001), tags = NULL)

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
                                    xbreaks = seq(0.001, 0.01, 0.001), tags = NULL)

#multiple plots

list_plot <- list(betaA_single_endo ,
               betaA_single_endo_no_dem ,
               sigmaA_single_endo ,
               sigmaA_single_endo_no_dem ,
               nuA_single_endo ,
               nuA_single_endo_no_dem)


mult_graphs_single_endo <- ggpubr::ggarrange(betaA_single_endo ,
                                             betaA_single_endo_no_dem ,
                                             sigmaA_single_endo ,
                                             sigmaA_single_endo_no_dem ,
                                             nuA_single_endo ,
                                             nuA_single_endo_no_dem,
                                         ncol = 2, nrow = 3,
                                         labels = c("A", "B", "C", "D", "E", "F"),
                                        common.legend = TRUE, legend = "right", vjust = 1)
mult_graphs_single_endo
annotate_figure(mult_graphs_single_endo,
                left = text_grob("porportion infected", rot = 90, vjust = 1, size = 12), top = text_grob("With demographics                                                                                                                   Without demographics",vjust = 1, size = 12, face = "bold"))

#export as 10 inch by 15 inch pdf


###two species graphs

two_species_betaA_eq_betaB1 <- area_plot(data = two_species_params_beta_variations_run,
                                      x_labs = expression(paste(beta[A]==beta[B])),
                                      y_labs = NULL,
                                      legend_labs = "# infected\n       AB",
                                      ind_var =  "betaA",
                                      def_param_plot = 0.0005,
                                      colouring = c("brown1","lavender","deepskyblue"),
                                      def_line_colour = "yellow",
                                      line_type = 1,
                                      legend_pos = NULL,
                                      defaults = c(betaB == betaA),
                                      xbreaks = seq(0.0001,0.001, 0.0001),
                                      titles = NULL, tags = NULL, other_param_plot = 0.0002, other_line_colour = "white")

two_species_betaA_over_betaB <- area_plot(data = two_species_params_beta_variations_run,
                                        x_labs = expression(paste(beta[A])),
                                        y_labs = NULL,
                                        legend_labs = "# infected\n       AB",
                                        ind_var =  "betaA",
                                        def_param_plot = 0.0005,
                                        colouring = c("brown1","lavender","deepskyblue"),
                                        def_line_colour = "yellow",
                                        line_type = 1,
                                        legend_pos = NULL,
                                        defaults = betaB ==0.0002,
                                        xbreaks = seq(0.0001,0.001, 0.0001),
                                        titles = NULL, tags = NULL, other_param_plot = 0.0002, other_line_colour = "white")


#sigma


two_species_sigmaA_eq_sigmaB_graph_UNEQ <- area_plot(data = two_species_params_sigmaA_eq_sigmaB_uneq_run,
                                                       x_labs = expression(paste(sigma[A]==sigma[B])),
                                                       y_labs = NULL,
                                                       legend_labs = "# infected\n       AB",
                                                       ind_var =  "sigmaB",
                                                       def_param_plot = 0.5,
                                                       colouring = c("brown1","lavender","deepskyblue"),
                                                       def_line_colour = "yellow",
                                                       line_type = 1,
                                                       legend_pos = NULL,
                                                       defaults = c(sigmaB ==sigmaA),
                                                       xbreaks = seq(0,1, 0.1),
                                                       titles = NULL, tags = NULL, other_param_plot = 0.2, other_line_colour = "white")


Two_species_sigmaAB_eq_BA_uneq <- area_plot(data = two_species_params_sigmaAB_run,
                                    x_labs = expression(paste(sigma[AB]==sigma[BA])),
                                        y_labs = NULL,
legend_labs = "# infected\n       AB",
ind_var =  "sigmaAB",
def_param_plot = 0.5,
colouring = c("brown1","lavender","deepskyblue"),
def_line_colour = "yellow",
line_type = 1,
legend_pos = NULL,
defaults = c(sigmaBA ==sigmaAB),
xbreaks = seq(0,1, 0.1),
titles = NULL, tags = NULL, other_param_plot = 0.2, other_line_colour = "white")


two_species_sigmaAB_eq_BA_eq <- area_plot(data = two_species_params_sigmaAB_run_eq,
                                            x_labs = expression(paste(sigma[AB]==sigma[BA])),
                                            y_labs = NULL,
                                            legend_labs = "# infected\n       AB",
                                            ind_var =  "sigmaAB",
                                            def_param_plot = 0.5,
                                            colouring = c("brown1","lavender","deepskyblue"),
                                            def_line_colour = "yellow",
                                            line_type = 1,
                                            legend_pos = NULL,
                                            defaults = c(sigmaBA ==sigmaAB),
                                            xbreaks = seq(0,1, 0.1),
                                            titles = NULL, tags = NULL, other_param_plot = 0.2, other_line_colour = "white")

two_species_sigmaAll_eq_BA_eq <- area_plot(data = two_species_params_sigmaAll_same_run_eq,
                                          x_labs = expression(sigma[AB]==sigma[BA]),
                                          y_labs = NULL,
                                          legend_labs = "# infected\n       AB",
                                          ind_var =  "sigmaAB",
                                          def_param_plot = 0.5,
                                          colouring = c("brown1","lavender","deepskyblue"),
                                          def_line_colour = "yellow",
                                          line_type = 1,
                                          legend_pos = NULL,
                                          defaults = c(sigmaBA ==sigmaAB),
                                          xbreaks = seq(0,1, 0.1),
                                          titles = NULL, tags = NULL, other_param_plot = 0.2, other_line_colour = "white")



mult_graphs_two_endo <- ggpubr::ggarrange(two_species_betaA_eq_betaB1 ,
                                             two_species_betaA_over_betaB,
                                          two_species_sigmaA_eq_sigmaB_graph,
                                             two_species_sigmaA_eq_sigmaB_graph_UNEQ,
                                             Two_species_sigmaAB_eq_BA_uneq,
                                          two_species_sigmaAll_eq_BA_eq,
                                             ncol = 2, nrow = 3,
                                             labels = c("A", "B", "C", "D", "E", "F"),
                                             common.legend = TRUE, legend = "right", vjust = 1)

mult_graphs_two_endo_4 <- ggpubr::ggarrange(two_species_betaA_eq_betaB1 ,
                                          two_species_betaA_over_betaB,
                                          Two_species_sigmaAB_eq_BA_uneq,
                                          two_species_sigmaAll_eq_BA_eq,
                                          ncol = 2, nrow = 2,
                                          labels = c("A", "B", "C", "D"),
                                          common.legend = TRUE, legend = "right", vjust = 1)
annotate_figure(mult_graphs_two_endo,
                left = text_grob("porportion infected", rot = 90, vjust = 1, size = 20))


annotate_figure(mult_graphs_two_endo,
                left = text_grob("porportion infected", rot = 90, vjust = 1, size = 12), top = text_grob("With demographics                                                                                                                   Without demographics",vjust = 1, size = 12, face = "bold"))

