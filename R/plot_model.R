#'Plot model
#'
#' @param results_file : results
#' @param plot_file : file name to save plots
#' @param colours : a vector specifying the colour ramp
#' @return PDF of plots of results over timesteps
#' @import ggplot2
#' @export
#'
#' @examples plot_model(results_file = ODE_results, plot_file = "ODE_Plot.pdf", colours = c("royalblue3","orange", "red"))

plot_model <- function(results_file = ODE_results,
                       colours = c("royalblue3","orange", "red"),
                       plot_file = "ODE_Plot.pdf"){

  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]
  if(model_det$endo_species == 2){
    k <- (model_det$endo_species+1)^(model_det$endo_no_per_sp)-1
    cols <-  grDevices::colorRampPalette(colours)(k)
  }else{
    k <- (model_det$endo_species+1)^(model_det$endo_no_per_sp)-1
    cols <-  grDevices::colorRampPalette(colours)(k)
  }
  pdf(file = plot_file)
      for(i in all_results){
    title <- bquote("Number of individuals carrying parasites over time")
    if(model_det$endo_species == 2){
    subtit <- bquote(list(K ==. (i$Parameters$K),
                          lambda ==. (i$Parameters$lambda),
                          mu==.(i$Parameters$mu),
                          beta[A]==.(i$Parameters$betAa),
                          beta[B]==.(i$Parameters$betaB),
                          sigma[A]==.(i$Parameters$sigmaA),
                          sigma[B]==.(i$Parameters$sigmaB),
                          sigma[AB]==.(i$Parameters$sigmaAB),
                          sigma[BA]==.(i$Parameters$betaBA),
                          nu[A]==.(i$Parameters$nuA),
                          nu[B]==.(i$Parameters$nuB))
                     )
    }else{
      subtit <- bquote(list(K ==. (i$Parameters$K),
                            lambda ==. (i$Parameters$lambda),
                            mu==.(i$Parameters$mu),
                            beta[A]==.(i$Parameters$betaA),
                            sigma[A]==.(i$Parameters$sigmaA),
                            nu[A]==.(i$Parameters$nuA))

      )
    }
    col <- c(colnames(i$Results[,-1]))
    initial_plot <- tidyr::pivot_longer(i$Results, all_of(col), names_to = "infection_status")
    initial_plot$infection_status <- factor(initial_plot$infection_status, levels = col)
    print(ggplot(initial_plot)+
      geom_line(aes(x=time, y = value,
                    colour = infection_status,
                    group = infection_status,
                    linetype = infection_status), size = 1.5)+
      ylab(label="# Host species")+
      xlab(label="Timesteps")+
      labs(colour = "infection status \n         A  B",
           linetype = "infection status \n         A  B")+
      scale_colour_manual(values=c("black",  cols))+
      theme_classic()+
      theme(axis.text = element_text(size = 20))+
      theme(axis.title = element_text(size =20))+
      theme(legend.text = element_text(size = 20))+
      theme(legend.key.size = unit(3, "lines"))+
      theme(legend.title = element_text(size = 20))+
      ggtitle(bquote(atop(bold(.(title)),atop(bold(.(subtit)))))))
  }
  while (!is.null(dev.list()))
    dev.off()
}
