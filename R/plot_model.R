#'Plot model
#'
#' @param results_file : results
#' @param plot_file : file name to save plots
#' @param Colors :  vector of colours to be included in the graph
#' @return PDF of plots of results over timesteps
#' @import ggplot2
#' @export
#'
#' @examples plot_model(results_file = ODE_results, plot_file = "ODE_Plot.pdf", Colors = c("blue","orange", "red"))

plot_model <- function(results_file = ODE_results,
                       plot_file = "ODE_Plot.pdf", Colors = c("blue","orange", "red")){

  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]
  if(model_det$endo_species == 2){
    k <- (model_det$endo_species+1)^(model_det$endo_no_per_sp)-1
    cols <-  grDevices::colorRampPalette(c(Colors))(k)
  }else{
    k <- (model_det$endo_species+1)^(model_det$endo_no_per_sp)-1
    cols <-  grDevices::colorRampPalette(c(Colors))(k)
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
    #col <- c(colnames(i$Results[,-1]))

    #recreate new ordered columns
    #species A
    specA <- rev(seq(1, endo_no_per_sp))
    colA <- rep(NA, length(specA))
    for(i in 1:length(specA)){
      colA[i] <- paste0("N", specA[i], 0)
    }

    #species B
    specB <- seq(1, endo_no_per_sp)
    colB <- rep(NA, length(specB))
    for(i in 1:length(specB)){
      colB[i] <- paste0("N", specB[i], 0)
    }

    #co-infected

    #coiA
    coiA1 <- rev(seq(1:endo_no_per_sp))
    coiA2 <- seq(1:endo_no_per_sp)
    coiA <- c(coiA1, coiA2)

    coiB <- c(coiA2, coiA1)

    coiAll <- rep(seq(1, endo_no_per_sp, 1), endo_no_per_sp*2)

    for(i in 1:length(endo_no_per_sp)){
    coiAll2 <- rep(1)
    }



    coiAll <- rep(NA, length(coiB))
     for(i in 1:length(coiAll)){
       coiAll[i] <- paste0("N",coiA[i],coiB[i])
     }

    #coiB
    coiB <- rev(coiA)

    coiAllB <- rep(NA, length(coiB))

    for(i in 1:length(coiAllB)){
      coiAllB[i] <- paste0("N",coiA[i],1)
    }


    #concatenate vectors

    col <- c("N00", specA, coiA, coidub, coiB, specB)


    initial_plot <- tidyr::pivot_longer(i$Results, all_of(col), names_to = "infection_status")
    initial_plot$infection_status <- factor(initial_plot$infection_status, levels = col)
    print(ggplot(initial_plot)+
      geom_line(aes(x=time, y = value,
                    colour = infection_status,
                    group = infection_status), size = 1.5)+
      ylab(label="# Host species")+
      xlab(label="Timesteps")+
      labs(colour = "infection status \n         A  B")+
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
