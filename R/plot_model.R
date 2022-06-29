#'Plot model
#'
#' @param results_file : results
#' @param plot_file : file name to save plots
#' @param Colors :  vector of colours to be included in the graph
#' @return PDF of plots of results over time steps
#' @import ggplot2
#' @export
#'
#' @examples plot_model(results_file = ODE_results, plot_file = "ODE_Plot.pdf", Colors = c("blue","orange", "red"))

plot_model <- function(results_file = ODE_results,
                       plot_file = "ODE_Plot.pdf",
                       Colors = c("blue", "orange", "red")) {
  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]
  if (model_det$endo_species == 2) {
    k <- (model_det$endo_species + 1) ^ (model_det$endo_no_per_sp) - 1
    cols <-  grDevices::colorRampPalette(c(Colors))(k)
  } else{
    k <- (model_det$endo_species + 1) ^ (model_det$endo_no_per_sp) - 1
    cols <-  grDevices::colorRampPalette(c(Colors))(k)
  }
  pdf(file = plot_file)
  for (i in 1:nrow(all_results)) {
    title <-
      bquote("Number of individuals carrying parasites over time")
    if (model_det$endo_species == 2) {
      subtit <- bquote(
        list(
          K == . (all_results[i, ]$K),
          lambda == . (all_results[i, ]$lambda),
          mu == .(all_results[i, ]$mu),
          beta[A] == .(all_results[i, ]$betAa),
          beta[B] == .(all_results[i, ]$betaB),
          sigma[A] == .(all_results[i, ]$sigmaA),
          sigma[B] == .(all_results[i, ]$sigmaB),
          sigma[AB] == .(all_results[i, ]$sigmaAB),
          sigma[BA] == .(all_results[i, ]$betaBA),
          nu[A] == .(all_results[i, ]$nuA),
          nu[B] == .(all_results[i, ]$nuB)
        )
      )
    } else{
      subtit <- bquote(list(
        K == . (all_results[i, ]$K),
        lambda == . (all_results[i, ]$lambda),
        mu == .(all_results[i, ]$mu),
        beta[A] == .(all_results[i, ]$betaA),
        sigma[A] == .(all_results[i, ]$sigmaA),
        nu[A] == .(all_results[i, ]$nuA)
      ))
    }

    col <- Reorder(res = all_results, mod_det = model_det)
    col_labs <- gsub("\\N","",col)


    initial_plot <-
      tidyr::pivot_longer(i$Results, all_of(col), names_to = "infection_status")
    initial_plot$infection_status <-
      factor(initial_plot$infection_status, levels = col)
    print(
      ggplot(initial_plot) +
        geom_line(
          aes(
            x = time,
            y = value,
            colour = infection_status,
            group = infection_status
          ),
          size = 1.5
        ) +
        ylab(label = "# Host species") +
        xlab(label = "Timesteps") +
        labs(colour = "infection status \n         A  B") +
        scale_colour_manual(values = c("black",  cols),labels =  col_labs) +
        theme_classic() +
        theme(axis.text = element_text(size = 20)) +
        theme(axis.title = element_text(size = 20)) +
        theme(legend.text = element_text(size = 20)) +
        theme(legend.key.size = unit(3, "lines")) +
        theme(legend.title = element_text(size = 20)) +
        ggtitle(bquote(atop(
          bold(.(title)), atop(bold(.(subtit)))
        )))
    )
  }
  while (!is.null(dev.list()))
    dev.off()
}
