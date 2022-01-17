#'Plot model stats
#'
#' @param stats_file : results
#' @param plot_file : file name to save plots
#' @param param_comb : parameters you would like to be graphed against
#' @param stastic : the the statistic you would like to be plotted
#' @param Colors :  vector of colours to be included in the graph
#' @return PDF of plots of results over timesteps
#' @import ggplot2
#' @export
#'
#' @examples plot_model(results_file = ODE_results, plot_file = "ODE_Plot.pdf", Colors = c("blue","orange", "red"))

# stats plots

stats_plot <- function(stats_file = NA,
                       param_comb = NA,
                       statistic = NA,
                       plot_file = NA,
                       Colors = NA){



 Plots <- ggplot(stats_file, )
}
