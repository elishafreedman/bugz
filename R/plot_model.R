#' Create  area plots for simulations generated in "bugz"
#'
#' @param data dataset
#' @param y_labs Label for y axis
#' @param Legend_labs Label for legend
#' @param ind_var desired independent variable
#' @param def_param_plot default parameters to plot
#' @param colouring vector for the colour gradient desired in the plot
#' @param def_line_colour the colour of the line representing the default parameter
#' @param line_type line type of default parameter
#' @param legend_pos legend position
#' @param defaults the default parameter values for the dataset
#' @param xbreaks the breaks for the tick marks and number on y axis
#' @param ybreaks the breaks for the tick marks and numbers on the x axis
#' @param titles title for the graph
#' @param second_def_line to add another line displaying a second default parameter to the plot
#' @param second_line_colour colour for this second default parameter line
#' @param tags tags for the plot
#' @return a 100% stacked area plot of results
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import forcats
#' @examples
#' @export


plot_model <- function(data = NULL ,
                       x_labs = expression(paste(sigma[A])),
                       y_labs = "% infected",
                       legend_labs = "Infection\n status \n       AB",
                       ind_var =  "betaA",
                       def_param_plot = NULL,
                       colouring = c("dodgerblue4", "darkcyan", "cadetblue1"),
                       alpha = 1,
                       uninfected_colour = "white",
                       def_line_colour = "grey",
                       line_type = 1,
                       legend_pos = "NULL",
                       defaults = c(sigmaA == 0 & nuA == 0.005),
                       xbreaks = NULL,
                       ybreaks = NULL,
                       titles = NULL,
                       second_def_line = NULL,
                       second_line_colour = NULL,
                       tags  = NULL,
                       tag_pos = NULL,
                       ticks = NULL) {
  #controlling which x axis ticks are labeled
  # label_at <- function(n)
  #   function(x) ifelse(x %% n == 0, x, "")

  model_det <- data[["simulation_details"]]
  all_results <- data[["simulations"]]

  k <- (model_det$endo_no_per_sp + 1) ^ (model_det$endo_species)
  colours <- colorRampPalette(colouring)((k))

  colours <- c(uninfected_colour, colours)

  #xbreaks <- dplyr::enquo(xbreaks)

  #defaults <- dplyr::enquo(defaults)
  if (is.na(defaults == TRUE)) {
    all_res <- all_results
  } else{
    #all_res <-  all_results |> dplyr::filter(!!defaults)
    all_res <- subset(all_results, eval(parse(defaults)))
  }
  col <- Reorder(res = all_res, mod_det = model_det)


  col_labs <- gsub("\\N", "", col)
  if (model_det$endo_species == 1) {
    col_labs <- gsub('.{1}$', '', col_labs)
  }
  datasum <-
    all_res |> tidyr::pivot_longer(all_of(col)) |>
    group_by(.data[[ind_var]]) |>
    summarise(sum = sum(value))


  data2 <- cbind(all_res, sum = datasum$sum)

  data2 |> tidyr::pivot_longer(all_of(col)) |>
    mutate(name = forcats::fct_relevel(name, all_of(col))) |>
    ggplot(aes(x = .data[[ind_var]], y = value / sum, fill = name)) +

    geom_area(size = 0.5,
              colour = "black",
              alpha = alpha) +
    scale_fill_manual(values = colours, labels = col_labs) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = ybreaks) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = xbreaks, labels = ticks) +
    geom_vline(
      xintercept = def_param_plot,
      linetype = line_type,
      colour = def_line_colour,
      size = 1
    ) +
    geom_vline(
      xintercept = second_def_line,
      linetype = line_type,
      colour = second_line_colour,
      size = 1
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 15),
      axis.title.y = element_text(angle = 90, size = 20),
      axis.title.x = element_text(size = 20),
      plot.title = element_text(
        color = "Black",
        size = 20,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = margin(0.7, 0.7, 0.7, 0.7, "cm"),
      legend.position = legend_pos,
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 20),
      plot.tag.position = tag_pos,
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 20),
      panel.grid = element_blank()
    ) +
    xlab(x_labs) +
    ylab(y_labs) +
    labs(fill = legend_labs, tag = tags) +
    ggtitle(titles)
}
