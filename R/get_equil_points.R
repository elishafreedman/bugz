#' Title
#'
#' @param results_file
#' @param eq_threshold
#' @param eq_var
#'
#' @return
#' @export
#'
#' @examples
get_equi_points<- function(results_file = OD_results, eq_threshold = 0.5,
                           eq_var = coef){
  #coefficient of variance
  coef <- function(X) {
    (sd(X) / mean(X)) * 100
  }

  #percent of variance

  perVar <- function(X) {
    abs((X[2] - X[1])/X[1])
  }

  # finding the time system reaches stable equilibrium  and fill in the raw results



  print(paste("Finding timepoint where stable equilibrium is reached"))
  for (i in 1:length(to_test)){
    s <- to_test[[i]][["Results"]]
    p <- to_test[[i]][["Parameters"]]
    for (k in 1:(nrow(s) - 1)){
      eq <- as.matrix(s[k:(k + 1),])
      VarC <- apply(eq[, -1], 2, FUN = eq_var)
      #print(VarC)
      if (all(VarC <= eq_threshold, na.rm = T)){
        eq_raw[i, ] <- c(p, s[k + 1, ])
        break
      }
    }
    print(paste("search for equilbrium at combination", i, "done"))
    if (all(is.na(eq_raw))) {
      warning(
        "system did not reach a stable equilbrium at maximum timestep
          (",
        model_det$max_timesteps,
        ")! \n statistics at equilibrium will not be calculated"
      )
    }
  }
  #prep file for saving
  if (any(cor_param_method == "ld")) {
    return(list(
      model_det = model_det,
      "equilibrum_stats" = list(
        "Results" = eq_raw,
        "Proportion" = eq_av_prop ,
        "correlations" = eq_param_cor,
        "ld" = eq_ld
      )
    ))
  } else{
    return(list(
      "model_det" = model_det,
      "equilibrum_stats" = list(
        "Results" = eq_raw,
        "Proportion" = eq_av_prop ,
        "correlation" = eq_param_cor
      )
    ))
  }

}
}
}
