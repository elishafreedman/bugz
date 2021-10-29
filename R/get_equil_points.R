#' get_equil_points
#'
#' @param results_file : file containing results.
#' @param eq_threshold : the threshold of variation between each time point at different.
#' @param eq_var : variance test for comparison.
#' @return
#' @export
#'
#' @examples params <- set_parameters(sigmaA = 0.1, sigmaAB = 0.1, sigmaB = 0.1)
#' ODE_results <- run_model()
#' equi <- function(results_file = ODE_results, eq_threshold = 0.5,eq_var = coef)
get_equi_points<- function(results_file = ODE_results, eq_threshold = 0.5,
                           eq_var = coef){
  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]
  #coefficient of variance
  coef <- function(X){
    (sd(X) / mean(X)) * 100
  }

  #percent of variance

  perVar <- function(X){
   abs((X[2] - X[1])/X[1])
  }

  # finding the time system reaches stable equilibrium  and fill in the raw results

  eq_raw <- data.frame()

  print(paste("Finding timepoint where stable equilibrium is reached"))
  for(i in 1:length(all_results)){
    s <- to_test[[i]][["Results"]]
    p <- to_test[[i]][["Parameters"]]
    for (k in 1:(nrow(s) - 1)){
      eq <- as.matrix(s[k:(k + 1),])
      VarC <- apply(eq[, -1], 2, FUN = eq_var)
      if (all(VarC <= eq_threshold, na.rm = T)){
        eq_raw[i, ] <- c(p, s[k + 1, ])
        break
      }
    }
    print(paste("search for equilbrium at combination", i, "done"))
    if (all(is.na(eq_raw))){
      warning(
        "system did not reach a stable equilbrium at maximum timestep
          (",
        model_det$max_timesteps,
        ")! \n statistics at equilibrium will not be calculated"
      )
    }
  }
  #prep file for saving
    return(list(
      model_det = model_det,
      parameters = parameters,
      Equilbrium = eq_raw
      )
    )
  }




