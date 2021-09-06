#' get_equil_points
#'
#' @param results_file : name of file where raw results have been stored. This file must be uploaded to the environment
#' @param eq_threshold: the percent threshold for coefficient of variance to determine the point of stable equilibrium
#' @param eq_var : the method for calculating variance through time.   choose from coefficient of variance (coef), percent variance (perVar), or the functions in "cor" function from the stats package.
#'
#' @return
#' @export
#'
#' @examples get_equil_points(results_file = ODE_results, eq_threshold = 0.01, eq_var = "perVar")

get_equil_points<-function(results_file = ODE_results, eq_threshold = 0.5,
                           eq_var = coef){

  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]


  colength <- rep(NA, ncol(all_results[[1]][["Results"]]) + ncol(parameters))
  eq_raw <- as.data.frame(matrix(NA, nrow = length(all_results),
                                 ncol = length(colength)))
  colnames(eq_raw) <- c(colnames(parameters),
                        colnames(all_results[[1]]$Results))

  #coefficient of variance
  coef <- function(X){
    (sd(X) / mean(X)) * 100
  }

  #percent of variance

  perVar <- function(x){
  x1 <- round(x[1], 2)
  x2 <- round(x[2], 2)
  max(abs(x1-x2)/x1, na.rm = TRUE) <eq_threshold
  }

  # finding the time system reaches stable equilibrium  and fill in the raw results



  print(paste("Finding timepoint where stable equilibrium is reached"))
  for (i in 1:length(all_results)){
    s <- all_results[[i]][["Results"]]
    p <- all_results[[i]][["Parameters"]]
    for (k in 1:(nrow(s) - 1)){
      eq <- as.matrix(s[k:(k + 1),])
      VarC <- apply(eq[, -1], 2, FUN = eq_var)
      #print(VarC)
      if(all(VarC == TRUE)){
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
        ")!"
      )
    }
  }

    return(list(
      "simulation_details" = model_det,
      "equilibrium_results" = eq_raw
    ))
  }


