get_equi_points<- function(results = Res,param = parameters, eq_t = eq_threshold){


  isEqui <- function(x, eq_thresh = eq_t){
    x1 <- round(x[1], 2)
    x2 <- round(x[2], 2)
    max(abs(x1-x2)/x1, na.rm = TRUE)<= eq_thresh
  }
  # finding the time system reaches stable equilibrium  and fill in the raw results

  eq_raw <- data.frame(matrix(nrow = nrow(param),
                              ncol = ncol(param)+ncol(all_results[[1]][["Results"]])))

  colnames(eq_raw) <- c(colnames(param), colnames(all_results[[1]][["Results"]]))

  print(paste("Finding timepoint where stable equilibrium is reached"))

  for(i in 1:length(results)){
    s <- results[i,]
    p <- param[i,]
    for (k in 1:(nrow(s) - 1)){
      eq <- as.matrix(s[k:(k + 1),])
      VarC <- apply(eq[, -1], 2, FUN = isEqui)
      if (all(VarC) == TRUE){
        eq_raw[i, ] <- c(p, s[k + 1, ])
        break
      }
    }
  }

    # if (all(is.na(eq_raw))){
    #   warning(
    #     "system did not reach a stable equilbrium at maximum timestep
    #       (",
    #     model_det$max_timesteps,
    #     ")! \n statistics at equilibrium will not be calculated"
    #   )
    # }
  }




