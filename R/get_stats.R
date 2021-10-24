#' get_stats
#'

#' @param results_file : name of file where raw results have been stored. This file must be uploaded to the environment
#'
#' @param test_parameters : vector of parameters tested.  if NA, all parameters will be tested apart from designated baseline.
#' additionally, parameters where values are to be matched during analysis can be specified with "etc == etc". Parameters that are able to do defined in this way are: nu, sigma, and beta.
#' @param set_baseline : vector of parameters to be kept at baseline, if NA, all parameters combinations will be considered.
#' @param cor_param_method : methods of correlation calculations, options include ld + Dprime and r2, and or
#' any method from the "cor" function in the stats package.if ld is selected Dprime and r2 will be automatically calculated.
#' @param eq_threshold : the percent threshold for coefficient of variance to determine the point of stable equilibrium
#'
#' @return A list containing the following: details of the model including the date of
#' simulation and parameter combinations, raw results at equilibrium, correlation tests at equilibrium for each test parameter, and if ld is selected linkage disequilibrium at each parameter combination
#' @export
#'
#'
#' @examples get_stats(results_file = ODE_results,test_parameters = c("sigmaBA = sigmaAB"),set_baseline = c(K = 200, lambda = 1, mu = 0.5),cor_param_method = c("ld", "pearson"), eq_threshold = 0.5,core_spec = 2)


get_stats <- function(results_file = ODE_results,
                      test_parameters = c("sigmaBA = sigmaAB", "nuA = nuB"),
                      set_baseline = c(K = 200, lambda = 1, mu = 0.5),
                      cor_param_method = c("ld", "pearson")){
  #split for ease of indexing
  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  all_results <- results_file[["simulations"]]


  #subset  the baseline parameters
  if (length(set_baseline >= 1) && is.na(set_baseline) == FALSE){
    to_testB <-
      lapply(all_results, function(x)
        x[all(set_baseline %in% x$Parameters)])
    to_testB <- to_testB[lapply(to_testB, length) > 0]
  } else{
    to_testB <- all_results
  }
  if (length(to_testB) == 0){
    warning("Baseline values for parameters you would like to test were not simulated!")
  }
  #isolating parameter combinations of interest from the dataset
    if(model_det$endo_species >= 2){
    #subset the matched parameters from the test_parameters arguments
    #create vector to  match
    dplyr::case_when(test_parameters == "sigmaAB = sigmaBA" ~ "sigmaAB",
      test_parameters == "sigmaBA = sigmaAB" ~ "sigmaAB",
      test_parameters == "sigmaA = sigmaB" ~ "sigmaA",
      test_parameters == "sigmaB = sigmaA" ~ "sigmaA",
      test_parameters == "betaB = betaA" ~ "betaA",
      test_parameters == "betaA = betaB" ~ "betaA",
      test_parameters == "nuA = nuB" ~ "nuA",
      test_parameters == "nuB = nuA" ~ "nuA")
    }
    print(test_parameters)
    to_test <- lapply(to_testB, function(x)
        x[unique(x$Parameters[, test_parameters]) %in% x$Parameters])
    to_test <- to_test[lapply(to_test, length) > 0]
  } else{
    to_test <- to_testB
  }
  if (length(to_test) == 0){
    warning("Parameter combinations you would like to test were not simulated!")
  }

  # stats at equilibrium #

  # raw results for each parameter at equilibrium
  colength <- rep(NA, ncol(to_test[[1]][["Results"]]) + ncol(parameters))
  eq_raw <- as.data.frame(matrix(NA, nrow = length(to_test),
                                 ncol = length(colength)))
  colnames(eq_raw) <- c(colnames(parameters),
                        colnames(to_test[[1]]$Results))

  #fill in parameters and and  infection status
  if (model_det$endo_species >= 2){
    eq_infect <- c("N00", "A", "B", "A_plus", "B_plus", "coinf")
  } else{
    eq_infect <- c("N0", "A", "A_plus")
  }



  # function to calculate the proportion at equilibrium
  prop <- function(X){
    X / (p$K * p$mu)
  }


  #proportion co-infected
  if (model_det$endo_species >= 2){
    coinf <-
      c(rowSums(eq_raw[, grep("^[^0slnbtKm]*$", colnames(eq_raw))]))
    #create vector that will be used later on for the pivot longer
    #prop single B
    B <- prop(eq_raw$N01)
    # double B
    #string pattern needed to recognise columns for species B co-infection
    if (model_det$endo_species >= 2){
      patB <- rep(NA, model_det$endo_no_per_sp)
      for (i in 2:model_det$endo_no_per_sp){
        patB[i] <- c(paste0("N0", i))
      }
      patB <- na.omit(patB)
      doubleB_prop <- c(prop(rowSums(eq_raw[eq_raw %in% patB])))
    }
  }
  #prop single A
  A <- prop(eq_raw$N10)


  # single species co-infections

  #string pattern needed to recognise columns for species A co-infections
  patA <- rep(NA, model_det$endo_no_per_sp)
  if (model_det$endo_species >= 2) {
    for (i in 2:model_det$endo_no_per_sp) {
      patA[i] <- c(paste0("N", i, "0"))
    }
  } else{
    for (i in 2:model_det$endo_no_per_sp)
      patA[i] <- c(paste0("N", i))
  }
  patA <- na.omit(patA)
  doubleA_prop <- c(prop(rowSums(eq_raw[eq_raw %in% patA])))

  #prop_uninfected

  if (model_det$endo_species >= 2) {
    N00 <- c(prop(eq_raw[, "N00"]))
    eq_dat <- cbind(
      eq_raw,
      data.frame(
        "N00" = N00,
        "A" = A,
        "B" = B,
        "A_plus" = doubleA_prop,
        "B_plus" = doubleB_prop,
        "coinf" = coinf
      )
    )
  } else{
    N0 <- c(prop(eq_raw[, "N0"]))
    eq_dat <- cbind(eq_raw, data.frame(
      "N0" = N0,
      "A" = A,
      "A_plus" = doubleA_prop
    ))
  }



  # equilibrium proportions
  if (all(is.na(test_parameters)) == TRUE) {
    par <- colnames(parameters)
    eq_met <-tidyr::pivot_longer(eq_dat, all_of(par), names_to = "parameter")
  } else{
    eq_met <-tidyr::pivot_longer(eq_dat, all_of(test_parameters), names_to = "parameter")
  }

  eq_met <- tidyr::pivot_longer(eq_met,
                         eq_infect,
                         names_to = "infection_status",
                         values_to = "proportion")
  eq_av_prop <- eq_met[, c("parameter", "value", "infection_status", "proportion")]
  eq_av_prop <- eq_av_prop[order(eq_av_prop$parameter, decreasing = FALSE),]

  #calculating correlation coefficients


  # linkage disequilibrium

  ld <- function(coinf,
                 N0,
                 A,
                 B,
                 A_plus = 0,
                 B_plus = 0){
      D <-coinf * N0 - sum(A, A_plus) * sum(B, B_plus)

      #D prime

      if(D >= 0){
        Dmax <- min(A*B, -(1-A)*(1-B))
        Dp <- D/Dmax
      }
      if(D < 0){
        Dmin <- min(A*(1-B), B*(1-A)*B)
        Dp <- D/Dmin
      }

      #r2
     r <- D^2/A*(1-A)*B*(1-B)

      #phi coefficient
     ph <- coinf*N0-A-B/sqrt(A+A_plus*B+B_plus*N0)
     c(D, Dp, r, ph)
    }


  #add levels to dataset

  eq_av_prop$parameter <- factor(eq_av_prop$parameter)
  eq_av_prop$infection_status <- factor(eq_av_prop$infection_status)


  print(paste("Calculating correlations at equilibrium"))

  #calculate correlation coefficients.
  if (length(cor_param_method) >= 1){

    #linkage disequilibrium related stats
    if (any(cor_param_method == "ld")){
      corT <- cor_param_method[grep("^[^ld]", cor_param_method)]
      eq_param_cor <-eq_av_prop |> dplyr::group_by(parameter, infection_status) |> dplyr::summarise(cor = cor(value, proportion, method = corT))
      eq_ld <- data.frame("ld" = NA, "Dprime", "r2", "phi_coef")
      cols <- c(colnames(to_test[[1]]$Parameters))
      eq_p <-data.frame(matrix(ncol = length(to_test[[1]]$Parameters)))
      colnames(eq_p) <- cols
      for (i in 1:length(to_test)) {
        eq_p[i, ] <- to_test[[i]][["Parameters"]]
        if (model_det$endo_species >= 2) {
          eq_ld[i, ] <- ld(
            coinf = eq_dat$coinf[i],
            N0 = eq_dat$N00[i],
            A = eq_dat$A[i],
            B = eq_dat$B[i],
            A_plus = eq_dat$A_plus[i],
            B_plus = eq_dat$B_plus[i]
          )

        } else{
          eq_ld[i, ]<- ld(
            coinf = eq_dat$coinf[i],
            N0 = eq_dat$N0[i],
            A = eq_dat$A[i],
            B = eq_dat$B[i],
            A_plus = eq_dat$A_plus[i],
            B_plus = eq_dat$B_plus[i]
          )
        }
      }
      eq_ld <- cbind(eq_p, eq_ld)
    } else{
      eq_param_cor <- eq_av_prop |> dplyr::group_by(parameter, infection_status) |> dplyr::summarise(cor = cor(value, proportion, method = cor_param_method))
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
