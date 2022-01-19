#' get_stats
#'

#' @param results_file : name of file where raw results have been stored. This file must be uploaded to the environment
#'
#' @param test_parameters : vector of parameters tested.  if NA, all parameters will be tested apart from designated baseline.
#' additionally, parameters where values are to be matched during analysis can be specified with "etc == etc". parameters that are able to do defined in this way are: nu, sigma, and beta.
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
#' @examples get_stats(results_file = ODE_results,test_parameters = c("sigmaBA = sigmaAB"),set_baseline = c(K = 200, lambda = 1, mu = 0.5),cor_param_method = c("ld", "pearson"),core_spec = 2)


get_stats <- function(results_file = ODE_eq,
                      test_parameters = c("sigmaBA" == "sigmaAB", "nuA" == "nuB"),
                      set_baseline = c(K = 200, lambda = 1, mu = 0.5),
                      cor_param_method = c("ld", "pearson")){

  #split for ease of indexing
  model_det <- results_file[["model_det"]]
  parameters <- results_file[["parameters"]]
  to_test <- results_file[["equilibrium"]]

  #subset  the baseline parameters


  if (length(set_baseline >= 1) && is.na(set_baseline) == FALSE){

     subset(parameters, parameters %in% set_baseline)
     subset(to_test, to_test %in% set_baseline)
   }

  if (length(to_test) == 0){
    warning("Baseline values for parameters you would like to test were not simulated!")
    #isolating parameter combinations of interest from the dataset
  }
  if(length(test_parameters) >= 1 && is.na(test_parameters) == FALSE){

    subset(parameters, parameters %in% test_parameters)
    subset(to_test, to_test %in% test_parameters)

  }else{
    to_test <- to_test
  }

  print("test parameters subsetted")

  if (length(to_test) == 0){
    warning("Parameter combinations you would like to test were not simulated!")
  }

  eq_infect <- c("N00", "A", "B", "A_plus", "B_plus", "coinf")

  #proportion of all coinfected

  coinf <- c(rowSums(to_test[, grep("^[^0slnbtKm]*$", colnames(to_test))]))

  #prop single B
    B <- apply(to_test,1, function(x)  x["N01"] / (x["K"] * x["mu"]))

  # double B

  #string pattern needed to recognise columns for species B co-infection

    patB <- rep(NA, model_det$endo_no_per_sp)

    for (i in 2:model_det$endo_no_per_sp){
      patB[i] <- c(paste0("N0", i))
    }

    patB <- na.omit(patB)
    dubB <- to_test |> dplyr::select(patB)
    dubB <- rowSums(dubB)
    dubBparam <- cbind(parameters, dubB)
    B_plus <- apply(dubBparam, 1, function(x) x["dubB"] / (x["K"] * x["mu"]))


     #prop single A
    A <- apply(to_test,1, function(x)  x["N10"] / (x["K"] * x["mu"]))

    # single species co-infections

    #string pattern needed to recognise columns for species A co-infections
    patA <- rep(NA, model_det$endo_no_per_sp)

    for (i in 2:model_det$endo_no_per_sp){
    patA[i] <- c(paste0("N", i, "0"))
    }

    patA <- na.omit(patA)
    dubA <- to_test |> dplyr::select(patA)
    dubA <- rowSums(dubA)
    dubAparam <- cbind(parameters, dubA)
    A_plus <- apply(dubAparam, 1, function(x) x["dubA"] / (x["K"] * x["mu"]))

    #prop_uninfected

    N00 <- apply(to_test, 1, function(x)  x["N00"] / (x["K"] * x["mu"]))

    eq_dat <- cbind(
      parameters,
      N00,
      A,
      B,
      A_plus,
      B_plus,
      coinf
    )

 # change test parameters so it's subsetable

      par <- colnames(parameters)
      eq_met <-tidyr::pivot_longer(eq_dat, all_of(par), names_to = "parameter")

    eq_met <- tidyr::pivot_longer(eq_met,
                                  all_of(eq_infect),
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
                   A_plus,
                   B_plus){
      D <-coinf * N0 - sum(A, A_plus) * sum(B, B_plus)

      #D prime

      if(D > 0){
        Dmax <- min(sum(A, A_plus) * (1-sum(B, B_plus)),
                    (1- sum(A, A_plus)) * sum(B, B_plus))
        Dp <- D/Dmax

      }else if (D < 0){
        Dmin <-  max(-(sum(A, A_plus) * sum(B, B_plus)),
                     -(1-sum(A, A_plus)) * (1-sum(B, B_plus)))
        Dp <- D/Dmin
      }

      #r2
      r <- D^2/sum(A, A_plus)*(1-sum(A, A_plus)) *
        sum(B, B_plus) * (1-sum(B, B_plus))

      #phi coefficient
      ph <- coinf*N0-sum(A, A_plus) * sum(B, B_plus)/sqrt(sum(coinf, sum(A, A_plus))*
                                                          sum(sum(B, B_plus), N0)*
                                                          sum(sum(A, A_plus), N0)*
                                                          sum(coinf, sum(B, B_plus)))

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
        eq_ld <- data.frame("ld" = NA, "Dprime" = NA, "r2" = NA, "phi_coef" = NA)
        cols <- c(colnames(parameters))
        eq_p <-data.frame(matrix(ncol = length(parameters)))
        colnames(eq_p) <- cols
        for (i in 1:nrow(to_test)){
          eq_p[i, ] <- parameters[i,]
            eq_ld[i, ] <- ld(
              coinf = eq_dat$coinf[i],
              N0 = eq_dat$N00[i],
              A = eq_dat$A[i],
              B = eq_dat$B[i],
              A_plus = eq_dat$A_plus[i],
              B_plus = eq_dat$B_plus[i]
            )
        }
        eq_ld <- cbind(eq_p, eq_ld)
      }else{
      eq_param_cor <- eq_av_prop |> dplyr::group_by(parameter, infection_status) |> dplyr::summarise(cor = cor(value, proportion, method = cor_param_method))
    }

      #prep file for saving
      if (any(cor_param_method == "ld")){
        return(list(
          model_det = model_det,
          "equilibrum_stats" = list(
            "Results" = to_test,
            "Proportion" = eq_dat,
            "correlations" = eq_param_cor,
            "ld" = eq_ld
          )
        ))
      } else{
        return(list(
          "model_det" = model_det,
          "equilibrum_stats" = list(
            "Results" = to_test,
            "Proportion" = eq_dat ,
            "correlation" = eq_param_cor
          )
        ))
      }
    }
}
