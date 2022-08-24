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
  model_det <- results_file[["simulation_details"]]
  parameters <- results_file[["param_combos"]]
  to_test <- results_file[["simulations"]]

  #subset  the baseline parameters


  if (length(set_baseline >= 1) && is.na(set_baseline) == FALSE){

    subset(parameters, parameters %in% set_baseline)
    subset(to_test, to_test %in% set_baseline)
  }

  if (length(to_test) == 0){
    warning("Baseline values for parameters you would like to test were not simulated!")
    #isolating parameter combinations of interest from the data set
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

  #total proportions

  test_sans_param <- to_test[, grep("^[^slnbtKm]*$", colnames(to_test))]

  prop_data <- data.frame(rowSums(test_sans_param))
  for(i in 1:length(prop_data)){
    fin_data <- apply(test_sans_param, c(1,2), function(x) x/prop_data[i,])
  }
  fin_data <- as.data.frame(fin_data)
  # proportion of all coinfected
  coinf <- c(rowSums(fin_data[, grep("^[^0slnbtKm]*$", colnames(fin_data))]))

  #proportion uninfected
  N00 <- fin_data[,"N00"]

  # proportion A
  A <- fin_data[,"N10"]
  #print(A)

  # proportion B
  B <- fin_data[,"N01"]


  # proportion double B

  #string pattern needed to recognise columns for species B co-infection

  patB <- rep(NA, model_det$endo_no_per_sp)

  for (i in 2:model_det$endo_no_per_sp){
    patB[i] <- c(paste0("N0", i))
  }

  patB <- na.omit(patB)
  dubB <- fin_data |> dplyr::select(patB)
  B_plus <- data.frame(B_plus = rowSums(dubB))
  # dubBparam <- cbind(parameters, dubB)


  # proportion double A

  #string pattern needed to recognise columns for species A co-infections
  patA <- rep(NA, model_det$endo_no_per_sp)

  for (i in 2:model_det$endo_no_per_sp){
    patA[i] <- c(paste0("N", i, "0"))
  }

  patA <- na.omit(patA)
  dubA <- fin_data |> dplyr::select(patA)
  A_plus <-  data.frame(A_plus = rowSums(dubA))
  print(class(A_plus))
  # print(A_plus)
  # dubAparam <- cbind(parameters, dubA)


  #proportions+parameters

  eq_infect <- c("N00", "A", "B", "A_plus", "B_plus", "coinf")

  eq_dat <- cbind(
    parameters,
    N00,
    A,
    B,
    A_plus,
    B_plus,
    coinf
  )
  print(eq_dat)

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
    D <-coinf*N0 - sum(A, A_plus) * sum(B, B_plus)

    #D prime

if(D > 0){
  Dmin <- min(sum(A, A_plus, coinf) * (1-sum(B, B_plus, coinf)),
              (1- sum(A, A_plus, coinf)) * sum(B, B_plus, coinf))
  Dp <- D/Dmin

}else{
  Dmax <-  max(-(sum(A, A_plus, coinf) * sum(B, B_plus, coinf)),
               -(1-sum(A, A_plus, coinf)) * (1-sum(B, B_plus, coinf)))
  Dp <- D/Dmax
}


    #r2
    r <- D^2/sum(A, A_plus)*(1-sum(A, A_plus)) *
      sum(B, B_plus) * (1-sum(B, B_plus))

    #phi coefficient
    ph <- coinf*N0-sum(A, A_plus) * sum(B, B_plus, coinf)/sqrt(sum(coinf, sum(A, A_plus))*
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
      for (i in 1:nrow(eq_dat)){
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

