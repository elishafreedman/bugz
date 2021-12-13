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
    for(i in 1:length(set_baseline)){
      subset(parameters, parameters %in% set_baseline[i])
      subset(to_test, to_test %in% set_baseline[i])
      print(paste(i, "baseline finished"))
    }
  }
  else{
    to_test <- all_results
    print("Baseline subsetted")
  }

  if (length(to_test) == 0){
    warning("Baseline values for parameters you would like to test were not simulated!")
    #isolating parameter combinations of interest from the dataset
  }
  if(length(test_parameters) >= 1 && is.na(test_parameters) == FALSE){
    for(i in 1:length(test_parameters)){
      subset(parameters, parameters %in% test_parameters[i])
      subset(to_test, to_test %in% test_parameters[i])
      print(paste(i, "parameters finished"))
    }

  }else{
    to_test <- to_test
  }
  # return(params)
  # return(to_test)
  print("test parameters subsetted")

  if (length(to_test) == 0){
    warning("Parameter combinations you would like to test were not simulated!")
  }

  # stats at equilibrium #

  # raw results for each parameter at equilibrium
  # colength <- rep(NA, ncol(to_test) + ncol(parameters))
  # to_test <- as.data.frame(matrix(NA, nrow = length(to_test),
  #                                ncol = length(colength)))
  # colnames(to_test) <- c(colnames(parameters),
  #                       colnames(to_test))

  #fill in parameters and and  infection status
  #if (model_det$endo_species >= 2){
  eq_infect <- c("N00", "A", "B", "A_plus", "B_plus", "coinf")
  #} else{
  #eq_infect <- c("N0", "A", "A_plus")
  #}



  # function to calculate the proportion at equilibrium



  #proportion co-infected
  #if (model_det$endo_species >= 2){




    prop <- function(x){
      x / (x$K * x$mu)
    }
    coinf <- c(rowSums(to_test[, grep("^[^0slnbtKm]*$", colnames(to_test))]))
    #create vector that will be used later on for the pivot longer
    #prop single B
    B <- apply(to_test,1, function(x)  x["N01"] / (x["K"] * x["mu"]))
    # double B
    #string pattern needed to recognise columns for species B co-infection
    #if (model_det$endo_species >= 2){
    patB <- rep(NA, model_det$endo_no_per_sp)
    for (i in 2:model_det$endo_no_per_sp){
      patB[i] <- c(paste0("N0", i))
    }
    patB <- na.omit(patB)
    dubB <- rowSums(to_test[to_test %in% patB])
    dubBparam <- cbind(parameters, dubB)
    B_plus <- apply(dubBparam, 1, function(x) x["dubB"] / (x["K"] * x["mu"]))
    #prop single A
    A <- apply(to_test,1, function(x)  x["N10"] / (x["K"] * x["mu"]))

    # single species co-infections

    #string pattern needed to recognise columns for species A co-infections
    patA <- rep(NA, model_det$endo_no_per_sp)
    # if (model_det$endo_species >= 2) {
    #   for (i in 2:model_det$endo_no_per_sp) {
    #     patA[i] <- c(paste0("N", i, "0"))
    #   }
    # } else{
    for (i in 2:model_det$endo_no_per_sp)
      patA[i] <- c(paste0("N", i))
    #}
    patA <- na.omit(patA)
    dubA <- rowSums(to_test[to_test %in% patA])
    dubAparam <- cbind(parameters, dubA)
    A_plus <- apply(dubAparam, 1, function(x) x["dubA"] / (x["K"] * x["mu"]))

    #prop_uninfected

    #if (model_det$endo_species >= 2){
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
 print(eq_dat)

 # change test parameters so it's subsetable


    #if (all(is.na(test_parameters)) == TRUE){
      par <- colnames(parameters)
      eq_met <-tidyr::pivot_longer(eq_dat, all_of(par), names_to = "parameter")

    # } else{
    #   eq_met <-tidyr::pivot_longer(eq_dat, all_of(test_parameters), names_to = "parameter")
    # }
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
        cols <- c(colnames(parameters))
        eq_p <-data.frame(matrix(ncol = length(parameters)))
        colnames(eq_p) <- cols
        for (i in 1:length(to_test)) {
          eq_p[i, ] <- parameters
          #if (model_det$endo_species >= 2) {
            eq_ld[i, ] <- ld(
              coinf = eq_dat$coinf[i],
              N0 = eq_dat$N00[i],
              A = eq_dat$A[i],
              B = eq_dat$B[i],
              A_plus = eq_dat$A_plus[i],
              B_plus = eq_dat$B_plus[i]
            )

          # } else{
          #   eq_ld[i, ]<- ld(
          #     coinf = eq_dat$coinf[i],
          #     N0 = eq_dat$N0[i],
          #     A = eq_dat$A[i],
          #     B = eq_dat$B[i],
          #     A_plus = eq_dat$A_plus[i],
          #     B_plus = eq_dat$B_plus[i]
          #   )
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
            "Results" = to_test,
            "Proportion" = eq_av_prop ,
            "correlations" = eq_param_cor,
            "ld" = eq_ld
          )
        ))
      } else{
        return(list(
          "model_det" = model_det,
          "equilibrum_stats" = list(
            "Results" = to_test,
            "Proportion" = eq_av_prop ,
            "correlation" = eq_param_cor
          )
        ))
      }
    }

