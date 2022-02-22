#' run_model
#'
#' @param endo_species : number of endosymbiont species within the system
#' @param endo_number : number of endosymbionts per species within the system
#' @param parameters : a data frame containing all parameter combinations to run the model on
#' @param tmax : The maximum number of time steps
#' @param core_spec : If the number of parameter combinations is large it may be wise to assign multiple cores for speed. if NA, number of cores used in process will be set at half the number of cores available for use on the computer.
#' @param eq_threshold : the threshold of variance between different work
#' @return list containing the details of the model,
#' a data frame of all parameter combinations,
#' and a list of model results for each parameter combination.
#' @export
#'
#' @examples
#' params <- set_parameters(two_species = TRUE,K = 200,lambda = 1,mu = 0.5,betaA =0.001,betaB = 0.001,sigmaA = 0.1,sigmaB = 0.1,sigmaAB = 1,sigmaBA = seq(0, 1, 0.1),nuA = 0.01,nuB = 0.01)
#' run_model(endo_species = 2,endo_number = 2,parameters = params,tmax = 1000,core_spec = NA, outfile = "ODE_results.rda")



run_model <- function(endo_species = 1, endo_number = 2,
                      parameters = params,
                      tmax = 1000,
                      eq_threshold = 0.1,
                      core_spec = NA,
                      kmax = NA) {
  #Build the equations

  ODE <- build_equations(endo_no = endo_number, endo_s = endo_species)

  ## importing the model function details ##

  endo_species <- ODE[["endo_s"]]
  endo_number <- ODE[["endo_no_per_sp"]]
  eqn <- ODE[["equations"]]
  ins <- ODE[["states"]]



  # collate parameters, initial states, and equation to run model


  sim_details <- list(
    simulation_ran = Sys.time(),
    endo_species = endo_species,
    endo_no_per_sp = endo_number,
    max_timesteps = tmax
  )

  times <- seq(0, tmax, 1)

  if (is.na(kmax)) {
    kmax <- parameters$K[1] * parameters$mu[1]
  }

  # create initial states vector #

  ini_state <- c(rep(0, length(ins)))
  names(ini_state) <- c(ins)
  ini_state <- dplyr::case_when(
    names(ini_state) == "N00" ~  kmax,
    names(ini_state) == "N01"  ~  1,
    names(ini_state) == "N10"  ~ 1
  )

  ini_state[is.na(ini_state)] <-  0
  names(ini_state) <- c(ins)
  print(ini_state)
  print(paste("simulation start time", Sys.time()))



  isEqui <- function(x, eq_t = eq_threshold){
    x1 <- round(x[1], 2)
    x2 <- round(x[2], 2)
    max(abs(x1 - x2) / x1, na.rm = TRUE) <= eq_t

  }

  ode_calc <- function(x){
    Res <- data.frame(deSolve::ode(ini_state, times, eqn, x))
    for (k in 1:(nrow(Res) - 1)){
      eq <- as.matrix(Res[k:(k + 1), ])
      VarC <- apply(eq[, -1], 2, FUN = isEqui)
      if (all(VarC) == TRUE){
        eq_raw <- data.frame(c(x, Res[k + 1,]))
        break
      }
    }
    return(eq_raw)
  }


  if (is.na(core_spec) == TRUE) {
    ncore <- parallel::detectCores() / 2
  } else{
    ncore <- core_spec
  }
  print(paste("simulation using", ncore, "cores"))


  clust <- parallel::makeCluster(ncore, "PSOCK")
  parallel::clusterExport(
    cl = clust,
    varlist = c(
      "ode_calc",
      "ini_state",
      "times",
      "eqn",
      "parameters",
      "ncore",
      "eq_threshold",
      "isEqui"
    ),
    envir = environment()
  )

  sims <- pbapply::pbapply(parameters, 1, ode_calc, cl = clust)
  sims <- plyr::ldply(sims, rbind, .id = NULL)

  parallel::stopCluster(clust)


  return(
    list(
      "simulation_details" = sim_details,
      "param_combos" = parameters,
      "simulations" = sims
    )
  )
  print(paste("simulatione end time", Sys.time()))
}
