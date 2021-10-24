#' run_model
#'
#' @param endo_species : number of endosymbiont species within the system
#' @param endo_number : number of endosymbionts per species within the system
#' @param parameters : a data frame containing all parameter combinations to run the model on
#' @param tmax : The maximum number of time steps
#' @param core_spec : If the number of parameter combinations is large it may be wise to assign multiple cores for speed. if NA, number of cores used in process will be set at half the number of cores available for use on the computer.
#'@param host_dem: if TRUE, host demographics (speciation, and extinction) will be included in the simulation.
#' @return list containing the details of the model,
#' a data frame of all parameter combinations,
#' and a list of model results for each parameter combination.
#' @export
#'
#' @examples
#' params <- set_parameters(two_species = TRUE,K = 200,lambda = 1,mu = 0.5,betaA =0.001,betaB = 0.001,sigmaA = 0.1,sigmaB = 0.1,sigmaAB = 1,sigmaBA = seq(0, 1, 0.1),nuA = 0.01,nuB = 0.01)
#' run_model(endo_species = 2,endo_number = 2,parameters = params,tmax = 1000,core_spec = NA, outfile = "ODE_results.rda", host_dem = FALSE)



run_model <- function(endo_species = 2,
                      endo_number = 2,
                      parameters = params,
                      tmax = 1000,
                      core_spec = NA, host_dem = TRUE){
  #Build the equations

  ODE <- build_equations(endo_s = endo_species, endo_no = endo_number, host = host_dem)

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
  results <- list()

  # create initial states vector #
 if(host_dem == TRUE){
  ini_state <- c(rep(0, length(ins)))
  names(ini_state) <- c(ins)
  ini_state <- dplyr::case_when(
    names(ini_state) == "N0"~ parameters$K[1] * parameters$mu[1],
    names(ini_state) == "N00"~parameters$K[1] * parameters$mu[1],
    names(ini_state) == "N1"~1,
    names(ini_state) == "N01"~1,
    names(ini_state) == "N10"~1
  )

  ini_state[is.na(ini_state)] <-  0
  names(ini_state) <- c(ins)
}else{
  ini_state <- c(rep(0, length(ins)))
  names(ini_state) <- c(ins)
  ini_state <- dplyr::case_when(
    names(ini_state) == "N0"~parameters$K[1],
    names(ini_state) == "N00"~parameters$K[1],
    names(ini_state) == "N1"~ 1,
    names(ini_state) == "N01"~1,
    names(ini_state) == "N10"~1
  )
  ini_state[is.na(ini_state)] <-  0
  names(ini_state) <- c(ins)
}

  print(paste("simulation start time", Sys.time()))


  ode_calc <- function(x) {
    list(Parameters = data.frame(t(x)), Results = data.frame(deSolve::ode(ini_state, times, eqn, x)))
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
      "ncore"
    ),
    envir = environment()
  )
  results <- pbapply::pbapply(parameters, 1, ode_calc, cl = clust)
  parallel::stopCluster(clust)

  #remove rows in
  results  <- lapply(results$simulation, function(x) x$Results[x$Results, all(cols(x$Results) !=0)])

  return(
    list(
      simulation_details = sim_details,
      param_combos = parameters,
      simulations = results
    )
  )
  print(paste("simulatione end time", Sys.time()))
}
