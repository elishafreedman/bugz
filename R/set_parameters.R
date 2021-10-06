#' set_parameters
#'
#' @param two_species : If TRUE,  parameters for two sets of species will be created.
#' @param dems :
#' @param K : carrying capacity of number of host species within the system.
#' @param lambda : host speciation rate
#' @param mu : loss rate  of host spcies
#' @param betaA :  transmission rate of species A
#' @param betaB : transmission rate of species B
#' @param sigmaA : priority effects of species A
#' @param sigmaB : priority effects of species B
#' @param sigmaAB : priority effects of co-infections
#' @param sigmaBA :priority effects of co-infections
#' @param nuA : loss rate of species A
#' @param nuB : loss rate of species B
#'
#' @return data frame of all parameter combinations
#' @export
#'
#' @examples set_parameters(two_species = TRUE, host_dem = TRUE, K = 200,lambda = 1,mu = 0.5,betaA = seq(0.0001,0.001, 0.0001),betaB = seq(0.0001, 0.001, 0.0001),sigmaA = seq(0, 1, 0.1),sigmaB = seq(0,1,0.1),sigmaAB = seq(0,1, 0.1),sigmaBA  = seq(0,1,0.1),nuA = 0.01,nuB = 0.01)

set_parameters <- function(two_species = TRUE, host_dem = TRUE,
                           K = 200,
                           lambda = 1,
                           mu = 0.5,
                           betaA = seq(0.0001,0.001, 0.0001),
                           betaB = seq(0.0001, 0.001, 0.0001),
                           sigmaA = seq(0, 1, 0.1),
                           sigmaB = seq(0,1,0.1),
                           sigmaAB = seq(0,1, 0.1),
                           sigmaBA  = seq(0,1,0.1),
                           nuA = 0.01,
                           nuB = 0.01){


  #creating the parameter list

  if(two_species == FALSE && host_dem == TRUE){
    sigmaB_range = 0
    nuB = 0
    beta_Brange = 0
    expand.grid(K= K,
                lambda = lambda,
                mu = mu,
                betaA= betaA,
                sigmaA= sigmaA,
                nuA = nuA)
    }else if(host_dem == FALSE && two_species == FALSE){
      expand.grid(K = K,
                  betaA= betaA,
                  sigmaA= sigmaA,
                  nuA = nuA)
    }


 if(host_dem == TRUE && two_species == TRUE){
    expand.grid(K= K,
                lambda = lambda,
                mu = mu,
                betaA= betaA,
                betaB = betaB,
                sigmaA= sigmaA,
                sigmaB = sigmaB,
                sigmaAB = sigmaAB,
                sigmaBA = sigmaBA,
                nuA = nuA,
                nuB = nuB)
 }else if(host_dem == FALSE && two_species == TRUE){
   expand.grid(K = K,
               betaA= betaA,
               betaB = betaB,
               sigmaA= sigmaA,
               sigmaB = sigmaB,
               sigmaAB = sigmaAB,
               sigmaBA = sigmaBA,
               nuA = nuA,
               nuB = nuB)
   }
}

