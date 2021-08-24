#N = the total number of host individuals

#K = carrying capacity

#k = the number of endosymbionts

#beta = the baseline transmission rate

#sigma = how strongly transmission rates decline based on which type of
#endosymbiont there is


# Function creating the ODE equations

#single transmission rate model added here represented by the argument "endosymbiont type"

build_equations <- function(endo_s = endo_species, endo_no = endo_number){
  k <- endo_no+1 #number of endosymbionts (including the uninfected state)


  #max types = 2 at the moment
  if(endo_s == 1){
    Ns <- rep("", k)
    for(i in 0:(k-1)){
      Ns[i + 1] <- paste0("N", i)
    }
    NN <- paste(Ns, collapse = "+")

    # left sides of equations:
    leftSides <- paste0("d", Ns)

    # matrix of terms for demography:
    matD <- matrix(rep(NA, k), nrow = 1)
    for(i in 0:(k-1)) {
      matD[i+1] <- paste0("lambda*N", i, "*(1 - (", NN, ")/K) - mu*N",i)
    }
    dim(matD) <- NULL

    # matrix of terms for the  influx force of infection

    # matrix for transmission events:


    N_A <- paste(Ns[-(0:(k-1)*k + 1)], collapse = "+")


    # total number of species with symbiont A:

    matrix0123 <- c(0:(k-1), k)

    pA <- paste0("sigmaA^", matrix0123)

                 #matrix of transmission
                 matT <- rep("", k)
                 for(i in 0:(k-1)){
                   if(i>0){ # influx of A symbionts
                     matT[i+1] <- paste0("+betaA*(", N_A,")*", pA[i], "*N", i-1)
                   }

                   if(i<(k-1)) {#outflux of A endosymbionts
                     matT[i+1] <- paste0(matT[i+1], "-betaA*(", N_A,")*", pA[i+1], "*N", i)
                   }
                 }




    matL <- matrix(rep(NA, k), nrow = 1)
    for(i in 0:(k-1)){
      matL[i+1] <- paste0("(-nuA)*N", i)
    }

    # matrix of terms for influx due to loss of symbionts:
    matG <- matrix(rep("", k), nrow = 1)
    for(i in 0:(k-1)) {
      if (i<(k-1))
        matG[i+1] <- paste0("nuA*N", i+1)
      else
        matG[i+1] <- "0"
    }
  }





   ###### two  types ########
  if(endo_s == 2){

    # string for total population size:
    Ns <- rep("", k^2)
    for(i in 0:(k-1)){
      for(j in 0:(k-1)){
        Ns[k*j + i + 1] <- paste0("N", i, j)
      }
    }
    NN <- paste(Ns, collapse = "+")
    NNb <- paste(Ns[-1], collapse = "+")

    # left sides of equations:
    leftSides <- paste0("d", Ns)

    # matrix of terms for demography:
    matD <- matrix(rep(NA, k^2), nrow = k)
    for(i in 0:(k-1)) {
      for(j in 0:(k-1)) {
        matD[i+1, j+1] <- paste0("lambda*N", i, j, "*(1 - (", NN, ")/K) - mu*N", i, j)
      }
    }
    dim(matD) <- NULL
    # matrix of terms for the  influx force of infection

    # matrix for transmission events:

    matT <- matrix(rep(NA, k^2), nrow = k)

    # total number of species with symbiont A:
    N_A <- paste(Ns[-(0:(k-1)*k + 1)], collapse = "+")
    # total number of species with symbiont B:
    N_B <- paste(Ns[-(1:k)], collapse = "+")

    # calculate pA an pB matrices (inhibition of transmission caused by pre-existing symbionts)
    matrix0123 <- matrix(rep(0:(k-1), k), nrow = k)

    pA <- matrix(paste0("sigmaA^", matrix0123, "*sigmaBA^", t(matrix0123)), nrow = k)
    pB <- matrix(paste0("sigmaB^", t(matrix0123), "*sigmaAB^", matrix0123),  nrow = k)

    #matrix of transmission
    matT <- matrix(rep("", (k)^2), nrow = k, byrow = TRUE)
    for(i in 0:(k-1)){
      for(j in 0:(k-1)){
        if(i>0){ # influx of A symbionts
          matT[i+1, j+1] <- paste0("+betaA*(", N_A,")*", pA[i, j+1], "*N", i-1, j)
        }
        if (j>0){# influx of B symbionts
          matT[i+1, j+1] <- paste0(matT[i+1, j+1], "+betaB*(", N_B,")*", pB[i+1, j], "*N", i, j-1)
        }
        if(i<(k-1)) {#outflux of A endosymbionts
          matT[i+1, j+1] <- paste0(matT[i+1, j+1], "-betaA*(", N_A,")*", pA[i+1, j+1], "*N", i, j)
        }
        if(j<(k-1)){# outflux of B symbionts
          matT[i+1, j+1] <- paste0(matT[i+1, j+1], "-betaB*(", N_B,")*", pB[i+1, j+1], "*N", i, j)
        }
      }
    }
    # matrix of terms for outflux due to loss of symbionts:
    matL <- matrix(rep(NA, k^2), nrow = k)
    for(i in 0:(k-1)){
      for(j in 0:(k-1)) {
        matL[i+1, j+1] <- paste0("(-", i, "*nuA - ", j, "*nuB)*N", i, j)
      }
    }

    # matrix of terms for influx due to loss of symbionts:
    matG <- matrix(rep("", k^2), nrow = k)
    for(i in 0:(k-1)) {
      for(j in 0:(k-1)) {
        if (i<(k-1))
          matG[i+1, j+1] <- paste0(i+1, "*nuA*N", i+1, j)
        else
          matG[i+1, j+1] <- "0"
        if (j<(k-1))
          matG[i+1, j+1] <- paste0(matG[i+1, j+1], " + ", j+1, "*nuB*N", i, j+1)
      }
    }
  }
    equations <- paste0("    ", leftSides, " = ", matD, matT, " + ", matL, " + ", matG, "\n")
    equations <- paste0(equations, collapse = "\n")
    functionFront <- "ODEsystem <- function(times, states, parameters){\n  with(as.list(c(states, parameters)),{\n"
    functionEnd <- paste0("\n    return(list(c(", paste0(leftSides, collapse = ", "), ")))\n  })\n}")
    completeFunction <- paste0(functionFront, equations, functionEnd)
    cat(completeFunction)

    ODEsystem <- eval(parse(text = completeFunction))

  #create vector of initial states for the run_model function to use
    ODE <- list(endo_s = endo_s,
                endo_no_per_sp = endo_no,
                equations = ODEsystem,
                states = Ns)
    return(ODE)
}

