#test second model


  ODEsysten <- function(times, states, parameters){
    with(as.list(c(states, parameters)),{
      d0 <- lambda*N0*(1-(N0+A+B+AB)/K) - mu*N0                                               - f(0, c(A, B, AB), beta, sigma)*N0 + nu(1, nuB)*A
      dA <- lambda*A*(1-(N0+A+B+AB)/K) - mu*A + f(0, c(A, B, AB), beta, sigma)*N0 - f(1, c(A, B, AB), beta, sigma)*A + nu(2, nuB)*B - nu(1, nuB)*A
      dB <- lambda*B*(1-(N0+A+B+AB)/K) - mu*B + f(1, c(A, B, AB), beta, sigma)*A - f(2, c(A, B, AB), beta, sigma)*B + nu(3, nuB)*AB - nu(2, nuB)*B
      dAB<- lambda*AB*(1-(N0+A+B+AB)/K) - mu*AB + f(2, c(A, B, AB), beta, sigma)*B - f(3, c(A, B, AB), beta, sigma)*AB + nu(4, nuB)*N4 - nu(3, nuB)*AB
      return(list(c(d0, dA, dB, dAB, )))
    })
  }


parameters <-c(K = 200, beta = 0.0001,  sigma = 0.1, lambda = 1, mu = 0.5, nuB = 0.01, rhoA = 1, rhoB = 1, thetaA =, thetaB = )
iniState<-c(0 = 100, A = 1, B= 1, AB = 0)

out1 <- as.data.frame(ode(iniState, times, ODEsysten, parameters))

#plot

