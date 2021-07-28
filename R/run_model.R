#' run_model
#'
#' @param endo_species : number of endosymbiont species within the system
#' @param endo_number : number of endosymbionts per species within the system
#' @param parameters : a data frame containing all parameter combinations to run the model on
#' @param tmax : The maximum number of time steps
#' @param int_plot : if TRUE launches a shiny app for quick glances of results
#' @param save_results : if TRUE saves the results as a list with the selected file name.
#'
#' @return Shiny app of all results and /or  a list containing the details of the model,
#' a data frame of all parameter combinations,
#' and a list of model results for each parameter combination.
#' @export
#'
#' @examples
run_model <- function(endo_species = 2,
                      endo_number = 2,
                      parameters = params,
                      tmax = 1000,
                      int_plot = TRUE, save_results = c(TRUE,  outfile = "ODE_results.rda")){


  #Build the equations

  ODE <- build_equations(endo_s = endo_species, endo_no = endo_number)

  ## importing the model function details ##

  endo_species <- ODE[["endo_s"]]
  endo_number <- ODE[["endo_no_per_sp"]]
  eqn <- ODE[["equations"]]
  ins <- ODE[["states"]]



  # collate parameters, initial states, and equation to run model


  sim_details <- list(simulation_ran = Sys.time(),
                      endo_species = endo_species,
                      endo_no_per_sp = endo_number,
                      max_timesteps = tmax)

  times <- seq(0, tmax, 1)
  results <- list()

  # create initial states vector #

ini_state <- c(rep(0, length(ins)))
names(ini_state) <- c(ins)
    ini_state <-  case_when(names(ini_state) == "N0"  ~  parameters$K[1]*parameters$mu[1],
              names(ini_state) == "N00" ~  parameters$K[1]*parameters$mu[1],
              names(ini_state) == "N1"   ~  1,
              names(ini_state) == "N01"  ~  1,
              names(ini_state) == "N10"  ~ 1)

ini_state[is.na(ini_state)] <-  0
names(ini_state) <- c(ins)
 d <- ini_state


#lapply
 saveR <- function(save_r  = save_results,
                   param = parameters,
                   inistate = ini_state,
                   equation = eqn,
                   time = times){
 print(paste("simulation start time", Sys.time()))

   if(any(save_r == TRUE)){
     ode_calc <- function(x){
       list(Parameters = data.frame(x), Results = data.frame(ode(inistate, time, equation,x)))
     }

     ncore <- detectCores()-1
     clust <- makeCluster(ncore, "PSOCK")
    clusterExport(cl = clust, varlist = c("ode_calc", "ode", "save_r", "inistate",
                  "time", "equation", "param"),
                  envir = environment())
    results <- pbapply(param, 1, ode_calc, cl = clust)
    stopCluster(clust)
   }
   ODE_results <- list(simulation_details = sim_details,
                       param_combos = parameters,
                       simulations = results)
   save(ODE_results, file = save_results["outfile"])

   load(save_results["outfile"], .GlobalEnv)
   print(paste("simulatione end time", Sys.time()))

 }

  if(int_plot == FALSE){
   P <-  saveR()
   return(P)
  }



  #running shiny app plot
  if(int_plot == TRUE){
    mins <- head(parameters, 1)
    maxs <- tail(parameters, 1)
    if(endo_species == 2){
      k <- (endo_species+1)^(endo_number)-1
      cols <-  colorRampPalette(c("royalblue3","orange", "red"))(k)
    }else{
      k <- (endo_species+1)^(endo_number)-1
      cols <-  colorRampPalette(c("royalblue3","blue"))(k)
    }
    plots <- function(inistat = ini_state, time = times,
                      equation = eqn, K = mins$K,
                      lambda = mins$lambda,
                      betaA = mins$betaA,
                      betaB = mins$betaB,
                      sigmaA = mins$sigmaA,
                      sigmaB = mins$sigmaB,
                      sigmaAB = mins$sigmaAB,
                      sigmaBA = mins$sigmaBA,
                      mu = mins$mu,
                      nuA = mins$nuA,
                      nuB = mins$nuB,
                      tmax = tmax){
      if(endo_species == 2){
        int_param <- c(K=K,
                       lambda = lambda,
                       betaA = betaA,
                       betaB = betaB,
                       sigmaA = sigmaA,
                       sigmaB = sigmaB,
                       sigmaAB = sigmaAB,
                       sigmaBA = sigmaBA,
                       mu = mu,
                       nuA = nuA,
                       nuB = nuB)
      }else{
        sigmaB = NA
        sigmaAB = NA
        sigmaBA = NA
        nuB = NA
        betaB = NA
        int_param <- c(K=K,
                       lambda = lambda,
                       betaA = betaA,
                       sigmaA = sigmaA,
                       mu = mu,
                       nuA = nuA)
      }
      #run model
      initial_plot <- data.frame(ode(inistat, time, equation, int_param))


      #reformat
      col <- c(colnames(initial_plot[,-1]))
      initial_plot <- pivot_longer(initial_plot, all_of(col), names_to = "infection_status")
      initial_plot$infection_status <- factor(initial_plot$infection_status, levels = col)

      ggplot(initial_plot)+
        geom_line(aes(x=time, y = value,
                      colour = infection_status,
                      group = infection_status,
                      linetype = infection_status), size = 1.5)+
        ylab(label="# Host species")+
        xlab(label="Timesteps")+
        labs(colour = "infection status \n         A  B",
             linetype = "infection status \n         A  B")+
        scale_colour_manual(values=c("black",  cols))+
        theme_classic()+
        theme(axis.text = element_text(size = 20))+
        theme(axis.title = element_text(size =20))+
        theme(legend.text = element_text(size = 20))+
        theme(legend.key.size = unit(3, "lines"))+
        theme(legend.title = element_text(size = 20))
    }

    if(endo_species >= 2){
      server<-function(input, output){
        output$Main_plot <- renderPlot({plots(betaA = input$betaA,
                                              betaB = input$betaB,
                                              sigmaA = input$sigmaA,
                                              sigmaB = input$sigmaB,
                                              sigmaAB = input$sigmaAB,
                                              sigmaBA = input$sigmaBA,
                                              nuA = input$nuA,
                                              nuB = input$nuB,
                                              tmax = input$tmax)

        })
        onStop(function() P<-saveR())
      }
      ui<-fluidPage(
        titlePanel("Endosymbiont coinfections"),
        sidebarLayout(
          sidebarPanel(sliderInput(inputId = "betaA", label = "Transmission rate: \u03B2 A", min = mins$betaA, max = maxs$betaA, value = mins$betaA, step = 0.0001, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "betaB", label = "Transmission rate: \u03B2 B", min = mins$betaB, max = maxs$betaB, value = mins$betaB, step = 0.0001, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "sigmaA", label = "Decline in transmission \u03C3 A", min = mins$sigmaA, max = maxs$sigmaA, value = mins$sigmaA, step = 0.1, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "sigmaB", label = " Decline in transmission \u03C3 B", min = mins$betaB, max = maxs$sigmaB, value = mins$betaB, step = 0.1, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "sigmaAB", label = " Decline in transmission \u03C3 AB", min = mins$sigmaAB, max = maxs$sigmaAB, value = mins$sigmaAB, step = 0.1,  animate = animationOptions(interval = 1000, loop = FALSE)),
                       sliderInput(inputId = "sigmaBA", label = " Decline in transmission \u03C3 BA", min = mins$sigmaBA, max = maxs$sigmaBA, value = mins$sigmaBA, step = 0.1, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "nuA",label = " baseline loss rate of species A \u03BD A", min = mins$nuA, max = maxs$nuA, value = mins$nuA, step = 0.001,  animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "nuB",label = "baseline loss rate of species B \u03BD B", min = mins$nuB, max = maxs$nuB, value = mins$nuB, step = 0.001, animate = animationOptions(interval = 100, loop = FALSE)),
                       selectInput(inputId = "tmax", label= "Time to be simulated", choices = seq(0, tmax, tmax/4), selected=tmax), width = 3, length = 4),
          mainPanel(plotOutput(outputId = "Main_plot", height = "800px"))
        )
      )
    }else{
      server<-function(input, output, session){
        output$Main_plot <- renderPlot({plots(betaA = input$betaA,
                                              sigmaA = input$sigmaA,
                                              nuA = input$nuA,
                                              tmax = input$tmax)

        })
        onStop(function() P <- saveR())
      }
      ui<-fluidPage(
        titlePanel("Endosymbiont coinfections"),
        sidebarLayout(
          sidebarPanel(sliderInput(inputId = "betaA", label = "	Transmission rate: \u03B2 A", min = mins$betaA, max = maxs$betaA, value = mins$betaA, step = 0.0001, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "sigmaA", label = "Decline in transmission \u03C3 A", min = mins$sigmaA, max = maxs$sigmaA, value = mins$sigmaA, step = 0.1, animate = animationOptions(interval = 100, loop = FALSE)),
                       sliderInput(inputId = "nuA",label = " baseline loss rate of species A \u03BD A", min = mins$nuA, max = maxs$nuA, value = mins$nuA, step = 0.001, animate = animationOptions(interval = 100, loop = FALSE)),
                       selectInput(inputId = "tmax", label= "Time to be simulated", choices = seq(0, tmax, tmax/4), selected=tmax), width = 3, length = 4),
          mainPanel(plotOutput(outputId = "Main_plot", height = "800px"))
        )
      )
    }
    return(shinyApp(ui=ui, server=server))
  }
 }

