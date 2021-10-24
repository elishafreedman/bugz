#' int_plot
#'
#' @param parameters : parameter combinations that will be simulated
#' @param endo_species : number of endosymbiont species no
#' @param endo_number : number of individuals per species
#' @param tmax : max number of timesteps
#' @param host_dem : if TRUE, host demographics (speciation, and extinction) will be included in the simulation.
#' @param Colors :  vector of colours to be included in the graph
#' @return A shiny application line graph of number of hosts with intracellular endosymbionts within the system at each timestep at each parameter combination.
#' @export
#' @import shiny
#' @import ggplot2
#' @examples params <- set_parameters(two_species = TRUE,K = 200,lambda = 1,mu = 0.5,betaA =0.001,betaB = 0.001,sigmaA = 0.1,sigmaB = 0.1,sigmaAB = 1,sigmaBA = seq(0, 1, 0.1),nuA = 0.01,nuB = 0.01)
#' int_plot(parameters = params, endo_species = 2, endo_number = 1)
int_plot <-function(parameters = params, endo_species = 2, endo_number = 1, tmax = 2500, host_dem = TRUE, Colors = c("blue","orange", "red")){
  eqn <- build_equations(endo_s = endo_species, endo_no = endo_number, host = host_dem)
  ins <- eqn[["states"]]
  eqn_s <- eqn[["equations"]]
  times <- seq(0,tmax, 1)
  mins <- head(parameters, 1)
  maxs <- tail(parameters, 1)

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


  ini_state[is.na(ini_state)] <-  0
  names(ini_state) <- c(ins)
  plots <- function(inistat = ini_state,
                    time = times,
                    equation = eqn_s,
                    K = mins$K,
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
                    tmax = tmax,
                    endo_sp = endo_species,
                    endo_num = endo_number
                    ){
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
    initial_plot <- data.frame(deSolve::ode(inistat, time, equation, int_param))
# print(ncol(initial_plot))


    if(endo_sp == 2){
      # k <- (endo_sp+1)^(endo_num)-1
     cols <- grDevices::colorRampPalette(Colors)(ncol(initial_plot)-2)
    }else{
      # k <- endo_num-1
      cols <-  grDevices::colorRampPalette(Colors)(ncol(initial_plot)-2)
    }

    #reformat
    col <- c(colnames(initial_plot[,-1]))
    initial_plot <- tidyr::pivot_longer(initial_plot, all_of(col), names_to = "infection_status")
    initial_plot$infection_status <- factor(initial_plot$infection_status, levels = col)
    ggplot2::ggplot(initial_plot)+
      geom_line(aes(x=time, y = value,
                    colour = infection_status,
                    group = infection_status), size = 1.5)+
      ylab(label="# Host species")+
      xlab(label="Timesteps")+
      labs(colour = "infection status \n         A  B")+
      scale_colour_manual(values=c("black", cols))+
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
    }
    ui<- fluidPage(
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

