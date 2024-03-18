#install.packages(deSolve)
library(deSolve)
#install.packages(ggplot2)
library(ggplot2)
#install.packages(dplyr)
library(dplyr)
#install.packages(shiny)
library(shiny)

ebola.app <- shinyApp(
  # This creates the User Interface (UI)
  ui <- pageWithSidebar(
    headerPanel("A Two-Patch Stochastic Model for the Seasonal Spillover of Ebola Virus"),
    sidebarPanel(
      sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(0,1)
      ),
      sliderInput("ATR", "Amplitude Transmission Rate:",
                  min = 0, max = 0.5, value = 0.5
      ),
      sliderInput("ADR", "Amplitude Dispersal Rate:",
                  min = 0, max = 0.5, value = 0.5
      ),
      sliderInput("AV_TR_A", "Avg. Transmission Rate (Patch A):",
                  min = 0, max = 60, value = 60
      ),
      sliderInput("AV_TR_B", "Avg. Transmission Rate (Patch B):",
                  min = 0, max = 5, value = 5
      ),
      sliderInput("AV_DR_AB", "Avg. Dispersal Rate (AB):",
                  min = 0, max = 0.4, value = 0.4
      ),
      sliderInput("AV_DR_BA", "Avg. Dispersal Rate (BA):",
                  min = 0, max = 0.04, value = 0.04
      )
      
    ),
    mainPanel(
      tabsetPanel(
        
        tabPanel("Two-Patch Model",
                 #         withMathJax(
                 
                 #          helpText("This is a Two-Patch$\hat{x}$ ODE EBOLA Model"))
                 #),
                 tabPanel("", plotOutput("plot1"),
                          helpText("")),
                 tabPanel("", plotOutput("plot2"),
                          helpText(""))
                 
        ),
        
        tabPanel("Model Description",
                 withMathJax(
                   helpText("Here we provide an interective two-patch SIR stochastic epidemic model for the spillover of Ebola virus with non-seasonal, seasonal and demographic variability."),
                   helpText("The model is expressed as a system of coupled ordinary differential equations (ODE) with periodic disease transmission and dispersal between the PATCH A, an endemic enviroment for the virus, e.g. intact forest, and the PATCH B, an human modified enviromental, e.g. village or human settlement. Each population patch is partitioned into Susceptible (S), Infected (I), and Recovered (R) sub-populations."),
                   helpText("The total population size (N) in each patch, for j = A, B, is: $$N_j = S_j +I_j +R_j$$  For the current model: $$N_A = 100$$ $$N_B = 1000$$"),
                   helpText("We also defined natural birth and death rate at patch A and B:"),
                   helpText("$$ \\mu_A=0.06$$"),
                   helpText("$$ \\mu_B=0.02$$"),
                   helpText("And recovery rate for patch A and B, for j = A,B as:"),
                   helpText("$$\\gamma_j=12$$"),
                   helpText("The disease transmission rate and the disease dispersal rate are assumed to be time-periodic functions due to the seasonality of the diseases,where h = S,I,R represents the movement rate for susceptible, or infected, or recovered individuals. "),
                   helpText("Disease transmission rate: $$\\beta_{j}(t)=\\beta_{j}(t+P)$$"),
                   helpText("The disease transmission rate: $$\\phi_{jk}^{h}(t)=\\phi_{jk}^{h}(t+P)$$"),
                   helpText("We assume these functions are continuous for t ∈ (−∞, ∞) with period P > 0."),
                   helpText("The system of ODEs for the two-patch SIR model takes the following form:"),
                   helpText("Patch A:"),
                   helpText("$$ \\dot{S}_A = \\mu_A N_A - \\beta_A(t)\\frac{I_A}{N_A}S_A - (\\mu_A+ \\phi_{AB}^{S}(t))S_A + \\phi_{BA}^{S}(t) S_B$$"),
                   helpText("$$ \\dot{I}_A = \\beta_A(t)\\frac{I_A}{N_A}S_A -(\\mu_A + \\gamma_A + \\phi_{AB}^{I}(t))I_A + \\phi_{BA}^{I}(t) I_B  $$"),
                   helpText("$$ \\dot{R}_A = \\gamma_A I_A - (\\mu_A + \\phi_{AB}^{R}(t))R_A + \\phi_{BA}^{R}(t)R_B  $$"),
                   helpText("Patch B:"),
                   helpText("$$ \\dot{S}_B = \\mu_B N_B - \\beta_B(t)\\frac{I_B}{N_B}S_B - (\\mu_B+ \\phi_{BA}^{S}(t))S_B + \\phi_{AB}^{S}(t) S_A$$"),
                   helpText("$$ \\dot{I}_B = \\beta_B(t)\\frac{I_B}{N_B}S_B -(\\mu_B + \\gamma_B + \\phi_{BA}^{I}(t))I_B + \\phi_{AB}^{I}(t) I_A  $$"),
                   helpText("$$ \\dot{R}_B = \\gamma_B I_B - (\\mu_B + \\phi_{BA}^{R}(t))R_B + \\phi_{AB}^{R}(t)R_A  $$"),
                   helpText("We assume a simple sinusoidal functions form for the disease transmission rates and dispersal rates with period P for the two-patch system."),
                   helpText("Disease transmission rate: $$\\beta_{j}(t)=\\overline{\\beta_j} (1+σ_j cos(\\frac{2\\pi t}{P}))$$ "),
                   helpText("Disease dispersal rate:$$\\phi_{jk}^{h}(t)=\\overline{\\phi_{jk}^{h}} (1+σ_{jk}^{h} cos(\\frac{2\\pi t}{P}))$$"),
                   helpText("Our Two-Patch Model utilizes the following initial parameters values, which can be interactively changed to investigate the spillover of Ebola in different scenarios:"),
                   helpText("Amplitude transmission rate:"),
                   helpText("$$σ_j=0.5$$"),
                   helpText("Amplitude dispersal rate for h=S,I,R:"),
                   helpText("$$σ_{jk}^{h}=0.5$$"),
                   helpText("Average transmission rate at patch A:"),
                   helpText("$$\\overline{\\beta_A}=60$$"),
                   helpText("Average transmission rate at patch B:"),
                   helpText("$$\\overline{\\beta_B}=5$$"),
                   helpText("Average dispersal rate at patch AB:"),
                   helpText("$$\\overline{\\beta_j}=0.4$$"),
                   helpText("Average dispersal rate at patch BA:"),
                   helpText("$$\\overline{\\phi_{jk}^{h}}=0.04$$"),
                   helpText("The time range can be chosen between 0 and 100 years."),
                   helpText(" You are welcome to changed these parameters and investigate the spillover of Ebola in different scenarios!")
                 )
        ))
    )
  ),
  
  # This creates the 'behind the scenes' code (Server)
  server <- function(input, output) {
    
    textsize <- 14
    
    #        mytheme <- theme(legend.position = "top",
    #                        axis.title = element_text(size = textsize + 2),
    #                       axis.text = element_text(size = textsize),
    #                      legend.title = element_text(size = textsize + 2),
    #                     legend.text = element_text(size = textsize))
    
    two_patch_SIR <- function(time, variables, parameters, phi_BA, phi_AB, beta_A, beta_B) {
      with(as.list(c(variables, parameters)), {
        
        dS_A <- (mu_A * N_A) - ((beta_A(time))*(I_A/N_A)*(S_A)) - ((mu_A + phi_AB(time))*S_A) + (phi_BA(time)*S_B) # Correct.
        
        dI_A <- (beta_A(time)*(I_A/N_A)*S_A) - ((mu_A + gamma_A + phi_AB(time)) * I_A) + (phi_BA(time)*I_B) # Correct.
        
        dR_A <- (gamma_A*I_A)-((mu_A + phi_AB(time))*R_A) + (phi_BA(time)*R_B)  # Correct.
        
        #-------------------------------------------------------------------#
        dS_B <-(mu_B * N_B) - ((beta_B(time))*(I_B/N_B)*(S_B)) - ((mu_B + phi_BA(time))*S_B) + (phi_AB(time)*S_A) # Correct.
        
        dI_B <-(beta_B(time)*(I_B/N_B)*S_B) - ((mu_B + gamma_B + phi_BA(time))*I_B) + (phi_AB(time)*I_A) # Correct.
        
        dR_B <-(gamma_B*I_B)-((mu_B + phi_BA(time))*R_B) + (phi_AB(time)*R_A) # Correct.
        
        return(list(c(dS_A, dI_A, dR_A, dS_B, dI_B, dR_B)))
      })
    }
    
    output$plot1 <- renderPlot({
      
      #####################################################
      # Disease Transmission Rate by (t) with Seasonality #
      #####################################################
      
      beta_A<-function(time){
        return(AV_TR_A*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
      }
      
      beta_B<-function(time){
        
        return(AV_TR_B*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
      }
      ##########################
      # Disease Dispersal Rate #  
      ##########################
      
      phi_AB<-function(time){
        return(AV_DR_AB*(1 + (ADR*cos((2*pi*time/1))))) # Correct.
      }
      
      phi_BA<-function(time){
        
        return(AV_DR_BA*(1 + (ADR*(cos((2*pi*time)/1))))) # Correct.
      }
      
      # Set up initial values for the variables
      initial_values <- c(
        S_A = 95, I_A = 5, R_A = 0,   # Initial values for patch A
        S_B = 995, I_B = 5, R_B = 0   # Initial values for patch B
      )
      
      # Define parameters for the system
      parameters_values <- list(
        N_A = 100,
        N_B = 1000,
        mu_A = 0.06,
        mu_B = 0.02,
        gamma_A = 12,
        gamma_B = 12,
        phi_BA = phi_BA,
        phi_AB = phi_AB,
        beta_A = beta_A,
        beta_B = beta_B
        
      )
      
      #########################################
      #            Dynamical parameters       #
      #---------------------------------------#
      
      time <- seq(input$T[1], input$T[2], by = 1/100)
      ATR<-(input$ATR)
      ADR<-(input$ADR)
      AV_TR_A<-(input$AV_TR_A)
      AV_TR_B<-(input$AV_TR_B)
      AV_DR_AB<-(input$AV_DR_AB)
      AV_DR_BA<-(input$AV_DR_BA)
      
      #########################################
      #            TWO PATCH ODE              #
      #---------------------------------------#
      out <- data.frame(ode(
        y = initial_values,
        times = time,
        func = two_patch_SIR,
        parms = parameters_values
      ))
      
      out %>% ggplot() +
        geom_line(aes(x = time, y = S_A, colour = "aS"), size = 1, alpha = 0.7) +
        geom_line(aes(x = time, y = I_A, colour = "dR"), size = 1, alpha = 0.7) +
        geom_line(aes(x = time, y = R_A , colour = "bE"), size = 1, alpha = 0.7) +
        scale_x_continuous("Years") + scale_y_continuous("Population") +
        scale_colour_discrete("Individuals", labels = c("S", "I", "R")) + ggtitle("PATCH A") + theme(legend.position = "right")
      
    })
    
    output$plot2 <- renderPlot({
      
      #####################################################
      # Disease Transmission Rate by (t) with Seasonality #
      #####################################################
      seasonality=0.5  #
      ##################
      
      beta_A<-function(time){
        return(AV_TR_A*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
      }
      
      beta_B<-function(time){
        
        return(AV_TR_B*(1+(ATR*cos((2*pi*time)/1)))) # Correct.
      }
      ##########################
      # Disease Dispersal Rate #  
      ##########################
      
      phi_AB<-function(time){
        return(AV_DR_AB*(1 + (ADR*cos((2*pi*time/1))))) # Correct.
      }
      
      phi_BA<-function(time){
        
        return(AV_DR_BA*(1 + (ADR*(cos((2*pi*time)/1))))) # Correct.
      }
      
      # Set up initial values for the variables
      initial_values <- c(
        S_A = 95, I_A = 5, R_A = 0,   # Initial values for patch A
        S_B = 995, I_B = 5, R_B = 0   # Initial values for patch B
      )
      
      # Define parameters for the system
      parameters_values <- list(
        N_A = 100,
        N_B = 1000,
        mu_A = 0.06,
        mu_B = 0.02,
        gamma_A = 12,
        gamma_B = 12,
        phi_BA = phi_BA,
        phi_AB = phi_AB,
        beta_A = beta_A,
        beta_B = beta_B
      )
      
      #########################################
      #            Dynamical parameters       #
      #---------------------------------------#
      time <- seq(input$T[1], input$T[2], by = 1/100)
      ATR<-(input$ATR)
      ADR<-(input$ADR)
      AV_TR_A<-(input$AV_TR_A)
      AV_TR_B<-(input$AV_TR_B)
      AV_DR_AB<-(input$AV_DR_AB)
      AV_DR_BA<-(input$AV_DR_BA)
      
      #########################################
      #            TWO PATCH ODE              #
      #---------------------------------------#
      out <- data.frame(ode(
        y = initial_values,
        times = time,
        func = two_patch_SIR,
        parms = parameters_values
      ))
      
      out %>% ggplot() +
        geom_line(aes(x = time, y = S_B, colour = "aS"), size = 1, alpha = 0.7) +
        geom_line(aes(x = time, y = I_B, colour = "dR"), size = 1, alpha = 0.7) +
        geom_line(aes(x = time, y = R_B , colour = "bE"), size = 1, alpha = 0.7) +
        scale_x_continuous("Years")+scale_y_continuous("Population") + ggtitle("PATCH B")+
        scale_colour_discrete("Individuals", labels = c("S", "I", "R")) + theme(legend.position = "right")
      
    })
  }
)

