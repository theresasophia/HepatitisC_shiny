######################################################################
## "Dynamic modelling of hepatitis C transmission among injecting drug 
## users: A tool to support WHO elimination targets"
##
##
## Author: Theresa Stocks <http://www.su.se/english/profiles/tstoc-1.219526>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 13/09/2018
######################################################################



library(rsconnect)
shinyUI(bootstrapPage(
  titlePanel("Dynamic modelling of hepatitis C transmission among injecting drug users "),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      h4("Input model parameters"),
      br(),     br(),
      sliderInput(inputId = "acute",label = HTML("<p>Y<sub>A</sub> (diagnosed acute cases per year)") ,value=129, min =100, max = 150, step = 1),
      selectInput(inputId =  "acute.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(),
      sliderInput(inputId = "chronic",label =  HTML("<p>Y<sub>C</sub> (diagnosed chronic cases per year)"), value=1129, min =1000, max = 1250, step = 1),
      selectInput(inputId =  "chronic.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(),
      sliderInput(inputId = "immigrants",label =  HTML("<p>Y<sub>I</sub> (diagnosed imported chronic cases per year)"), value=184, min =150, max = 250, step = 1),
      selectInput(inputId =  "immigrants.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(), br(),
      sliderInput(inputId = "pop",label = "N (active IDU population size in endemic state)", value=26000, min =8000, max = 40000, step = 1000),
      selectInput(inputId =  "pop.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(),   br(),br(),
      sliderInput(inputId = "prev",label = "prev (HCV prevalence among active IDU before interventions (%)):", value=0.817, min =0.2, max = 0.99, step = 0.05),
      selectInput(inputId =  "prev.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(), br(),
      sliderInput(inputId = "cess",label = "c (permanent IDU cessation rate)", value=17.8, min =10, max = 20, step = 1),
      selectInput(inputId =  "cess.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(), br(),
      sliderInput(inputId = "death",label = HTML("&mu; (IDU mortality rate)"), value=50, min =40, max = 82, step = 1),
      selectInput(inputId =  "death.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(),     br(), br(),
      sliderInput(inputId = "deathchronic",label = HTML("&mu;+&rho;  (IDU with HCV mortality rate)"), value=50, min =10, max = 50, step = 1),  
      selectInput(inputId =  "deathchronic.uc", label = h5("Standard error"), 
                  choices = list("0 %" = 0,"1 %" = 0.01, "2 %" = 0.02,"3 %" = 0.03,"4 %" = 0.04,"5 %" = 0.05,"6 %" = 0.06,"7 %" = 0.07," 8 %" = 0.08,"9 %" = 0.09,"10 %" = 0.1, "20%" = 0.2,"30 %" = 0.3, "40%" = 0.4,"50 %" = 0.5, "60 %" = 0.6, "70%" = 0.7,"80 %"=0.8, "90 %"=0.9, "100 %"=1), 
                  selected = 0.01),
      br(), br(),
      br(), br(),
      radioButtons( "sample", label = h3("Sample size"), choices = list("10" = 10, "100" = 100, "250" = 250, "1000"=1000), selected = 100),
      radioButtons( "years_to_simulate", label = h3("Years to simulate"), choices = list("5" = 5, "10" = 10, "25" = 25, "50" = 50), selected = 10)),
  
    
 mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Information",  h3("Objective"),
                           p("This project is investigating the burden of undiagnosed hepatitis C (HCV)
                             cases, incidence, average time until diagnosis, reproduction number and the impact of treatment for countries where 
                             the main route of HCV transmission is via injecting drug users sharing contaminated needles in order to 
                             support the WHO goal of viral HCV elimination by 2030. This is a collaborative work between
                             Stockholm University and the Public Health Agency of Sweden and is based on the manuscript (link to follow)."),
                           h3("Abstract"),
                           p("To reach the WHO goal of hepatitis C elimination, it is essential to identify the  number of people 
                              unaware of their infection, the true incidence and to investigate the effect of interventions on the 
                              disease transmission dynamics. In the developed countries, one of the primary routes of hepatitis C 
                              virus (HCV)  transmission is via contaminated needles shared by injecting drug users (IDUs). However, 
                              substantial underreporting combined with  high uncertainty regarding the size of this difficult to reach 
                              population makes it challenging to estimate the core indicators recommended by the WHO. In order to 
                              enable countries to monitor their progress towards the elimination goal, we present a novel, deterministic 
                              multi-layered dynamic transmission model for HCV transmission within an IDU community. 
                              First, based on this model, and using routine surveillance data, we estimate the number of undiagnosed 
                              IDUs, the true incidence, the  reproduction numbers, the average time until diagnosis and associated 
                              uncertainties. Second, we examine the impact of two interventions on disease dynamics: 1) direct-acting
                              antiviral drug treatment, and 2) needle exchange programs. The novelty of our approach is the explicit 
                              distinction into diagnosed and undiagnosed cases and active and former IDUs which enables new insights 
                              on the transmission dynamics and intervention strategies. 
                              We illustrate our approach for Swedish data, however, our model can  be  easily applied to countries 
                              with a similar disease burden and risk groups. To make the model accessible for relevant users and to 
                              support communication of our results to public health decision makers, the model and its output are
                              made available through a shiny app."),
                           h3("Schematic model representation"),
                           br(),br(),
                           img(src = "new_fig.pdf", height = 500, width = 800),
                           p(HTML("Schematic representation of the base transmission model and the  needle exchange program (NEP). 
                              We distinguish between susceptible (S<sub>{d,e,f}</sub>) and two  stages  of infection (acute (A<sub>{d,e,f}</sub>) 
                              and chronic (C<sub>{d,e,f}</sub>)). The indices indicate whether or not (0=no and 1=yes) d = diagnosed, 
                              e = enrolled in needle exchange program and f = former IDU  and super-scripted with I if imported. 
                              The rates on the arrows are explained in the complementing manuscript."))
                           ),
                  
                  tabPanel("Parameter distributions",h3("Distributions"), br(),
                           plotOutput(outputId = "YA_dist", height = 200, width = 300),   
                           br(),     br(), br(),
                           plotOutput(outputId = "YC_dist", height = 200, width = 300),
                           br(),     br(), br(),
                           plotOutput(outputId = "YI_dist", height = 200, width = 300),
                           br(),     br(), br(),
                           plotOutput(outputId = "pop_dist", height = 200, width = 300), 
                           br(),     br(), 
                           plotOutput(outputId = "prev_dist", height = 200, width = 300), 
                           br(),     br(), 
                           plotOutput(outputId = "cess_dist", height = 200, width = 300),
                           br(),     br(), 
                           plotOutput(outputId = "death_dist", height = 200, width = 300), 
                           br(),     br(), 
                           plotOutput(outputId = "deathchronic_dist", height = 200, width = 300)), 
  
                  
                  
                  tabPanel("Outcomes base model",   h3("Key outcomes"),
                           p("Here we present the results assuming that the system is in equilibrium."),
                           h4("Estimated Reproduction Numbers"),
                           textOutput("R0_in"),
                           plotOutput(outputId = "r0_in", height = 200, width = 300),
                           textOutput("R0_out"),
                           plotOutput(outputId = "r0_out", height = 200, width = 300),
                           h4("Estimated true yearly incidence"),
                           textOutput("incidence_txt"),
                           plotOutput(outputId = "incidence_pl", height = 200, width = 300),
                           h4("Estimated time until diagnosis"),
                           textOutput("time_till_diagnosis_txt"),
                           plotOutput(outputId = "time_till_diagnosis_pl", height = 200, width = 300),
                           h4("Estimated number of  undiagnosed"),
                           p("The table shows the number of undiagnosed cases stratified by status of IDU and disease progression (CI)."),
                           tableOutput(outputId = "table_undia"),
                           textOutput("total_undiagnosed"),
                           h4("Sensitivity analysis"),
                           tableOutput(outputId = "table_sens")),


                            
                  tabPanel("Interventions",h3("Treatment and needle exchange programs"),
                           
                           sidebarLayout(
                             
                             # Sidebar panel for inputs ----
                             sidebarPanel(
                               h4("Input intervention parameters"),
                               br(),
                               h4("Treatment in base model"),
                               br(),
                               sliderInput(inputId = "trea",label =  HTML("<p>&tau;<sub>a</sub>  (treatment rate active IDUs)") ,value=0, min =0, max = 2, step = 0.1),
                               sliderInput(inputId = "trea_c",label = HTML("<p>&tau;<sub>f</sub>  (treatment rate former IDUs)") ,value=0, min =0, max = 2, step = 0.1),
                               h4("Treatment with NEP"),
                               br(),
                               sliderInput(inputId = "e",label = HTML("&epsilon;  (enrollment rate into NEPs)") ,value=0, min =0, max = 2, step = 0.1),
                               sliderInput(inputId = "d",label = HTML("&delta;  (drop-out rate NEPs)")  ,value=0, min =0, max = 2, step = 0.1),
                               sliderInput(inputId = "trea_ne",label = HTML("<p>&tau;<sub>e</sub>  (treatment rate of IDUs in NEP)") ,value=0, min  =0, max = 1, step = 0.1),
                               sliderInput(inputId = "rE",label =  HTML("<p>r<sub>e</sub>  (diagnosis rate in NEP)")  ,value=1/13, min =0, max = 5, step = 0.5),
                               selectInput(inputId =  "alpha", label = HTML("&alpha; (relative risk of acquiring/spreading HCV)"), 
                                           choices = list("1" = 1,"0.9" = 0.9, "0.8" = 0.8,"0.7" = 0.7, "0.6" = 0.6,"0.5" = 0.5, "0.4" = 0.4, "0.3" = 0.3,"0.2"=0.2, "0.1"=0.1, "0"=0), 
                                           selected = 1)),
                             
                           
                             mainPanel(
                               br(),
                               h4("Impact of direct-acting antiviral (DAA) drug treatment and needle exchange programs"),
                               plotOutput(outputId = "sim_plot"),                         
                               tableOutput(outputId = "table_trea"),
                               br(),br(),br(),br())
                           ))
                  )))
      
    ))
    
   
  
  


