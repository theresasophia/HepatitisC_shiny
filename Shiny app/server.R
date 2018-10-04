######################################################################
## R code generating the output for the shiny app for the manuscript
## "Dynamic modelling of hepatitis C transmission among injecting drug 
## users: A tool to support WHO elimination targets"
##
##
## Author: Theresa Stocks <http://www.su.se/english/profiles/tstoc-1.219526>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 13/09/2018
######################################################################

#load all reqired packages
library(shiny)
library(pomp)
library(magrittr)
library(plyr)
library(reshape2)
library(ggplot2)
require(nleqslv)
library(scales)
library(foreach)
library(truncnorm)
library(doParallel)
library(readr)
library(dplyr)
library(mice)
library(ggplot2)
require(mgcv)
library(tidyverse)
library(lubridate)
library(magrittr)
registerDoParallel()



cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)


shinyServer(function(input, output) {
   
  curEndemic  <- reactive({ endemic(input$pop, input$pop.uc, 
                                    input$death, input$death.uc,
                                    input$deathchronic, input$deathchronic.uc, 
                                    input$cess, input$cess.uc,
                                    input$prev, input$prev.uc,
                                    input$acute, input$acute.uc, 
                                    input$chronic, input$chronic.uc,
                                    input$immigrants, input$immigrants.uc, 
                                    input$sample,input$trea,input$trea_c,
                                    input$trea_ne,input$alpha,input$e,input$d,input$rE) })
  

  # Distibution of gamma
  output$g_dist <- renderPlot({
    gamma<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"g"]), aes(curEndemic( )[["box"]][,"g"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",g)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(gamma)
  })
  
  # Distibution of total IDU population
  output$pop_dist <- renderPlot({
    pop<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"N"]), aes(curEndemic( )[["box"]][,"N"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",N)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(pop)
  })
  
  # Distibution of death rate
  output$death_dist <- renderPlot({
    inv<- 1/curEndemic( )[["box"]][,"mu"]
    death<-  ggplot(as.data.frame(inv), aes(inv)) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",mu)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(death)
  })
  
  # Distibution of deathrate with HCV
  output$deathchronic_dist <- renderPlot({
    inv<- 1/curEndemic( )[["box"]][,"d"]
    death_HCV<-  ggplot(as.data.frame(inv), aes(inv)) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",mu + rho)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(death_HCV)
  })
  
  # Distibution of cessation rate 
  output$cess_dist <- renderPlot({
    inv<- 1/curEndemic( )[["box"]][,"ce"]
    cess<-  ggplot(as.data.frame(inv), aes(inv)) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",c)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(cess)
  })
  
  # Distibution of clearance rate
  output$clear_dist <- renderPlot({
    clearance<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"p"]), aes(curEndemic( )[["box"]][,"p"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",p)))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(clearance)
  })
  
  # Distibution of prevalence among active IDU
  output$prev_dist <- renderPlot({
    prevalence<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"pr"]), aes(curEndemic( )[["box"]][,"pr"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution prev")))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(prevalence)
  })
  
  # Distibution of new infected active
  output$YA_dist <- renderPlot({
    YA<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"YA"]), aes(curEndemic( )[["box"]][,"YA"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",Y[A])))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(YA)
  })
  
  # Distibution of new infected chronic
  output$YC_dist <- renderPlot({
    YC<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"YC"]), aes(curEndemic( )[["box"]][,"YC"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",Y[C])))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(YC)
  })
  
  # Distibution of new infected imported 
  output$YI_dist <- renderPlot({
    YI<-  ggplot(as.data.frame(curEndemic( )[["box"]][,"YI"]), aes(curEndemic( )[["box"]][,"YI"])) +
      geom_histogram(bins=25,aes(y=..count../sum(..count..)))+ 
      xlab(expression(paste("Distribution ",Y[I])))+      
      scale_y_continuous(labels = scales::percent)+ylab("Proportion")
    
    print(YI)
  })
  
  
  ##### RESULT TAB #######
  #distribution of R_in
   output$r0_in <- renderPlot({
      r0_in<-  ggplot(as.data.frame(curEndemic( )[["r0_in"]]), aes(curEndemic( )[["r0_in"]])) +
        geom_histogram(bins=25)+ xlab(expression('Reproduction number R'^"in"))
     
       print(r0_in)
   })
  
  #distribution of R_out
  output$r0_out <- renderPlot({
    r0_out<-  ggplot(as.data.frame(curEndemic( )[["r0_out"]]), aes(curEndemic( )[["r0_out"]])) +
      geom_histogram(bins=25)+ xlab(expression('Reproduction number R'^"out"))
    
    print(r0_out)
  })
  
  #distribution of time until diagnosis
   output$time_till_diagnosis_pl <- renderPlot({
     time_till_diagnosis<-  ggplot(as.data.frame(curEndemic()[["time_till_diagnosis"]]), aes(curEndemic()[["time_till_diagnosis"]]))+
       geom_histogram(bins=25)+ xlab("Time until diagnosis (years)")
     print(time_till_diagnosis)
    })
   
   #distribution of incidence
   output$incidence_pl <- renderPlot({
     incidence<-  ggplot(as.data.frame(curEndemic()[["incidence"]]), aes(curEndemic()[["incidence"]]))+
       geom_histogram(bins=25)+  xlab(expression(paste("True incidence [", yrs^{-1},']')))
     print(incidence)
   })
   
   
   # active undiagnosed with acute
   output$A <- renderText({
     paste("The estimated average number of active IDUs living undiagnosed with acute hepatitis C is ", 
           round(mean(curEndemic()[["x0"]][,"A_0"])),"[",round(quantile(curEndemic()[["x0"]][,"A_0"],probs = c(0.025,0.975))[1]),",",
           round(quantile(curEndemic()[["x0"]][,"A_0"],probs = c(0.025,0.975))[2]),"]","per year.")
   })

   #avg time until diagnosis
   output$time_till_diagnosis <- renderText({
     paste( "The estimated average  time from infection to diagnosis is ", round(mean(curEndemic()[["time_till_diagnosis"]]),2),"[",
            round(quantile(curEndemic()[["time_till_diagnosis"]],probs = c(0.025,0.975))[1],2),",",
            round(quantile(curEndemic()[["time_till_diagnosis"]],probs = c(0.025,0.975))[2],2),"]"," years.")
   })
   
   # former undiagnosed acute
   output$A_C <- renderText({
      paste("The estimated average number of  former IDUs living undiagnosed with acute hepatitis C is ", 
            round(mean(curEndemic()[["x0"]][,"A_C_0"])),"[",round(quantile(curEndemic()[["x0"]][,"A_C_0"],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["x0"]][,"A_C_0"],probs = c(0.025,0.975))[2]),"]","per year.")
    })

   #active undiagnosed chronic
   output$C <- renderText({
     paste( "The estimated average of active IDUs living undiagnosed with chronic hepatitis C is ", 
            round(mean(curEndemic()[["x0"]][,"C_0"])),"[",round(quantile(curEndemic()[["x0"]][,"C_0"],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["x0"]][,"C_0"],probs = c(0.025,0.975))[2]),"]","per year.")
   })


   # fromer undiagnosed chronic
   output$C_C <- renderText({
     paste( "The estimated average of former IDUs living undiagnosed with chronic hepatitis C is ", 
            round(mean(curEndemic()[["x0"]][,"C_C_0"])),"[",round(quantile(curEndemic()[["x0"]][,"C_C_0"],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["x0"]][,"C_C_0"],probs = c(0.025,0.975))[2]),"]","per year.")
   })

   # active imported undiagosed chronic
   output$CI<- renderText({
     paste( "The estimated average of active immigrant IDUs living undiagnosed with chronic hepatitis C is ", 
            round(mean(curEndemic()[["x0"]][,"CI_0"])),"[",round(quantile(curEndemic()[["x0"]][,"CI_0"],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["x0"]][,"CI_0"],probs = c(0.025,0.975))[2]),"]","per year.")
   })

   # former imported undiagnosed chronic
   output$CI_C<- renderText({
     paste( "The estimated average of former immigrant IDUs living undiagnosed with chronic hepatitis C is ", 
            round(mean(curEndemic()[["x0"]][,"CI_C_0"])),"[",round(quantile(curEndemic()[["x0"]][,"CI_C_0"],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["x0"]][,"CI_C_0"],probs = c(0.025,0.975))[2]),"]","per year.")
   })

   # sum undiagnosed (former and active)
   output$total_undiagnosed<- renderText({
     paste( "The estimated total average of former and active undiagnosed with hepatitis C is ", round(mean(curEndemic()[["total_undiagnosed"]])),"[",
            round(quantile(curEndemic()[["total_undiagnosed"]],probs = c(0.025,0.975))[1]),",",
            round(quantile(curEndemic()[["total_undiagnosed"]],probs = c(0.025,0.975))[2]),"]","per year.")
   })

   
    output$bar_plot <- renderPlot({
         total_ud_pie <- data.frame(
         group = c("Acute (active)", "Acute (former)", "Chronic (active)", "Chronic (former)", "Immigration (active)", "Immigration (former)"),
         value = c(mean(curEndemic()[["x0"]][,"A_0"]), mean(curEndemic()[["x0"]][,"A_C_0"]),mean(curEndemic()[["x0"]][,"C_0"]),
                   mean(curEndemic()[["x0"]][,"C_C_0"]),mean(curEndemic()[["x0"]][,"CI_0"]),mean(curEndemic()[["x0"]][,"CI_C_0"])))
         bar_plot <- ggplot(total_ud_pie, aes(x="", y=value, fill=group))+ geom_bar(width = 1, stat = "identity", position="dodge") + xlab("") +ylab("Undiagnosed cases")
         print(bar_plot)
    })
   
  # R0_in
   output$R0_in<- renderText({
      paste("The estimated  reproduction number for IDU's that get infected in the country under consideration is", round(mean(curEndemic( )[["r0_in"]]),2),"[",round(quantile(curEndemic( )[["r0_in"]],probs = c(0.025,0.975))[1],2),
            round(quantile(curEndemic( )[["r0_in"]],probs = c(0.025,0.975))[2],2),"]",".")
    })
   
   # R0_out
   output$R0_out<- renderText({
     paste("The estimated  reproduction number for imported cases is", round(mean(curEndemic( )[["r0_out"]]),2),"[",round(quantile(curEndemic( )[["r0_out"]],probs = c(0.025,0.975))[1],2),
           round(quantile(curEndemic( )[["r0_out"]],probs = c(0.025,0.975))[2],2),"]",".")
   })
   
   #yearly incidence
   output$incidence_txt<- renderText({
     paste("The true yearly incidence is", round(mean(curEndemic( )[["incidence"]]),2),"[",round(quantile(curEndemic( )[["incidence"]],probs = c(0.025,0.975))[1],2),
           round(quantile(curEndemic( )[["incidence"]],probs = c(0.025,0.975))[2],2),"]",".")
   })
   
   output$time_till_diagnosis_txt<- renderText({
     paste("The average time until diagnosis is", round(mean(curEndemic( )[["time_till_diagnosis"]]),2),"[",round(quantile(curEndemic( )[["time_till_diagnosis"]],probs = c(0.025,0.975))[1],2),
           round(quantile(curEndemic( )[["time_till_diagnosis"]],probs = c(0.025,0.975))[2],2),"]","years.")
   })
   
   #table distribution undiagnosed
   output$table_undia<- renderTable({
     frame <- rbind(Acute=c(Active=paste_meansd(curEndemic()[["x0"]][,"A_0"],digits=0), 
                            Former=paste_meansd(curEndemic()[["x0"]][,"A_C_0"],digits=0),
                            Total= paste_meansd(curEndemic()[["x0"]][,"A_0"]+curEndemic()[["x0"]][,"A_C_0"],digits=0)
                    ),
                    Chronic=c(Active=paste_meansd(curEndemic()[["x0"]][,"C_0"]+curEndemic()[["x0"]][,"CI_0"],digits=0), 
                              Former=paste_meansd(curEndemic()[["x0"]][,"C_C_0"]+curEndemic()[["x0"]][,"CI_C_0"],digits=0),
                              Total = paste_meansd(curEndemic()[["x0"]][,"C_0"]+curEndemic()[["x0"]][,"CI_0"]+ curEndemic()[["x0"]][,"C_C_0"]+
                                                     curEndemic()[["x0"]][,"CI_C_0"],digits=0)
                    ), 
                    Total=c(Active=paste_meansd(curEndemic()[["x0"]][,"A_0"]+curEndemic()[["x0"]][,"C_0"]+ curEndemic()[["x0"]][,"CI_0"],digits=0), 
                            Former=paste_meansd(curEndemic()[["x0"]][,"A_C_0"]+curEndemic()[["x0"]][,"C_C_0"]+
                                                curEndemic()[["x0"]][,"CI_C_0"],digits=0),
                            Total= paste_meansd(curEndemic()[["x0"]][,"A_0"]+curEndemic()[["x0"]][,"C_0"]+ curEndemic()[["x0"]][,"CI_0"]+
                                                curEndemic()[["x0"]][,"A_C_0"]+curEndemic()[["x0"]][,"C_C_0"]+ curEndemic()[["x0"]][,"CI_C_0"],digits=0)
                    )
                    )
     frame.data <- as.data.frame(frame)
     frame.data
   }, rownames =TRUE, digits=0)
   
   
   
   ###############   SENSITIVITY #####################

   output$table_sens<- renderTable({
     
   #parameters for which we want to carry out the sensitivity analysis
   param_vec_sens<- cbind(input$pop, 
                    input$death,
                    input$deathchronic,  
                    input$cess, 
                    clear_sim=0.26,
                    input$prev,
                    input$acute,  
                    input$chronic,
                    input$immigrants)
   
   #name the parameters 
   colnames(param_vec_sens) <- c("pop", "death", "deathchronic",  "cess",  "clear", "prev", "acute", "chronic", "immigrants")

   ## We now want to calulate the effect that every parameter has on the output so we create a for loop which
   ## loops through all possible free paparmeters


   # initializing the loop
   res_sens<- c(0,0,0,0)
   # specifiying the sensitivity level
   sens_level<- 0.1
   
   for(name in colnames(param_vec_sens)){
     # store the results in res_sens
     res_sens <- rbind(res_sens, sensitivity(sens_level,name,param_vec_sens)[["percentage"]])
   }
   
   # delete the first row which was the initalization 
   res_sens <- as.data.frame(round(res_sens[-1,],0))
   
   # add % in each cell
   res_sens <- sapply(res_sens , paste, "%", sep="")
   
   # rename the rows of the result matrix so we see which of the rows denotes which effect
   rownames(res_sens) <- cbind(" Active IDU population size", " Mortality rate",
                                                 " Mortality rate for IDUs with HCV", " Permanent IDU cessation rate",
                                                 " Clearance probability", " Prevalence among active IDU", "Diagnosed acute cases",
                                                 " Diagnosed chronic cases", " Diagnosed imported chronic cases")
  
   print(res_sens)

    },rownames = TRUE, align='r', 
   caption = "Summary of sensitivity analysis showing percentage change in time until diagnosis, 
             $R^{\\text{in}}$, number of undiagnosed and true incidence as each quantity is changed 
             from 90% to 110% of the mean value.")

  

   
   ############TREATMENT TAB
   #TODO make years to simulate an input AND fix the ggplot with mean and confidence interval
   # Preparing the pomp object
   # global variable how many years into the future should be simulated
   years_to_simulate <- 20
   # bind parameters and endemic level (and round as pomp takes only integer initial values)
   params <- reactive({ cbind(curEndemic()[["parms"]],round(curEndemic()[["x0"]])) })
   
   #initialize pomp object with the first parameter vector 
   sir_object <- reactive({ pomp_object(params()[1,], years_to_simulate)  })
   
   #transform paramer vector to dataframe
   params_df <-  reactive({ as.data.frame(params()) })
   
   # simulate the trajectory for each parameter vector 
   sir_sim <- reactive({ foreach(guess=iter(params_df(),"row"),.packages='pomp', .options.multicore=mcopts) %dopar% {
     trajectory(sir_object()[["sir"]],params= c(unlist(guess)), as.data.frame=TRUE)
   }
   })

   
   output$sim_plot <- renderPlot({ 
     
     #extract solutions to a dataframe
     sir_sim_df  <-  ldply(sir_sim(), data.frame)
     #give an id for each year to each simulation
     sir_sim_df$id <- factor(rep(seq(1:input$sample), each = years_to_simulate ))
     
     sir_sim_df%>%
       subset(select=c(time,id, S, A,C,A_C,C_C,C_dia, A_dia,C_C_dia, A_C_dia, CI, CI_C,S_NE, A_NE,C_NE,C_dia_NE, A_dia_NE,CI_NE)) %>%
       melt(id=c("time","id")) %>%
       mutate(variable = factor(variable)) %>%
       mutate(variable = recode(variable, S = "S['0,0,0']")) %>%
       mutate(variable = recode(variable, A_C = "A['0,0,1']")) %>%
       mutate(variable = recode(variable, A = "A['0,0,0']")) %>%
       mutate(variable = recode(variable, C ="C['0,0,0']")) %>%
       mutate(variable = recode(variable, C_C ="C['0,0,1']")) %>%
       mutate(variable = recode(variable, C_dia = "C['1,0,0']")) %>%
       mutate(variable = recode(variable, A_dia = "A['1,0,0']")) %>%
       mutate(variable = recode(variable, C_C_dia = "C['1,0,1']")) %>%
       mutate(variable = recode(variable, A_C_dia = "A['1,0,1']")) %>%
       mutate(variable = recode(variable, CI = "C['0,0,0']^I")) %>%
       mutate(variable = recode(variable, CI_C = "C['0,0,1']^I")) %>%
       mutate(variable = recode(variable, S_NE = "S['0,1,0']")) %>%
       mutate(variable = recode(variable, A_NE = "A['0,1,0']")) %>%
       mutate(variable = recode(variable, C_NE = "C['0,1,0']")) %>%
       mutate(variable = recode(variable, CI_NE = "C['0,1,0']^I")) %>%
       mutate(variable = recode(variable, A_dia_NE = "A['1,1,0']")) %>%
       mutate(variable = recode(variable, C_dia_NE = "C['1,1,0']")) %>%
       ggplot(aes(x=time,y=value,    group=id))+ geom_line()+
       scale_alpha_discrete(breaks=c(`TRUE`=0.2,`FALSE`=0.2))+#geom_line(data=w,color="red",aes(x=time,y=value),inherit.aes = FALSE)+
       labs(color="")+ xlab("Time [yrs]") + ylab("Individuals")+
       guides(alpha=FALSE)+
       facet_wrap(~variable,ncol=3,scales="free_y", labeller = label_parsed)+theme_bw()->tra_plot
     
    print(tra_plot)

   })
   
##TO DO### Check the equations here
   
  output$table_trea<- renderTable({  
    
    #extract solutions to a dataframe
    sir_sim_df <- ldply(sir_sim(), data.frame)
    #give an id for each year to each simulation
    sir_sim_df$id <- factor(rep(seq(1:input$sample), each = years_to_simulate ))
    
    # gelp functions to get the first and last entry of the dataframe
    tail <- sir_sim_df%>%
            subset(time==years_to_simulate) 
    head <- sir_sim_df%>%
            subset(time==1) 
    
    prevalence_sce_b <- (rowSums(head[c("A","C","C_dia", "A_dia", "CI","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                         rowSums(head[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    population_sce_b <- rowSums(head[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    incidence_sce_b <- t(head[c("S")]) * params()[,"Beta"] *(rowSums(head[c("A","C","C_dia", "A_dia", "CI")])+ 
                       params()[,"alpha"]*rowSums(head[c("A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                       rowSums(head[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])+
                       t(head[c("S_NE")] )* params()[,"alpha"]* params()[,"Beta"] * (rowSums(head[c("A","C","C_dia", "A_dia", "CI")])+ 
                       params()[,"alpha"]*rowSums(head[c("A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                       rowSums(head[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    active_undia_sce_b <- rowSums(head[c("A","C","CI", "A_NE", "C_NE", "CI_NE")])
    former_undia_sce_b <- rowSums(head[c("A_C", "C_C",  "CI_C")])
    
    total_infected_sce_b  <- rowSums(head[c("A","C","C_dia", "A_dia", "A_C", "C_C", "C_C_dia", 
                                          "A_C_dia","CI", "CI_C","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
     
    base_scenario <- data.frame( Prevalence=mean(prevalence_sce_b),
                                 Incidence=mean(incidence_sce_b),
                                 Active=mean(active_undia_sce_b),
                                 Former=mean(former_undia_sce_b),
                                 Total=mean(total_infected_sce_b),
                                 Population=mean(population_sce_b))
  
    
    prevalence_sce <- (rowSums(tail[c("A","C","C_dia", "A_dia", "CI","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                       rowSums(tail[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    population_sce <- rowSums(tail[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    incidence_sce <- t(tail[c("S")]) * params()[,"Beta"] * (rowSums(tail[c("A","C","C_dia", "A_dia", "CI")])+ 
                     params()[,"alpha"]*rowSums(tail[c("A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                     rowSums(tail[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])+
                     t(tail[c("S_NE")] )* params()[,"alpha"]* params()[,"Beta"] * (rowSums(tail[c("A","C","C_dia", "A_dia", "CI")])+ 
                     params()[,"alpha"]*rowSums(tail[c("A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")]))/
                     rowSums(tail[c("S","A","C","C_dia", "A_dia", "CI", "S_NE","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    active_undia_sce <- rowSums(tail[c("A","C","CI", "A_NE", "C_NE", "CI_NE")])
    former_undia_sce <- rowSums(tail[c("A_C", "C_C",  "CI_C")])
    
    total_infected_sce  <- rowSums(tail[c("A","C","C_dia", "A_dia", "A_C", "C_C", "C_C_dia", 
                                          "A_C_dia","CI", "CI_C","A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE")])
    
    
    scenario_with_intervention <- data.frame( Prevalence=mean(prevalence_sce),
                                              Incidence=mean(incidence_sce),
                                              Active=mean(active_undia_sce),
                                              Former=mean(former_undia_sce),
                                              Total=mean(total_infected_sce),
                                              Population=mean(population_sce))
    
    # calcualte the change in the key putcomes relative to the base scenario
    scenario_with_intervention_change <- rbind(round((scenario_with_intervention/base_scenario-1)*100,0))
    # round prevalence different from the catrgories containing individuals
    base_scenario <- cbind(Prevalence= round(base_scenario[,1],2), round(base_scenario[,-1],0))  
    #add % in each cell
    scenario_with_intervention_change <- sapply(scenario_with_intervention_change, paste, "%", sep="")
    
    # bind the base with the change in key outcomes due to interventions
    scenario_comparison <- rbind(base_scenario, scenario_with_intervention_change)
    
    #name the rows of the table
    row.names(scenario_comparison) <- c("Baseline scenario","Percentage change with selected treatment scenario")

    print(scenario_comparison)
    
   }, include.rownames = TRUE,align='rrrrrr')
  
  
})


