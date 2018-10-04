######################################################################
## R code generating the functions needed for the manuscript
## "Dynamic modelling of hepatitis C transmission among injecting drug 
## users: A tool to support WHO elimination targets"
##
##
## Author: Theresa Stocks <http://www.su.se/english/profiles/tstoc-1.219526>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 13/09/2018
######################################################################



hepC <-function(x, par) {
  # This fuction is the ODE system 
  #
  # Args:
  #    x:
  #    par:
  #
  #
  # Returns:
  #     y:
  
  y <- rep(NA, length(x))
  rC<- par["YC"]/(x[2] + x[6])
  ci <- (par["YI"]- rC*x[9])/rC
  rA<- par["YA"]/(par["pr"]*par["N"] - x[2] - x[4] - x[3] - ci+ x[5])
  s <- par["N"] - par["pr"]*par["N"]
  a<-  par["pr"]*par["N"] - x[2] - x[4] - x[3] - ci
  y[1] <- x[1]*par["pr"]*s - a*(par["ce"] + par["mu"] + rA +par["g"])
  y[2] <-(1 - par["p"])*par["g"]*a - (par["ce"] + par["d"] + rC)*x[2] 
  y[3] <- rC*x[2]  - par["ce"]*x[3] - par["d"]*x[3] + (1 - par["p"])*par["g"]*x[4] +  rC*ci
  y[4] <- -par["g"]*x[4] +rA*a - par["ce"]*x[4] - par["mu"]*x[4] 
  y[5] <- par["ce"]*a - (par["g"] + par["mu"] + rA)*x[5] 
  y[6] <- (1 - par["p"])*par["g"]*x[5] + par["ce"]*x[2] - par["d"]*x[6] -rC*x[6] 
  y[7] <-  rC*x[6]  -  par["d"]*x[7] + (1 -  par["p"])* par["g"]*x[8] +  par["ce"]*x[3] + rC*x[9]
  y[8] <- par["ce"]*x[4] + rA*x[5] - par["g"]*x[8] - par["mu"]*x[8]
  y[9] <- (par["d"]+par["ce"])/par["d"]*(par["d"] + par["ce"] + rC)*ci + par["ce"]*ci - (par["d"] + rC)*x[9]
  return(y)
}
# x1=beta, x2=c, x3=c*, x4=a*, x5=a_c, x6=c_c, x7=c_c*, x8=a_c*, x9 = ci_c
# y1=beta, y2=c, y3=c*, y4=a* , y5=a_c ,y6=c_c , y7=c_c* , y8=a_c*, y9=ci_c 





help_func <- function(list,i){
  # Help function that extracts row x from the ith entry of a list k
  return(list[[i]]$x)
}


endemic <- function(pop, pop.uc, 
                    death, death.uc, 
                    deathchronic, deathchronic.uc, 
                    cess, cess.uc,  
                    prev, prev.uc,
                    acute, acute.uc, 
                    chronic, chronic.uc,
                    immigrants, immigrants.uc, 
                    sample_size, 
                    trea_a,trea_f,trea_e, 
                    alpha,e,d,rE,
                    gam=2, 
                    clear_sim=0.26, clear.uc_sim=0.001){
  
  # This function takes all parameters and associated uncertainties of the model as input and then draws a number of 
  # samples according to a truncated normal distribution.
  #
  # Args:
  #
  #
  # Returns:
  #     x0: endemic level
  #     parms: parameter vector including the estimates
  #     time_till_diagnosis: average time until diagnosis
  #     total_undiagnosed: total number of undiagnosed
  #     r0_in: reproduction number of original cases
  #     r0_out: reproduction number of imported cases
  #     box: drawn parameters 
  #     incidence: true annual incidence
 
  
  # convert sample size to numeric
  sample_size=as.numeric(sample_size)
  # bind the parameters and their uncertainties as a box
  box <- rbind(g=c(gam,0.01),
               N=c(pop,pop*as.numeric(pop.uc)), 
               mu=c(1/death,1/death*as.numeric(death.uc)), 
               d=c(1/deathchronic,1/deathchronic*as.numeric(deathchronic.uc)),
               ce=c(1/cess,1/cess* as.numeric(cess.uc)),
               p= c(clear_sim, clear.uc_sim), 
               pr=c(prev,prev*as.numeric(prev.uc)), 
               YA=c(acute,acute*as.numeric(acute.uc)), 
               YC=c(chronic,chronic*as.numeric(chronic.uc)), 
               YI=c(immigrants,immigrants*as.numeric(immigrants.uc)), 
               k=c((1/deathchronic+1/cess)/(1/deathchronic),0.001))
  
  # draws smaple_size vector combinations from a truncated normal distribution with mean and sd as specified in box
  params_sampled <- as.data.frame(apply(box,1,function(x)rtruncnorm(sample_size,a=0, mean=x[1],sd=x[2])))
  #convert df into matrix
  par <- data.matrix(params_sampled)
  
  # choose random staring value to solve systme of ODE equations
  xstart <- c(7,2416,4000,10,10,10,10,10,10)
  
  #solve the system of equations for every parameter vector seperately 
  res<- apply(par,1,nleqslv,x=xstart,fn=hepC, jac=NULL)
  #extract the solutions for each paramerter vector 
  res$x <- sapply(seq(1:sample_size), help_func, list=res)
  
  #calculate the endemic level
  x0 <-cbind(S_0 = as.numeric(par[,"N"]-par[,"pr"]*par[,"N"]),
             A_0 = as.numeric((par[,"pr"]*par[,"N"] - res$x[2,] - res$x[4,] - res$x[3,]-
                                 (par[,"YI"]- par[,"YC"]/(res$x[2,]+res$x[6,])*res$x[9,])/(par[,"YC"]/(res$x[2,]+res$x[6,])))),
             C_0 = res$x[2,],A_dia_0=res$x[4,],
             C_dia_0 = res$x[3,],A_C_0=res$x[5,],
             C_C_0 = res$x[6,],
             C_C_dia_0 = res$x[7,],
             A_C_dia_0 = res$x[8,],
             H_A_0 = 0,
             H_C_0 = 0,
             H_A_C_0 = 0,
             H_C_C_0 = 0, 
             CI_0 = as.numeric((par[,"YI"]- par[,"YC"]/(res$x[2,]+res$x[6,])*res$x[9,])/(par[,"YC"]/(res$x[2,]+res$x[6,]))) ,
             CI_C_0 = res$x[9,],
             H_CI_0 = 0, 
             H_CI_C_0 = 0,
             Pop_0 = 0,S_NE_0=0,A_NE_0=0,C_NE_0=0,CI_NE_0=0,A_dia_NE_0=0,C_dia_NE_0=0)
  
  # stop if the endemic level
  stopifnot(x0 >= 0)

  
  #parameter vector
  params <- cbind(gamma=as.numeric(par[,"g"]), #rate of leaving acute to chronic state
                  N=as.numeric(par[,"N"]), # number of active IDU
                  mu=as.numeric(par[,"mu"]), # IDU specific death rate 
                  delta=as.numeric(par[,"d"]), # increased  death rate due to HCV
                  od=0,  #overdispersion parameter
                  cess=as.numeric(par[,"ce"]), #IDU cessation rate 
                  p= as.numeric(par[,"p"]),
                  prev=as.numeric(par[,"pr"]), # HCV prevalence among active IDU in endemic state without intervention 
                  Beta=res$x[1,], # number of infectious contact
                  rA= as.numeric(par[,"YA"]/(x0[,"A_0"]+x0[,"A_C_0"])), #diagnosis rate of active IDU
                  rC= as.numeric(par[,"YC"]/(x0[,"C_0"]+x0[,"C_C_0"])), #diagnosis rate of former IDU
                  theta=as.numeric((x0[,"S_0"]*par[,"ce"] +x0[,"S_0"]*par[,"mu"]+
                                      x0[,"A_0"]*par[,"ce"]+x0[,"A_0"]*par[,"mu"]+
                                      x0[,"C_0"]*par[,"ce"] +x0[,"C_0"]*par[,"d"]+
                                      x0[,"C_dia_0"]*par[,"d"] + x0[,"C_dia_0"]*par[,"ce"]+
                                      x0[,"A_dia_0"]*par[,"ce"]+x0[,"A_dia_0"]*par[,"mu"]-
                                      (par[,"d"]+par[,"ce"]+par[,"YC"]/(x0[,"C_0"]+x0[,"C_C_0"]))*x0[,"CI_0"]
                                    +par[,"d"]*x0[,"CI_0"]+par[,"ce"]*x0[,"CI_0"])), # birth rate of active IDU
                  k=as.numeric(par[,"k"]), # proportionality constant immigration 
                  i=as.numeric((par[,"d"]+par[,"ce"]+par[,"YC"]/(x0[,"C_0"]+x0[,"C_C_0"]))*x0[,"CI_0"]), #immigration rate of active IDU
                  trea_a=trea_a, #treatmentrate active drug users not in NEP
                  trea_f=trea_f, #treatement rate former drug users
                  alpha=as.numeric(alpha), # fraction of reduced infectiousness and susceptibility
                  e=e, #enrollment rate NEP
                  trea_e=trea_e, #treatment rate NEP
                  cess_ne=as.numeric(par[,"ce"]),  #cessation rate NEP
                  rE=rE, #diagnosis rate NEP
                  d=d  #drop-out rate
  )
  
  #gives an error when negative values occur
  stopifnot(params >= 0)
  
  
  #calculate the basic reproduction number within counrty under consideration and outside
  r0_in <- c(params[,"Beta"]/(params[,"gamma"]+ params[,"cess"]+params[,"mu"])*(1+((1-params[,"p"])*params[,"gamma"])/
                                                                                   (params[,"cess"]+ params[,"delta"])))
  r0_out <- c(params[,"Beta"]/(params[,"cess"]+ params[,"delta"]))
  
  #caluclate the time since infection in each stage
  t_AI <- as.numeric(1/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"]) )
  t_AF <- as.numeric(1/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"]) + 
                       1/(params[,"gamma"]+params[,"mu"]+params[,"rA"]) )
  t_CI <- as.numeric(1/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"]) + 
                       1/(params[,"cess"]+params[,"delta"]+params[,"rC"])) 
  t_CF <- as.numeric(1/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"]) + 
                       (1-params[,"p"])*params[,"gamma"]/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"])*
                       (1/(params[,"cess"]+params[,"delta"]+params[,"rC"]) + 1/(params[,"delta"]+params[,"rC"]))+
                       params[,"cess"]/(params[,"cess"]+params[,"mu"]+ params[,"gamma"]+params[,"rA"])*
                       (1/(params[,"gamma"]+params[,"mu"]+params[,"rA"]) + 1/(params[,"delta"]+params[,"rC"])))
  
  t_CIm <- as.numeric(1/(params[,"delta"]+params[,"cess"]+params[,"rC"]))
  t_CFim <- as.numeric(1/(params[,"delta"]+params[,"cess"]+params[,"rC"]) + 1/(params[,"delta"]+params[,"rC"]))
  
  #total number of newly diagnosed each year (wihtour NEP)
  total_diagnosed <- as.numeric(params[,"rA"]*x0[,"A_0"] + params[,"rA"]*x0[,"A_C_0"] + params[,"rC"]*x0[,"C_0"] + 
                                params[,"rC"]*x0[,"C_C_0"] + params[,"rC"]*x0[,"CI_0"] + params[,"rC"]*x0[,"CI_C_0"])
  
  #caluclate the time until diagnosis
  time_till_diagnosis <- as.numeric(t_AI*(params[,"rA"]*x0[,"A_0"])/total_diagnosed + t_AF*(params[,"rA"]*x0[,"A_C_0"])/total_diagnosed + t_CI* (params[,"rC"]*x0[,"C_0"])/total_diagnosed +
                          t_CF *(params[,"rC"]*x0[,"C_C_0"])/total_diagnosed + t_CIm*(params[,"rC"]*x0[,"CI_0"])/total_diagnosed + t_CFim*(params[,"rC"]*x0[,"CI_C_0"])/total_diagnosed)
  
  #total number of undiagnosed
  total_undiagnosed<- x0[,"A_0"] + x0[,"A_C_0"] + x0[,"C_0"] + x0[,"C_C_0"] + x0[,"CI_0"] + x0[,"CI_C_0"]
  
  # true incidence
  incidence <- x0[,"S_0"] * params[,"Beta"]* (x0[,"A_0"] + x0[,"C_0"] +x0[,"A_dia_0"] +x0[,"C_dia_0"] +
                                                x0[,"CI_0"] )/(x0[,"S_0"] + x0[,"A_0"] + x0[,"C_0"] +x0[,"A_dia_0"] +x0[,"C_dia_0"] +x0[,"CI_0"] )
  #save output as a list
  list(x0=x0,parms=params,time_till_diagnosis=time_till_diagnosis, total_undiagnosed=total_undiagnosed, r0_in=r0_in,r0_out=r0_out, box=par, incidence= incidence)
}


#Helpful reporting functions      
paste_meansd <- function(x, digits = 2, na.rm = TRUE){
  paste0(round(mean(x, na.rm = na.rm), digits), " (", round(sd(x, na.rm = na.rm), digits), ")")  
}

paste_meanCI <- function(x, digits = 2, na.rm = TRUE){
  paste0(round(mean(x, na.rm = na.rm), digits), " [", round(quantile(x,probs = c(0.025,0.975))[1],digits),  ",",
         round(quantile(x,probs = c(0.025,0.975))[2],digits),"]")  
}



sensitivity  <- function(sens_level, name_param , param_vec){
  # Computes the change in key outcomes for a given parameter vector where all parameters are fixed
  # except for one which is varied plus/minus a sensitivity level
  #
  # Args:
  #   sens_level: the sensitiviy level i.e. how much the free parameter is varied (numeric)
  #   name_param: name of the free parameter that should be varied (character)
  #   param_vec: vector of all parameters at their assumed mean value (named vector)
  #
  # Returns:
  #   matrix: matrix that contains absolut values of  key outcomes for upper, mean and lower bound of parameter vector
  #   percentage: vector with the percentage change for each key outcomes (from lower to upper) if name_param is varied
  #   matrix_perc: concatinates matric and percentage
  
  
  # initializing the upper boundary vector
  value_up <- param_vec
  
  
  # we replace the free parameter name_param with value + sens_level*value
  value_up[,name_param] <- value_up[,name_param] + sens_level*value_up[,name_param]
  # we only want the mean so we fix uncertainty to as small as possible everywhere
  no_uncertainty <- 0.000001
  
  # we draw  two samples (because one does not work)
  sample_sens <- 2
  
  # we calulte the new result list for the upper boundary
  curEndemic_up <- endemic(value_up[,"pop"], no_uncertainty, 
                           value_up[,"death"], no_uncertainty, 
                           value_up[,"deathchronic"],no_uncertainty, 
                           value_up[,"cess"], no_uncertainty, 
                           value_up[,"prev"], no_uncertainty, 
                           value_up[,"acute"], no_uncertainty, 
                           value_up[,"chronic"],no_uncertainty, 
                           value_up[,"immigrants"], no_uncertainty,
                           sample_sens,0,0,0,0,0,0,0)
  
  
  # initialize the lower boundary vector
  value_lo <- param_vec
  # we replace the free parameter name_param with value - sens_level*value
  value_lo[,name_param] <- value_lo[,name_param] - sens_level*value_lo[,name_param]
  
  
  # we calulte the new reults for the lower boundary
  curEndemic_lo <- endemic(value_lo[,"pop"],no_uncertainty , 
                           value_lo[,"death"], no_uncertainty, 
                           value_lo[,"deathchronic"],no_uncertainty, 
                           value_lo[,"cess"], no_uncertainty, 
                           value_lo[,"prev"], no_uncertainty, 
                           value_lo[,"acute"], no_uncertainty, 
                           value_lo[,"chronic"], no_uncertainty, 
                           value_lo[,"immigrants"], no_uncertainty,
                           sample_sens,0,0,0,0,0,0,0)
  
  curEndemic_mean <- endemic(param_vec[,"pop"],no_uncertainty, 
                             param_vec[,"death"], no_uncertainty,
                             param_vec[,"deathchronic"], no_uncertainty, 
                             param_vec[,"cess"], no_uncertainty,
                             param_vec[,"prev"], no_uncertainty, 
                             param_vec[,"acute"], no_uncertainty, 
                             param_vec[,"chronic"], no_uncertainty, 
                             param_vec[,"immigrants"], no_uncertainty,
                             sample_sens,0,0,0,0,0,0,0)
  
  # result matrix cotaining the mean of the  key outcome if calucalted for mean, upper boundary and lower boundary vector 
  matrix <- rbind( lower=c(t=mean(curEndemic_lo$time_till_diagnosis),
                           r=mean(curEndemic_lo$r0_in),
                           u=mean(curEndemic_lo$total_undiagnosed),
                           i=mean(curEndemic_lo$incidence)),
                    mean =c(t=mean(curEndemic_mean$time_till_diagnosis),
                           r=mean(curEndemic_mean$r0_in),
                           u=mean(curEndemic_mean$total_undiagnosed),
                           i=mean(curEndemic_mean$incidence)),
                   upper=c(t=mean(curEndemic_up$time_till_diagnosis),
                           r=mean(curEndemic_up$r0_in),
                           u=mean(curEndemic_up$total_undiagnosed),
                           i=mean(curEndemic_up$incidence)))
  
  
  # we rename the cpolumes of the result matrix
  colnames(matrix) <- cbind("Time until diagnosis","Rin"," Undiagnosed","Incidence")
  
  # we calucalte the effect that the variation of the one parameter has on the quantities of interest in percentage (*100)
  # from lower to uppeer
  percentage <- 100* (matrix["upper",]/matrix["lower",]-1)
  
  # we concatinate the matrix with the percentages
  matrix_perc <- round(rbind(matrix, percent=percentage),2)
  
  # we return a list with the results: seperated the matrix, the percentages and a concatination of both
  list(matrix, percentage = percentage, matrix_perc)
}






pomp_object <- function(params, years_to_simulate){
  # Function takes a parameter vector and then creates a pomp object which contains a determinstic
  # and stochstic simulator for the HCV model as presented in the manuscript
  #
  # Args:
  #   params: vector of parameter (containing endemic level and parameters)
  #   years_to_simulate: number of years to be simulated into the future
  #
  # Returns: 
  #   sir: pomp object
  #
  #
  
  # measurement model 
  dmeas <- Csnippet("
                    if (ISNA(cases_A)) {
                    lik = (give_log) ? 0 : 1;
                    } else {
                    lik =  dnbinom_mu(cases_A, 1/od,  H_A + H_A_C, 1) + 
                    dnbinom_mu(cases_C, 1/od,  H_C + H_C_C, 1) + 
                    dnbinom_mu(cases_I, 1/od,  H_CI + H_CI_C, 1);
                    lik = (give_log) ? lik : exp(lik);}")
  
  rmeas <-  Csnippet("
                     cases_A = rnbinom_mu(1/od, H_A+H_A_C);
                     cases_C = rnbinom_mu(1/od, H_C+H_C_C);
                     cases_I = rnbinom_mu(1/od, H_CI+H_CI_C);
                     ")
  
  # process model is Markovian SIRS 
  sir.step <- Csnippet("double rate[70];
                       double dN[70];
                       
                       rate[0] = theta;
                       rate[1] = cess;
                       rate[2] = mu;
                       rate[3] = Beta*(A + C + A_dia + C_dia + CI+ alpha*(A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE))/
                       (S + A + C + A_dia + C_dia + CI+ S_NE + A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE);
                       rate[4] = e;
                       
                       rate[5] = p*gamma;
                       rate[6] = cess;
                       rate[7] = mu;
                       rate[8] = (1-p)*gamma;
                       rate[9] = rA;
                       rate[10] = e;
                       
                       rate[11] = cess;
                       rate[12] = delta;
                       rate[13] = rC;
                       rate[14] = e;
                       
                       rate[15] = cess;
                       rate[16] = delta;
                       rate[17] = trea_a;
                       rate[18] = e;
                       
                       rate[19] = cess;
                       rate[20] = mu;
                       rate[21] = (1-p)*gamma;
                       rate[22] = p*gamma;
                       rate[23] = e;
                       
                       rate[24] = mu;
                       rate[25] = p*gamma;
                       rate[26] = rA;
                       rate[27] = (1-p)*gamma;
                       
                       rate[28] = delta;
                       rate[29] = rC;
                       
                       rate[30] = delta;
                       rate[31] = trea_f;  
                       
                       rate[32] = mu;
                       rate[33] = p*gamma;
                       rate[34] = (1-p)*gamma;
                       
                       rate[35] = i;  
                       
                       rate[36] = delta;
                       rate[37] = cess;
                       rate[38] = rC;
                       rate[39] = e;
                       
                       rate[40] = k*i; 
                       
                       rate[41] = delta;
                       rate[42] = rC;
                       
                       rate[43] = mu;   
                       rate[44] = d;
                       rate[45] = alpha*Beta*(A + C + A_dia + C_dia + CI+ alpha*(A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE))/
                       (S + A + C + A_dia + C_dia + CI+ S_NE + A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE);
                       rate[46] = cess_ne; 
                       
                       rate[47] = p*gamma;
                       rate[48] = mu;
                       rate[49] = d;
                       rate[50] = (1-p)*gamma;
                       rate[51] = cess_ne;
                       rate[52] = rE;
                       
                       rate[53] = delta;
                       rate[54] = d;                      
                       rate[55] = rE;
                       rate[56] = cess_ne;
                       
                       rate[57] = cess_ne;
                       rate[58] = delta;
                       rate[59] = d;
                       rate[60] = rE;
                       
                       rate[61] = delta;
                       rate[62] = trea_e;
                       rate[63] = d;
                       rate[64] = cess_ne;
                       
                       rate[65] = cess_ne;
                       rate[66] = d;
                       rate[67] = (1-p)*gamma;
                       rate[68] = mu;
                       rate[69] = p*gamma;
                       
                       
                       dN[0] = rpois(rate[0]*dt); // births are Poisson
                       reulermultinom(4, S, &rate[1], dt, &dN[1]);
                       reulermultinom(6, A, &rate[5], dt, &dN[5]);
                       reulermultinom(4, C, &rate[11], dt, &dN[11]);
                       reulermultinom(4, C_dia, &rate[15], dt, &dN[15]);
                       reulermultinom(5, A_dia, &rate[19], dt, &dN[19]);
                       reulermultinom(4, A_C, &rate[24], dt, &dN[24]);
                       reulermultinom(2, C_C, &rate[28], dt, &dN[28]);
                       reulermultinom(2, C_C_dia, &rate[30], dt, &dN[30]);
                       reulermultinom(3, A_C_dia, &rate[32], dt, &dN[32]);
                       
                       reulermultinom(4, CI, &rate[36], dt, &dN[36]);
                       reulermultinom(2, CI_C, &rate[41], dt, &dN[41]);
                       
                       reulermultinom(4, S_NE, &rate[43], dt, &dN[43]);
                       reulermultinom(6, A_NE, &rate[47], dt, &dN[47]);
                       reulermultinom(4, C_NE, &rate[53], dt, &dN[53]);
                       reulermultinom(4, CI_NE, &rate[57], dt, &dN[57]);
                       reulermultinom(4, C_dia_NE, &rate[61], dt, &dN[61]);
                       reulermultinom(5, A_dia_NE, &rate[65], dt, &dN[65]);
                       
                       dN[35] = rpois(rate[35]*dt); // births are Poisson
                       dN[40] = rpois(rate[40]*dt); // births are Poisson
                       
                       S += dN[0] - dN[4] + dN[44] - dN[2] - dN[3] + dN[5] + dN[22] - dN[1] + dN[17];
                       A += dN[3] - dN[10] - dN[7] + dN[49] - dN[8] - dN[9] - dN[6] - dN[5];
                       C += dN[8] - dN[14] - dN[12] + dN[54] - dN[13] - dN[11]; 
                       
                       C_dia += dN[21] - dN[18] - dN[16] +  dN[63] + dN[13] + dN[38] - dN[17] - dN[15];
                       A_dia += -dN[22] - dN[20] - dN[23] + dN[9] - dN[21] + dN[66] - dN[19];
                       
                       A_C += dN[6] - dN[27] - dN[26] - dN[24] + dN[51] - dN[25];
                       C_C += dN[27] + dN[11] - dN[29] - dN[28] + dN[56]; //immigration
                       
                       C_C_dia += dN[34] + dN[29] + dN[15] +  dN[42] - dN[31] - dN[30] + dN[64];
                       A_C_dia += dN[26] + dN[19] - dN[34] - dN[32] + dN[65] - dN[33];
                       
                       CI += dN[35] - dN[39] - dN[36] + dN[59] - dN[38] - dN[37];  
                       CI_C += dN[40] + dN[37] - dN[42] - dN[41] + dN[57];
                       
                       S_NE += -dN[46] + dN[4] - dN[43] - dN[44] - dN[45] + dN[47] + dN[69] + dN[62];
                       A_NE += dN[45] + dN[10] - dN[48] - dN[49] - dN[50] - dN[52] - dN[51] - dN[47];
                       C_NE += dN[50] + dN[14] - dN[53] - dN[54] - dN[55] - dN[56]; 
                       
                       C_dia_NE += dN[18]    - dN[61] + dN[55] +  dN[60] - dN[62] - dN[63] - dN[64] + dN[67];
                       A_dia_NE += -dN[69]   + dN[23] - dN[68] + dN[52] - dN[67] - dN[66] - dN[65];
                       CI_NE += dN[39] - dN[58] - dN[59] - dN[60] - dN[57];  
                       
                       H_A += dN[9];
                       H_C += dN[13];
                       H_A_C += dN[26];
                       H_C_C += dN[29];
                       H_CI += dN[38];
                       H_CI_C += dN[42];
                       
                       Pop = S + A + C + A_dia + C_dia + CI + S_NE + A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE;
                       ")
  
  
  # # deterministic skeleton
  sir.skel <- Csnippet("double rate[70];
                        double term[70];
  
                        rate[0] = theta;
                        rate[1] = cess;
                        rate[2] = mu;
                        rate[3] = Beta*(A + C + A_dia + C_dia + CI+ alpha*(A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE))/
                                  (S + A + C + A_dia + C_dia + CI+ S_NE + A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE);
                        rate[4] = e;
  
                        rate[5] = p*gamma;
                        rate[6] = cess;
                        rate[7] = mu;
                        rate[8] = (1-p)*gamma;
                        rate[9] = rA;
                        rate[10] = e;
  
                        rate[11] = cess;
                        rate[12] = delta;
                        rate[13] = rC;
                        rate[14] = e;
  
                        rate[15] = cess;
                        rate[16] = delta;
                        rate[17] = trea_a;
                        rate[18] = e;
  
                        rate[19] = cess;
                        rate[20] = mu;
                        rate[21] = (1-p)*gamma;
                        rate[22] = p*gamma;
                        rate[23] = e;
  
                        rate[24] = mu;
                        rate[25] = p*gamma;
                        rate[26] = rA;
                        rate[27] = (1-p)*gamma;
  
                        rate[28] = delta;
                        rate[29] = rC;
  
                        rate[30] = delta;
                        rate[31] = trea_f;  
  
                        rate[32] = mu;
                        rate[33] = p*gamma;
                        rate[34] = (1-p)*gamma;
  
                        rate[35] = i;   
                        rate[36] = delta;
                        rate[37] = cess;
                        rate[38] = rC;
                        rate[39] = e;
  
                        rate[40] = k*i;   
                        rate[41] = delta;
                        rate[42] = rC;
  
                        rate[43] = mu;   
                        rate[44] = d;
                        rate[45] = alpha*Beta*(A + C + A_dia + C_dia + CI+ alpha*(A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE))/
                                  (S + A + C + A_dia + C_dia + CI+ S_NE + A_NE + C_NE + A_dia_NE + C_dia_NE + CI_NE);
                        rate[46] = cess_ne; 
  
                        rate[47] = p*gamma;
                        rate[48] = mu;
                        rate[49] = d;
                        rate[50] = (1-p)*gamma;
                        rate[51] = cess_ne;
                        rate[52] = rE;
  
                        rate[53] = delta;
                        rate[54] = d;                      
                        rate[55] = rE;
                        rate[56] = cess_ne;
  
                        rate[57] = cess_ne;
                        rate[58] = delta;
                        rate[59] = d;
                        rate[60] = rE;
  
                        rate[61] = delta;
                        rate[62] = trea_e;
                        rate[63] = d;
                        rate[64] = cess_ne;
  
                        rate[65] = cess_ne;
                        rate[66] = d;
                        rate[67] = (1-p)*gamma;
                        rate[68] = mu;
                        rate[69] = p*gamma;
      
  // compute the several terms
                        term[0] = rate[0];
  
                        term[1] = rate[1] * S;
                        term[2] = rate[2] * S;
                        term[3] = rate[3] * S;
                        term[4] = rate[4] * S;
  
                        term[5] = rate[5] * A;
                        term[6] = rate[6] * A;
                        term[7] = rate[7] * A;
                        term[8] = rate[8] * A;
                        term[9] = rate[9] * A;
                        term[10] = rate[10] * A;
  
                        term[11] = rate[11] * C;
                        term[12] = rate[12] * C;
                        term[13] = rate[13] * C;
                        term[14] = rate[14] * C;
  
                        term[15] = rate[15] * C_dia;
                        term[16] = rate[16] * C_dia;
                        term[17] = rate[17] * C_dia;
                        term[18] = rate[18] * C_dia;
  
                        term[19] = rate[19] * A_dia;
                        term[20] = rate[20] * A_dia;
                        term[21] = rate[21] * A_dia;
                        term[22] = rate[22] * A_dia;
                        term[23] = rate[23] * A_dia;
  
                        term[24] = rate[24] * A_C;
                        term[25] = rate[25] * A_C;
                        term[26] = rate[26] * A_C;
                        term[27] = rate[27] * A_C;
  
                        term[28] = rate[28] * C_C;
                        term[29] = rate[29] * C_C;
  
                        term[30] = rate[30] * C_C_dia;
                        term[31] = rate[31] * C_C_dia; 
  
                        term[32] = rate[32] * A_C_dia;
                        term[33] = rate[33] * A_C_dia;
                        term[34] = rate[34] * A_C_dia;
  
                        term[35] = rate[35];
  
                        term[36] = rate[36] * CI;
                        term[37] = rate[37] * CI;
                        term[38] = rate[38] * CI;
                        term[39] = rate[39] * CI;
  
                        term[40] = rate[40];
  
                        term[41] = rate[41] * CI_C;
                        term[42] = rate[42] * CI_C;
  
                        term[43] = rate[43] * S_NE;
                        term[44] = rate[44] * S_NE;
                        term[45] = rate[45] * S_NE;
                        term[46] = rate[46] * S_NE;
  
                        term[47] = rate[47] * A_NE;
                        term[48] = rate[48] * A_NE;
                        term[49] = rate[49] * A_NE;
                        term[50] = rate[50] * A_NE;
                        term[51] = rate[51] * A_NE;
                        term[52] = rate[52] * A_NE;
  
                        term[53] = rate[53] * C_NE;
                        term[54] = rate[54] * C_NE;
                        term[55] = rate[55] * C_NE;
                        term[56] = rate[56] * C_NE;
  
                        term[57] = rate[57] * CI_NE;
                        term[58] = rate[58] * CI_NE;
                        term[59] = rate[59] * CI_NE;
                        term[60] = rate[60] * CI_NE;
  
                        term[61] = rate[61] * C_dia_NE;
                        term[62] = rate[62] * C_dia_NE;
                        term[63] = rate[63] * C_dia_NE;
                        term[64] = rate[64] * C_dia_NE;
  
                        term[65] = rate[65] * A_dia_NE;
                        term[66] = rate[66] * A_dia_NE;
                        term[67] = rate[67] * A_dia_NE;
                        term[68] = rate[68] * A_dia_NE;
                        term[69] = rate[69] * A_dia_NE;
  
  
                        DS = term[0] - term[1] - term[2] - term[3] - term[4] + term[5] + term[22] + term[17] + term[44];
                        DA = term[3] - term[10] - term[7] + term[49] - term[8] - term[9]- term[6] -term[5];
                        DC = term[8] - term[14] - term[12] + term[54] - term[13] - term[11];
                        DC_dia = term[21] + term[63] - term[16] -  term[18] + term[13] + term[38]- term[17] - term[15];
                        DA_dia = -term[19]   + term[66] - term[22] - term[20] - term[23] + term[9] - term[21];
                        DA_C = -term[25] + term[6] - term[27] - term[26] - term[24] + term[51];
                        DC_C = term[27] + term[11] - term[29] - term[28] + term[56];
                        DC_C_dia = term[29]  + term[15] + term[42] - term[31]  - term[30] + term[64] + term[34];
                        DA_C_dia = term[26]   + term[19] - term[34] - term[32] + term[65] - term[33];
                        DCI = term[35] - term[39] - term[36] + term[59] - term[38] - term[37];
                        DCI_C = term[37] - term[42] - term[41] + term[57] + term[40];
  
                        DS_NE = term[62]  + term[4] - term[43] - term[44] - term[45] + term[47] + term[69] - term[46];
                        DA_NE = term[45] + term[10] - term[48] - term[49] - term[50] - term[52] - term[51] - term[47];
                        DC_NE = term[50] + term[14] - term[53] - term[54] - term[55] - term[56];
                        DCI_NE = term[39] - term[58] - term[59] - term[60] - term[57];
                        DA_dia_NE = -term[69] + term[23] - term[68] + term[52] - term[67] - term[66] - term[65];
                        DC_dia_NE = term[67] + term[18] - term[61] + term[55] + term[60] - term[62] - term[63] - term[64];
  
                        DH_A = term[9];
                        DH_C = term[13];
                        DH_A_C = term[26];
                        DH_C_C = term[29];
                        DH_CI =  term[38];
                        DH_CI_C =  term[42];
  
                        DPop = S + A + C + A_dia + C_dia + CI;
  ")
  

  # create an empty data frame for simulation
  data.frame(time=seq(1:years_to_simulate), cases_A=rep(NA,years_to_simulate),cases_C=rep(NA,years_to_simulate),cases_I=rep(NA,years_to_simulate))%>%
    arrange(time) -> dat
  
  
  init <- Csnippet("S = S_0;
                    A = A_0;
                    C = C_0;
                    C_dia = C_dia_0;
                    A_dia = A_dia_0 ;
                    A_C = A_C_0;
                    C_C = C_C_0;
                    A_C_dia = A_C_dia_0;
                    C_C_dia = C_C_dia_0;
                   
                    S_NE = S_NE_0;
                    A_NE = A_NE_0;
                    C_NE = C_NE_0;
                    CI_NE = CI_NE_0;
                    A_dia_NE = A_dia_NE_0;
                    C_dia_NE = C_dia_NE_0;
                   
                    H_A_C = H_A_C_0;
                    H_C_C = H_C_C_0;
                    H_A = H_A_0;
                    H_C = H_C_0;
                    CI = CI_0;
                    CI_C = CI_C_0;
                    H_CI = H_CI_0;
                    H_CI_C = H_CI_C_0;
                    Pop = Pop_0;
                    ")
  
  pomp(data = dat,
       times="time",
       t0=1,
       dmeasure = dmeas,
       rmeasure = rmeas,
       rprocess = euler.sim(step.fun = sir.step, delta.t = 1/50),
       statenames = c("S", "A", "C","C_dia", "A_dia", "A_C", "C_C","C_C_dia", "A_C_dia", 
                      "H_A_C","H_C_C",  "H_A","H_C", "CI", "CI_C","H_CI", "H_CI_C","Pop", "S_NE", "A_NE", "C_NE", "CI_NE", "A_dia_NE", "C_dia_NE"),
       paramnames = names(params),
       zeronames=c("H_A","H_C","H_A_C","H_C_C","H_CI", "H_CI_C"),
       skeleton=vectorfield(sir.skel),
       initializer=init,
       params = params
  ) -> sir
  
  list(sir=sir)
  
}








