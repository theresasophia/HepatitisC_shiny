######################################################################
## R code generating the number of annual observed acute and chonic 
## hepatitis cases with IDU as transmission route based on the accumulated 
## number of cases between 2005-2017 made availbale by the Public Health 
## Agency of Sweden. It is assumed that the disease was endemic during that 
## time period and no significant changes in disease epidemiology have occured.
## We chose "Alternative 2" accounting for missing values in the data 
## for  the calculations for the number of cases that we use in the manuscript 
## for demonstration. 
##
## Author: Theresa Stocks <http://www.su.se/english/profiles/tstoc-1.219526>
## Affiliation: Department of Mathematics, Stockholm University, Sweden
##
## Date: 4/10/2018
######################################################################


rm(list = ls())     # clear objects  
graphics.off()      # close graphics windows

# load necessary libraries
library(readxl)

# read_excel reads both xls and xlsx files
setwd("~/Dropbox/HepatitisC/R code")

hep_dat <- as.data.frame(read_excel("hep.xlsx"))
head(hep_dat)
class(hep_dat)
hep_dat$type

#ASSUMPTION: set the unkonwns to 5 cases
hep_dat[hep_dat=="<10"]<-5
hep_dat$number <- as.numeric(hep_dat$number)
sum(hep_dat$number)


# in the following we denote for the 
# 1st letter: A= Acute, C= Chronic, M= Missing;
# 2nd letter: S= Sweden, O= Other, M= Missing; 
# 3rd letter: I= IDU, O= Other, M= Missing

ASI_t <- subset(hep_dat, type == "Acute infection" & country == "Sweden" & route == "IDU")
ASI <- sum(ASI_t$number)

ASO_t <- subset(hep_dat, type == "Acute infection" & country == "Sweden" & route == "Other")
ASO <- sum(ASO_t$number)

ASM_t <- subset(hep_dat, type == "Acute infection" & country == "Sweden" & route == "Missing, unknown route")
ASM <- sum(ASM_t$number)

AOM_t <- subset(hep_dat, type == "Acute infection" & country == "Other" & route == "Missing, unknown route")
AOM <- sum(AOM_t$number)

AOI_t <- subset(hep_dat, type == "Acute infection" & country == "Other" & route == "IDU")
AOI <- sum(AOI_t$number)

AOO_t <- subset(hep_dat, type == "Acute infection" & country == "Other" & route == "Other")
AOO <- sum(AOO_t$number)

AMI_t <- subset(hep_dat, type == "Acute infection" & country == "Missing, unknown country" & route == "IDU")
AMI <- sum(AMI_t$number)

AMM_t <- subset(hep_dat, type == "Acute infection" & country == "Missing, unknown country" & route == "Missing, unknown route")
AMM <- sum(AMM_t$number)

MOI_t <- subset(hep_dat, type=="Missing, unknown type"& country == "Other" & route == "IDU")
MOI <- sum(MOI_t$number)

MSI_t <- subset(hep_dat, type=="Missing, unknown type"& country == "Sweden" & route == "IDU")
MSI <- sum(MSI_t$number)

MSM_t <- subset(hep_dat, type=="Missing, unknown type"& country == "Sweden" & route == "Missing, unknown route")
MSM <- sum(MSM_t$number)

MMI_t <- subset(hep_dat, type=="Missing, unknown type"& country == "Missing, unknown country" & route == "IDU")
MMI <- sum(MMI_t$number)

CSI_t <- subset(hep_dat, type == "Chronic infection" & country == "Sweden" & route == "IDU")
CSI <- sum(CSI_t$number)
 
CSO_t <- subset(hep_dat, type == "Chronic infection" & country == "Sweden" & route == "Other")
CSO <- sum(CSO_t$number)

CSM_t <- subset(hep_dat, type == "Chronic infection" & country == "Sweden" & route == "Missing, unknown route")
CSM <- sum(CSM_t$number) 

CMI_t <- subset(hep_dat, type == "Chronic infection" & country == "Missing, unknown country" & route == "IDU")
CMI <- sum(CMI_t$number)

COI_t <- subset(hep_dat, type == "Chronic infection" & country == "Other" & route == "IDU")
COI <- sum(COI_t$number)

COO_t <- subset(hep_dat, type == "Chronic infection" & country == "Other" & route == "Other")
COO <- sum(COO_t$number)

COM_t <- subset(hep_dat, type == "Chronic infection" & country == "Other" & route == "Missing, unknown route")
COM <- sum(COM_t$number)

CMM_t <- subset(hep_dat, type == "Chronic infection" & country == "Missing, unknown country" & route == "Missing, unknown route")
CMM <- sum(CMM_t$number)

MOM_t <- subset(hep_dat, type == "Missing, unknown type" & country == "Other" & route == "Missing, unknown route")
MOM <- sum(MOM_t$number)

MMM_t <- subset(hep_dat, type == "Missing, unknown type" & country == "Missing, unknown country" & route == "Missing, unknown route")
MMM <- sum(MMM_t$number)

MMO_t <- subset(hep_dat, type == "Missing, unknown type" & country == "Missing, unknown country" & route == "Other")
MMO <- sum(MMO_t$number)

MSO_t <- subset(hep_dat, type == "Missing, unknown type" & country == "Sweden" & route == "Other")
MSO <- sum(MSO_t$number)

MOO_t <- subset(hep_dat, type == "Missing, unknown type" & country == "Other" & route == "Other")
MOO <- sum(MOO_t$number)

AMO_t <- subset(hep_dat, type == "Acute infection" & country == "Missing, unknown country" & route == "Other")
AMO <- sum(AMO_t$number)

CMO_t <- subset(hep_dat, type == "Chronic infection" & country == "Missing, unknown country" & route == "Other")
CMO <- sum(CMO_t$number)

missing_t <- subset(hep_dat, type == "Missing, unknown type"  | country == "Missing, unknown country" | route == "Missing, unknown route")
missing <- sum(missing_t$number)

total_cases <- sum(hep_dat$number)

# fraction of missing cases 
round(missing/total_cases,3) 


# plotting the distribution of missing cases
slices <- c(MMI ,MMO , MMM , MOI , MOM , MSI, MSM , MSO , MOO ,
              CMO , COM , CSM , CMM  , CMI ,
              AMM  ,ASM , AOM , AMI ,AMO)
lbls <- c("MMI",  "MMO",  "MMM " ,"MOI", " MOM",  "MSI",  "MSM"  ,"MSO" , "MOO" ,
            "CMO" , "COM ", "CSM" , "CMM"  , "CMI" ,
            "AMM" , "ASM" , "AOM" , "AMI" , "AMO")
pie(slices,labels = lbls,main="Distribution of missing values")

# Alternative 1- ignoring all missing values, i.e taking ACI as accumulated number of acute 
# infections in Sweden with IDU as transmission route

# Alternative 2- proportional scaling of cases to account for missing values
acute_sweden_IDU <- ASI + 
                    ASI/(ASI+ASO)*ASM + ASI/(ASI+AOI)*AMI + ASI/(ASI+CSI)*MSI + 
                    ASI/(ASI+ ASO+AOI+AOO)*AMM + ASI/(ASI+CSI+ASO+CSO)*MSM + ASI/(ASI+AOI+CSI+COI)*MMI + 
                    ASI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM


# as it is assumed that the number of cases have been constant over the past 12 years we can divide 
# by the number of years (i.e. 12 years)
round(ASI/12)
round(acute_sweden_IDU/12)

acute_other_IDU <- AOI + 
                   AOI/(AOI+AOO)*AOM + AOI/(ASI+AOI)*AMI + AOI/(AOI+COI)*MOI + 
                   AOI/(ASI+ ASO+AOI+AOO)*AMM + AOI/(AOI+COI+AOO+COO)*MOM + AOI/(ASI+AOI+CSI+COI)*MMI + 
                   AOI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM

round(acute_other_IDU/12 )
round(AOI/12)

chronic_sweden_IDU <- CSI + 
                      CSI/(CSI+CSO)*CSM + CSI/(CSI+COI)*CMI + CSI/(ASI+CSI)*MSI + 
                      CSI/(CSI+ CSO+COI+COO)*CMM + CSI/(ASI+CSI+ASO+CSO)*MSM + CSI/(ASI+AOI+CSI+COI)*MMI + 
                      CSI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM

round(CSI/12)
round(chronic_sweden_IDU/12)
 
ASI/(ASI+CSI) 
AOI/(AOI+COI) 
acute_sweden_IDU/(acute_sweden_IDU+chronic_sweden_IDU) 

chronic_other_IDU <- COI + 
                     COI/(COI+COO)*COM + COI/(CSI+COI)*CMI + COI/(AOI+COI)*MOI + 
                     COI/(CSI+ CSO+COI+COO)*CMM + COI/(AOI+COI+AOO+COO)*MOM + COI/(ASI+AOI+CSI+COI)*MMI + 
                     COI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM

round(chronic_other_IDU/12)
round(COI/12)

chronic_other_IDU/(chronic_other_IDU+chronic_sweden_IDU)
COI/(COI+CSI)


raw_data <-matrix(c(round(ASI/12),round(CSI/12),round(AOI/12),round(COI/12)),ncol=2,byrow=TRUE)
colnames(raw_data) <- c("acute","chronic")
rownames(raw_data) <- c("Sweden","Abroad")
raw_data <- as.table(raw_data)
raw_data

data_miss <- matrix(c(round(acute_sweden_IDU/12),round(chronic_sweden_IDU/12),round(acute_other_IDU/12 ),
                      round(chronic_other_IDU/12)),ncol=2,byrow=TRUE)
colnames(data_miss) <- c("acute","chronic")
rownames(data_miss) <- c("Sweden","Abroad")
data_miss <- as.table(data_miss)
data_miss

A_M_t <- subset(hep_dat, agecat4_calc_jul=="18-29 yrs" & route == "Missing, unknown route")
A_M <- sum(A_M_t$number)


# proportion of IDU cases if proportinal scaling
(acute_sweden_IDU + acute_other_IDU + chronic_sweden_IDU + chronic_other_IDU)/sum(hep_dat$number)

# proportion of IDU cases if ignoring all missing cases
(ASI+AOI+CSI+COI)/sum(hep_dat$number)



###### Alternative 3 (all swedish unknown routes of transmission are IDU)
acute_sweden_IDU <- ASI + 
  ASM + ASI/(ASI+AOI)*AMI + ASI/(ASI+CSI)*MSI + 
  ASI/(ASI+AOI)*AMM + ASI/(ASI+CSI)*MSM + ASI/(ASI+AOI+CSI+COI)*MMI + 
  ASI/(ASI+AOI+CSI+COI)*MMM


acute_other_IDU <- AOI + 
  AOI/(AOI+AOO)*AOM + AOI/(ASI+AOI)*AMI + AOI/(AOI+COI)*MOI + 
  AOI/(ASI+ ASO+AOI+AOO)*AMM + AOI/(AOI+COI+AOO+COO)*MOM + AOI/(ASI+AOI+CSI+COI)*MMI + 
  AOI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM


chronic_sweden_IDU <- CSI + 
   CSM + CSI/(CSI+COI)*CMI + CSI/(ASI+CSI)*MSI + 
  CSI/(CSI+ COI)*CMM + CSI/(ASI+CSI)*MSM + CSI/(ASI+AOI+CSI+COI)*MMI + 
  CSI/(ASI+AOI+CSI+COI)*MMM

chronic_other_IDU <- COI + 
  COI/(COI+COO)*COM + COI/(CSI+COI)*CMI + COI/(AOI+COI)*MOI + 
  COI/(CSI+ CSO+COI+COO)*CMM + COI/(AOI+COI+AOO+COO)*MOM + COI/(ASI+AOI+CSI+COI)*MMI + 
  COI/(ASI+ASO+AOI+AOO+CSI+CSO+COI+COO)*MMM

explain <- acute_sweden_IDU + acute_other_IDU + chronic_sweden_IDU + chronic_other_IDU
explain/sum(hep_dat$number)

(chronic_sweden_IDU + chronic_other_IDU)/12

explain/sum(hep_dat$number)
round(chronic_other_IDU/12)
round(COI/12)

raw_data3 <-matrix(c(round(ASI/12),round(CSI/12),round(AOI/12),round(COI/12)),ncol=2,byrow=TRUE)
colnames(raw_data) <- c("acute","chronic")
rownames(raw_data) <- c("Sweden","Abroad")
raw_data <- as.table(raw_data)
raw_data
acute_raw <- round(ASI/12)+round(AOI/12)
acute_raw


data_miss3 <- matrix(c(round(acute_sweden_IDU/12),round(chronic_sweden_IDU/12),round(acute_other_IDU/12 ),round(chronic_other_IDU/12)),ncol=2,byrow=TRUE)
colnames(data_miss3) <- c("acute","chronic")
rownames(data_miss3) <- c("Sweden","Abroad")
data_miss <- as.table(data_miss)
data_miss3 <- as.table(data_miss3)

chronic_other_IDU/(chronic_other_IDU+chronic_sweden_IDU)

acute_raw <- round(ASI/12)+round(AOI/12)
acute_raw

chronic_raw <- round(CSI/12)+round(COI/12)
chronic_raw

acute_3 <- round(acute_sweden_IDU/12)+round(acute_other_IDU/12 )
acute_3

chronic_3 <- round(chronic_sweden_IDU/12)+round(chronic_other_IDU/12 )
chronic_3
