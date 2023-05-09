###Sensitivity of MCV1 timeliness###
##MCV1=9mnths,mcv2=18 mnths######mainmode3
##MCV1=10mnths,mcv2=18 mnths######mainmode4

###Sensitivity of the cutoff estimates###
##MCV1=9mnths,mcv2=18 mnths######mainmode5
##beta1<=13mnths,13=beta2<=20mnths
##MCV1=10mnths,mcv2=18 mnths######mainmode6
##beta1<=15mnths,12=beta2<=22mnths

rm(list=ls())
library(deSolve)
library(tidyverse)
library(rootSolve)
library(readxl)
library(MASS)
require(MCMCvis)
require(coda)
source("Newscripts/Plots.R")
source("Newscripts/Casedatafunc.R")

source("Scripts/modelfunctions.R")

##generate age-struc.monthly from 0-24 mnths,then yearly from 2yrs to 14 yrs
Age_struc<-as.data.frame(c(seq(0,24,by=1),seq(36,168,by=12)))
names(Age_struc)<-"lower_lim"
Age_struc=Age_struc %>% mutate(Upper_lim=c(seq(1,24,by=1),seq(35,180,by=12))) %>% 
  mutate(newagegrp=paste(lower_lim,"_",Upper_lim)) %>% 
  mutate(upper_bound=ifelse(lower_lim<24,Upper_lim,Upper_lim+1)) %>% 
  mutate(diff=upper_bound-lower_lim) %>% 
  mutate(age.out=1/diff) %>% filter(lower_lim!=168)

source("Scripts/modelfunctions.R")
df=Age_struc
noagegps=dim(df)[1]
agegp_len<-df$diff
state_prop2=c(agegp_len*100,rep(0,noagegps),rep(0,noagegps),rep(0,noagegps),rep(0,noagegps)) ##are we replacing all state_props or just ps

ageing=1/agegp_len
##Coverage estimates
##MCV1
MCV1_data<-read_excel("Data/Coverage.xlsx",sheet="MCV1")
##MCV2
##1st MCV2 is 2013
MCV2_data<-read_excel("Data/Coverage.xlsx",sheet="MCV2")
##Cases
Case_data<-read_excel("Data/Coverage.xlsx",sheet="Cases")
###Average cases over the entire period###
av_cases<-sum(Case_data$Cases)/nrow(Case_data)
##########Observed seroprevalence###
load("Data/obs_sero")
df_obs<-obs_sero
#############generate RC for MCV1=9mnths,MCV2=18mnths####
load("Data/mainmod3")
df_chain=mainmod3$chain[10000:50000,]
df_chain= df_chain%>% as.data.frame()
names(df_chain) =c("Beta1","Beta2","Beta3","omega")
df_chain_new=df_chain
mcmcMatrix=as.matrix(df_chain_new)
numSamples = 1000

outDf_RCall<- matrix(NA,nrow=numSamples, ncol = 4)
for (mm in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  
  Beta1hi <- mcmcMatrix[randomNumber,"Beta1"]
  Beta2hi <- mcmcMatrix[randomNumber, "Beta2"]
  Beta3hi <- mcmcMatrix[randomNumber,"Beta3"]
  omegahi <- mcmcMatrix[randomNumber,"omega"]
  
  ##########nationmode8
  cov96=MCV1_data[1,3] %>% pull("MCV1")
  cov96=cov96
  case96=Case_data[1,3] %>% pull("Cases")
  params.list96=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov96,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case96,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  
  model96=as.data.frame(ode(y=state_prop2,
                            times=seq(0,12,1),
                            parms=params.list96,
                            fun=SIS_ode))
  
  
  ##ode(1997)
  cov97=MCV1_data[2,3] %>% pull("MCV1")
  cov97=cov97
  case97=Case_data[2,3] %>% pull("Cases")
  start97<-as.data.frame(t(model96[13,-1])) %>% pull("13")
  params.list97=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov97,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case97,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model97=as.data.frame(ode(y=start97,
                            times=seq(0,12,1),
                            parms=params.list97,
                            fun=SIS_ode))
  
  
  ##ode(1998)
  cov98=MCV1_data[3,3] %>% pull("MCV1")
  case98=Case_data[3,3] %>% pull("Cases")
  cov98=cov98
  case98=Case_data[3,3] %>% pull("Cases")
  start98<-as.data.frame(t(model97[13,-1])) %>% pull("13")
  params.list98=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov98,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case98,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model98=as.data.frame(ode(y=start98,
                            times=seq(0,12,1),
                            parms=params.list98,
                            fun=SIS_ode))
  
  ##ode(1999)
  cov99=MCV1_data[4,3] %>% pull("MCV1")
  cov99=cov99
  case99=Case_data[4,3] %>% pull("Cases")
  start99<-as.data.frame(t(model98[13,-1])) %>% pull("13")
  params.list99=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov99,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case99,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model99=as.data.frame(ode(y=start99,
                            times=seq(0,12,1),
                            parms=params.list99,
                            fun=SIS_ode))
  
  ##ode(2000)
  cov00=MCV1_data[5,3] %>% pull("MCV1")
  cov00=cov00
  case00=Case_data[5,3] %>% pull("Cases")
  start00<-as.data.frame(t(model99[13,-1])) %>% pull("13")
  params.list00=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov00,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case00,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model00=as.data.frame(ode(y=start00,
                            times=seq(0,12,1),
                            parms=params.list00,
                            fun=SIS_ode))
  
  ##ode(2001)
  cov01=MCV1_data[6,3] %>% pull("MCV1")
  cov01=cov01
  case01=Case_data[6,3] %>% pull("Cases")
  start01<-as.data.frame(t(model00[13,-1])) %>% pull("13")
  params.list01=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov01,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case01,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model01=as.data.frame(ode(y=start01,
                            times=seq(0,12,1),
                            parms=params.list01,
                            fun=SIS_ode))
  
  
  ##ode(2002)
  ##1st SIA in June 2002 targetting 9m-14yrs with a cov of 94%
  cov02=MCV1_data[7,3] %>% pull("MCV1")
  cov02=cov02
  case02=Case_data[7,3] %>% pull("Cases")
  
  ###Nationmode8
  SIA02<-c(rep(0,8),rep(0.94,28))
  SIA02a<-c(rep(0,8),rep(0.94,3),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA02b<-c(rep(0,11),rep(0.94,7),rep(0,18))##vaccine failure of one dose given at 12-18months
  SIA02c<-c(rep(0,18),rep(0.94,18))##vaccine failure of one dose given at >18months
  
  
  start02b<-as.data.frame(t(model01[13,-1])) %>% pull("13")
  params.list02=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov02,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case02,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model02b=as.data.frame(ode(y=start02b,
                             times=seq(0,5,1), ##run it to end may
                             parms=params.list02,
                             fun=SIS_ode))
  
  myoutput02<-model02b[6,]
  startduring02<-as.data.frame(t(myoutput02[,-1])) %>% pull("6")
  
  ##during SIA
  ##94% SIA,68%mcv1 and 0mcv2
  
  
  S=startduring02[1:noagegps]-((startduring02[1:noagegps]*Beta1hi*SIA02a)+
                                 (startduring02[1:noagegps]*Beta2hi*SIA02b)+
                                 (startduring02[1:noagegps]*Beta3hi*SIA02c))
  NI=startduring02[(1:noagegps)+noagegps]
  MCV1 <-startduring02[(1:noagegps)+2*noagegps]
  MCV2 <-startduring02[(1:noagegps)+3*noagegps]
  SIA <- startduring02[(1:noagegps)+4*noagegps]+
    ((startduring02[1:noagegps]*Beta1hi*SIA02a)+
       (startduring02[1:noagegps]*Beta2hi*SIA02b)+
       (startduring02[1:noagegps]*Beta3hi*SIA02c))
  
  
  startafter02<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model02a=as.data.frame(ode(y=startafter02,
                             times=seq(6,12,1), ##run from july to Aug
                             parms=params.list02,
                             fun=SIS_ode))
  
  
  model02<-rbind.data.frame(model02b,model02a)
  
  ##ode(2003)
  cov03=MCV1_data[8,3] %>% pull("MCV1")
  cov03=cov03
  case03=Case_data[8,3] %>% pull("Cases")
  start03<-as.data.frame(t(model02[13,-1])) %>% pull("13")
  params.list03=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov03,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case03,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model03=as.data.frame(ode(y=start03,
                            times=seq(0,12,1),
                            parms=params.list03,
                            fun=SIS_ode))
  
  ##ode(2004)
  cov04=MCV1_data[9,3] %>% pull("MCV1")
  cov04=cov04
  case04=Case_data[9,3] %>% pull("Cases")
  start04<-as.data.frame(t(model03[13,-1])) %>% pull("13")
  params.list04=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov04,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case04,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model04=as.data.frame(ode(y=start04,
                            times=seq(0,12,1),
                            parms=params.list04,
                            fun=SIS_ode))
  
  ##ode(2005)
  cov05=MCV1_data[10,3] %>% pull("MCV1")
  cov05=cov05
  case05=Case_data[10,3] %>% pull("Cases")
  start05<-as.data.frame(t(model04[13,-1])) %>% pull("13")
  params.list05=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov05,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case05,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model05=as.data.frame(ode(y=start05,
                            times=seq(0,12,1),
                            parms=params.list05,
                            fun=SIS_ode))
  
  
  ##ode(2006)
  ##2nd SIA in JuLY 2002 targetting 9m-5yrs with a cov 90%
  cov06=MCV1_data[11,3] %>% pull("MCV1")
  cov06=cov06
  case06=Case_data[11,3] %>% pull("Cases")
  
  SIA06<-c(rep(0,8),rep(0.9,19),rep(0,9))
  SIA06a<-c(rep(0,8),rep(0.9,3),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA06b<-c(rep(0,11),rep(0.9,7),rep(0,18))##vaccine failure of one dose given at <=18months
  SIA06c<-c(rep(0,18),rep(0.9,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  start06b<-as.data.frame(t(model05[13,-1])) %>% pull("13")
  params.list06=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov06,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case06,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model06b=as.data.frame(ode(y=start06b,
                             times=seq(0,6,1),
                             parms=params.list06,
                             fun=SIS_ode))
  ##during SIA
  myoutput06<-model06b[7,]
  startduring06<-as.data.frame(t(myoutput06[,-1])) %>% pull("7")
  ##implement SIA
  ##90% SIA,67%mcv1 and 0mcv2
  
  S=startduring06[1:noagegps]-((startduring06[1:noagegps]*Beta1hi*SIA06a)+
                                 (startduring06[1:noagegps]*Beta2hi*SIA06b)+
                                 (startduring06[1:noagegps]*Beta3hi*SIA06c))
  NI=startduring06[(1:noagegps)+noagegps]
  MCV1 <-startduring06[(1:noagegps)+2*noagegps]
  MCV2 <-startduring06[(1:noagegps)+3*noagegps]
  SIA <- startduring06[(1:noagegps)+4*noagegps]+
    ((startduring06[1:noagegps]*Beta1hi*SIA06a)+
       (startduring06[1:noagegps]*Beta2hi*SIA06b)+
       (startduring06[1:noagegps]*Beta3hi*SIA06c))
  
  startafter06<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model06a=as.data.frame(ode(y=startafter06,
                             times=seq(7,12,1), ##run from july to Dec
                             parms=params.list06,
                             fun=SIS_ode))
  
  model06<-rbind.data.frame(model06b,model06a)
  
  ##ode(2007)
  cov07=MCV1_data[12,3] %>% pull("MCV1")
  cov07=cov07
  case07=Case_data[12,3] %>% pull("Cases")
  start07<-as.data.frame(t(model06[13,-1])) %>% pull("13")
  params.list07=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov07,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case07,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model07=as.data.frame(ode(y=start07,
                            times=seq(0,12,1),
                            parms=params.list07,
                            fun=SIS_ode))
  
  
  ##ode(2008)
  cov08=MCV1_data[13,3] %>% pull("MCV1")
  cov08=cov08
  case08=Case_data[13,3] %>% pull("Cases")
  case08= case08
  start08<-as.data.frame(t(model07[13,-1])) %>% pull("13")
  params.list08=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov08,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case08,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model08=as.data.frame(ode(y=start08,
                            times=seq(0,12,1),
                            parms=params.list08,
                            fun=SIS_ode))
  
  ##ode(2009)
  ##3rd SIA in September 2009 targeting 9m-5yrs with a cov 82%
  
  cov09=MCV1_data[14,3] %>% pull("MCV1")
  cov09=cov09
  case09=Case_data[14,3] %>% pull("Cases")
  case09= case09
  
  
  ###Nationmode8
  SIA09<-c(rep(0,8),rep(0.82,19),rep(0,9))
  SIA09a<-c(rep(0,8),rep(0.82,3),rep(0,25))##vaccine failure of one dose given at <11mnths
  SIA09b<-c(rep(0,11),rep(0.82,7),rep(0,18))##vaccine failure of one dose given at 10-18mnths
  SIA09c<-c(rep(0,18),rep(0.82,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start09b<-as.data.frame(t(model08[13,-1])) %>% pull("13")
  params.list09=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov09,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case09,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model09b=as.data.frame(ode(y=start09b,
                             times=seq(0,8,1),
                             parms=params.list09,
                             fun=SIS_ode))
  ##during SIA
  myoutput09<-model09b[9,]
  startduring09<-as.data.frame(t(myoutput09[,-1])) %>% pull("9")
  ##implement SIA
  ##82% SIA,77%mcv1 and 0mcv2
  
  S=startduring09[1:noagegps]-((startduring09[1:noagegps]*Beta1hi*SIA09a)+
                                 (startduring09[1:noagegps]*Beta2hi*SIA09b)+
                                 (startduring09[1:noagegps]*Beta3hi*SIA09c))
  NI=startduring09[(1:noagegps)+noagegps]
  MCV1 <-startduring09[(1:noagegps)+2*noagegps]
  MCV2 <-startduring09[(1:noagegps)+3*noagegps]
  SIA <- startduring09[(1:noagegps)+4*noagegps]+
    ((startduring09[1:noagegps]*Beta1hi*SIA09a)+
       (startduring09[1:noagegps]*Beta2hi*SIA09b)+
       (startduring09[1:noagegps]*Beta3hi*SIA09c))
  
  
  startafter09<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model09a=as.data.frame(ode(y=startafter09,
                             times=seq(9,12,1), ##run from july to Aug
                             parms=params.list09,
                             fun=SIS_ode))
  
  model09<-rbind.data.frame(model09b,model09a)
  
  
  ##ode(2010)
  cov10=MCV1_data[15,3] %>% pull("MCV1")
  cov10=cov10
  case10=Case_data[15,3] %>% pull("Cases")
  case10= case10
  start10<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  params.list10=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov10,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case10,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model10=as.data.frame(ode(y=start10,
                            times=seq(0,12,1),
                            parms=params.list10,
                            fun=SIS_ode))
  
  ##ode(2011)
  cov11=MCV1_data[16,3] %>% pull("MCV1")
  cov11=cov11
  case11=Case_data[16,3] %>% pull("Cases")
  case11=case11
  start11<-as.data.frame(t(model10[13,-1])) %>% pull("13")
  params.list11=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov11,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case11,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model11=as.data.frame(ode(y=start11,
                            times=seq(0,12,1),
                            parms=params.list11,
                            fun=SIS_ode))
  
  ##ode(2012)
  ##4th SIA in November 2012 targeting 9m-5yrs with a cov 92%
  cov12=MCV1_data[17,3] %>% pull("MCV1")
  cov12=cov12
  case12=Case_data[17,3] %>% pull("Cases")
  case12=case12
  
  ###Nationmode8
  SIA12<-c(rep(0,8),rep(0.92,19),rep(0,9))
  SIA12a<-c(rep(0,8),rep(0.92,3),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA12b<-c(rep(0,11),rep(0.92,7),rep(0,18))##vaccine failure of one dose given at 12-18mnths
  SIA12c<-c(rep(0,18),rep(0.92,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start12b<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  params.list12=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov12,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case12,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model12b=as.data.frame(ode(y=start12b,
                             times=seq(0,10,1),
                             parms=params.list12,
                             fun=SIS_ode))
  ##during SIA
  myoutput12<-model12b[11,]
  startduring12<-as.data.frame(t(myoutput12[,-1])) %>% pull("11")
  ##implement SIA
  S=startduring12[1:noagegps]-((startduring12[1:noagegps]*Beta1hi*SIA12a)+
                                 (startduring12[1:noagegps]*Beta2hi*SIA12b)+
                                 (startduring12[1:noagegps]*Beta3hi*SIA12c))
  NI=startduring12[(1:noagegps)+noagegps]
  MCV1 <-startduring12[(1:noagegps)+2*noagegps]
  MCV2 <-startduring12[(1:noagegps)+3*noagegps]
  SIA <- startduring12[(1:noagegps)+4*noagegps]+
    ((startduring12[1:noagegps]*Beta1hi*SIA12a)+
       (startduring12[1:noagegps]*Beta2hi*SIA12b)+
       (startduring12[1:noagegps]*Beta3hi*SIA12c))
  
  
  startafter12<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model12a=as.data.frame(ode(y=startafter12,
                             times=seq(11,12,1), ##run from july to Aug
                             parms=params.list12,
                             fun=SIS_ode))
  
  model12<-rbind.data.frame(model12b,model12a)
  
  ##ode(2013)
  cov13=MCV1_data[18,3] %>% pull("MCV1")
  cov13=cov13
  case13=Case_data[18,3] %>% pull("Cases")
  case13=case13
  start13<-as.data.frame(t(model12[13,-1])) %>% pull("13")
  params.list13=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov13,rep(0,27)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case13,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model13=as.data.frame(ode(y=start13,
                            times=seq(0,12,1),
                            parms=params.list13,
                            fun=SIS_ode))
  
  ##ode(2014)
  ##introduce 2nd dose
  cov14=MCV1_data[19,3] %>% pull("MCV1")
  cov14=cov14
  case14=Case_data[19,3] %>% pull("Cases")
  case14=case14
  cov2_14=MCV2_data[19,3] %>% pull("MCV2")
  cov2_14=cov2_14
  start14<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  params.list14=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov14,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_14,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case14,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model14=as.data.frame(ode(y=start14,
                            times=seq(0,12,1),
                            parms=params.list14,
                            fun=SIS_ode))
  
  
  ##ode(2015)
  cov15=MCV1_data[20,3] %>% pull("MCV1")
  cov15=cov15
  case15=Case_data[20,3] %>% pull("Cases")
  cov2_15=MCV2_data[20,3] %>% pull("MCV2")
  cov2_15=cov2_15
  start15<-as.data.frame(t(model14[13,-1])) %>% pull("13")
  params.list15=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov15,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_15,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case15,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model15=as.data.frame(ode(y=start15,
                            times=seq(0,12,1),
                            parms=params.list15,
                            fun=SIS_ode))
  
  
  ##ode(2016)
  ##5th SIA in May 2016 targeting 9m-14yrs with a cov of 95%
  cov16=MCV1_data[21,3] %>% pull("MCV1")
  cov16=cov16
  case16=Case_data[21,3] %>% pull("Cases")
  cov2_16=MCV2_data[21,3] %>% pull("MCV2")
  cov2_16=cov2_16
  
  ###Nationmode8
  SIA16<-c(rep(0,8),rep(0.95,28))
  SIA16a<-c(rep(0,8),rep(0.95,3),rep(0,25))##vaccine failure of one dose given at <1yr
  SIA16b<-c(rep(0,11),rep(0.95,7),rep(0,18))##vaccine failure of one dose given at 12-18months
  SIA16c<-c(rep(0,18),rep(0.95,18))##vaccine failure of one dose given at >18mnths
  
  start16b<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  params.list16=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov16,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_16,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case16,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model16b=as.data.frame(ode(y=start16b,
                             times=seq(0,4,1),
                             parms=params.list16,
                             fun=SIS_ode))
  ##during SIA
  myoutput16<-model16b[5,]
  startduring16<-as.data.frame(t(myoutput16[,-1])) %>% pull("5")
  ##implement SIA
  ##92% SIA,81%mcv1 and 0mcv2
  S=startduring16[1:noagegps]-((startduring16[1:noagegps]*Beta1hi*SIA16a)+
                                 (startduring16[1:noagegps]*Beta2hi*SIA16b)+
                                 (startduring16[1:noagegps]*Beta3hi*SIA16c))
  NI=startduring16[(1:noagegps)+noagegps]
  MCV1 <-startduring16[(1:noagegps)+2*noagegps]
  MCV2 <-startduring16[(1:noagegps)+3*noagegps]
  SIA <- startduring16[(1:noagegps)+4*noagegps]+
    ((startduring16[1:noagegps]*Beta1hi*SIA16a)+
       (startduring16[1:noagegps]*Beta2hi*SIA16b)+
       (startduring16[1:noagegps]*Beta3hi*SIA16c))
  startafter16<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model16a=as.data.frame(ode(y=startafter16,
                             times=seq(5,12,1), ##run from july to Aug
                             parms=params.list16,
                             fun=SIS_ode))
  
  model16<-rbind.data.frame(model16b,model16a)
  
  ##ode(2017)
  cov17=MCV1_data[22,3] %>% pull("MCV1")
  cov17=cov17
  case17=Case_data[22,3] %>% pull("Cases")
  cov2_17=MCV2_data[22,3] %>% pull("MCV2")
  cov2_17=cov2_17
  start17<-as.data.frame(t(model16[13,-1])) %>% pull("13")
  params.list17=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov17,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_17,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case17,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model17=as.data.frame(ode(y=start17,
                            times=seq(0,12,1),
                            parms=params.list17,
                            fun=SIS_ode))
  
  ##4th plot-pple in S,MCV1,MCV2,SIA,SIA2
  ##ode(2018)
  cov18=MCV1_data[23,3] %>% pull("MCV1")
  cov18=cov18
  case18=Case_data[23,3] %>% pull("Cases")
  cov2_18=MCV2_data[23,3] %>% pull("MCV2")
  cov2_18=cov2_18
  start18<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  params.list18=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov18,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_18,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case18,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model18=as.data.frame(ode(y=start18,
                            times=seq(0,12,1),
                            parms=params.list18,
                            fun=SIS_ode))
  
  ##ode(2019)
  cov19=MCV1_data[24,3] %>% pull("MCV1")
  cov19=cov19
  case19=Case_data[24,3] %>% pull("Cases")
  cov2_19=MCV2_data[24,3] %>% pull("MCV2")
  cov2_19=cov2_19
  start19<-as.data.frame(t(model18[13,-1])) %>% pull("13")
  params.list19=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov19,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_19,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case19,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model19=as.data.frame(ode(y=start19,
                            times=seq(0,12,1),
                            parms=params.list19,
                            fun=SIS_ode))
  
  
  
  ##ode(2020)
  cov20=MCV1_data[25,3] %>% pull("MCV1")
  cov20=cov20
  case20=Case_data[25,3] %>% pull("Cases")
  cov2_20=MCV2_data[25,3] %>% pull("MCV2")
  cov2_20=cov2_20
  start20<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  params.list20=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov20,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_20,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case20,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model20=as.data.frame(ode(y=start20,
                            times=seq(0,12,1),
                            parms=params.list20,
                            fun=SIS_ode))
  
  ##ode(2021)
  cov21=MCV1_data[26,3] %>% pull("MCV1")
  cov21=cov21
  case21=Case_data[26,3] %>% pull("Cases")
  cov2_21=MCV2_data[26,3] %>% pull("MCV2")
  cov2_21=cov2_21
  start21<-as.data.frame(t(model20[13,-1])) %>% pull("13")
  params.list21=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,8),cov21,rep(0,27)),
    alpha_v2=c(rep(0,17),cov2_21,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case21,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model21=as.data.frame(ode(y=start21,
                            times=seq(0,12,1),
                            parms=params.list21,
                            fun=SIS_ode))
  
  
  ###Generate the likelihood
  ##The likelihood is a combination of datasets in years;
  ##2009
  res09<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  myres09=res09 %>% matrix(nrow=noagegps) 
  rownames(myres09)=1:36
  colnames(myres09)=c("S","NI","MCV1","MCV2","SIA")
  ##2011
  res11<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  myres11=res11 %>% matrix(nrow=noagegps) 
  rownames(myres11)=1:36
  colnames(myres11)=c("S","NI","MCV1","MCV2","SIA")
  ##2013
  res13<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  myres13=res13 %>% matrix(nrow=noagegps) 
  rownames(myres13)=1:36
  colnames(myres13)=c("S","NI","MCV1","MCV2","SIA")
  ##2015
  res15<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  myres15=res15 %>% matrix(nrow=noagegps) 
  rownames(myres15)=1:36
  colnames(myres15)=c("S","NI","MCV1","MCV2","SIA")
  ##2017
  res17<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  myres17=res17 %>% matrix(nrow=noagegps) 
  rownames(myres17)=1:36
  colnames(myres17)=c("S","NI","MCV1","MCV2","SIA")
  ##2019
  res19<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  myres19=res19 %>% matrix(nrow=noagegps) 
  rownames(myres19)=1:36
  colnames(myres19)=c("S","NI","MCV1","MCV2","SIA")
  ##2021
  res21<-as.data.frame(t(model21[13,-1])) %>% pull("13")
  myres21=res21 %>% matrix(nrow=noagegps) 
  rownames(myres21)=1:36
  colnames(myres21)=c("S","NI","MCV1","MCV2","SIA")
  
  ##combined pred
  df_pred=rbind(myres09,myres11,myres13,myres15,myres17,myres19,myres21)
  df_pred=df_pred+1e-6
  
  df_RCall=as.data.frame(df_pred)
  rownames(df_RCall)=1:252
  df_RCall=df_RCall %>%
    mutate(Year=rep(seq(2009,2021,2),each=36)) %>%
    mutate(Agegroup=rep(1:36,7)) %>%
    mutate(Age_year=
             case_when(Agegroup<= 6 ~ "0.5",
                       Agegroup<= 36 ~ "0"))
  
  
  df_RCall$Age_year=factor(df_RCall$Age_year,levels=unique(df_RCall$Age_year))
  df_RCall$MI=ifelse(df_RCall$Age_year==0.5,df_RCall$S,0)
  df_RCall$S=ifelse(df_RCall$Age_year==0.5,0,df_RCall$S)
  
  df_RCall2=df_RCall %>%dplyr::select(!c(Agegroup,Age_year,Year)) %>%
    #group_by(Year) %>%
    summarise_all(sum) %>%ungroup()
  
  predc_sero.MCV1=df_RCall2[,"MCV1"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.MCV2=df_RCall2[,"MCV2"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.SIA=df_RCall2[,"SIA"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.NI=df_RCall2[,"NI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.S=df_RCall2[,"S"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.MI=df_RCall2[,"MI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  
  predc_sero.sample=cbind(predc_sero.MCV1,predc_sero.MCV2,predc_sero.SIA,predc_sero.NI)
  predc_sero.sample2=as.vector(predc_sero.sample) %>% unlist()
  
  outDf_RCall[mm,] <-  predc_sero.sample2
}

## for each row in the matrix get quantiles
quantileMatrix <- matrix(NA,nrow=ncol(outDf_RCall), ncol = 3)
for(jj in 1:ncol(outDf_RCall)){
  quantiles <- outDf_RCall[,jj] %>% quantile(probs=c(.5,.025,.975))
  quantileMatrix[jj,] <- quantiles
}

####plots for predicted and observed data
###predicted data
predRC_allgps=quantileMatrix*100
predRC_allgps=as.data.frame(predRC_allgps)
colnames(predRC_allgps)=c("midx","lox","hix")
predRC_allgps=predRC_allgps %>%
  mutate(Immunity_status=c("MCV1","MCV2","SIAs","Natural Infection"))

predRC_allgps$Immunity_status=factor(predRC_allgps$Immunity_status,
                                     levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_MCV1sens9=predRC_allgps
#Save the data
save(predRC_MCV1sens9,file = "Data/predRC_MCV1sens9")

##generate RC for MCV1=10mnths,mcv2=18mnths#########
load("Data/mainmod4")
df_chain=Nationmod9$chain[10000:50000,]
df_chain= df_chain%>% as.data.frame()
names(df_chain) =c("Beta1","Beta2","Beta3","omega")
df_chain_new=df_chain
mcmcMatrix=as.matrix(df_chain_new)
numSamples = 1000

outDf_RCall<- matrix(NA,nrow=numSamples, ncol = 4)
for (mm in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  
  Beta1hi <- mcmcMatrix[randomNumber,"Beta1"]
  Beta2hi <- mcmcMatrix[randomNumber, "Beta2"]
  Beta3hi <- mcmcMatrix[randomNumber,"Beta3"]
  omegahi <- mcmcMatrix[randomNumber,"omega"]
  
  ##########nationmode9
  cov96=MCV1_data[1,3] %>% pull("MCV1")
  cov96=cov96
  case96=Case_data[1,3] %>% pull("Cases")
  params.list96=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov96,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case96,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  
  model96=as.data.frame(ode(y=state_prop2,
                            times=seq(0,12,1),
                            parms=params.list96,
                            fun=SIS_ode))
  
  
  ##ode(1997)
  cov97=MCV1_data[2,3] %>% pull("MCV1")
  cov97=cov97
  case97=Case_data[2,3] %>% pull("Cases")
  start97<-as.data.frame(t(model96[13,-1])) %>% pull("13")
  params.list97=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov97,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case97,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model97=as.data.frame(ode(y=start97,
                            times=seq(0,12,1),
                            parms=params.list97,
                            fun=SIS_ode))
  
  
  ##ode(1998)
  cov98=MCV1_data[3,3] %>% pull("MCV1")
  case98=Case_data[3,3] %>% pull("Cases")
  cov98=cov98
  case98=Case_data[3,3] %>% pull("Cases")
  start98<-as.data.frame(t(model97[13,-1])) %>% pull("13")
  params.list98=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov98,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case98,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model98=as.data.frame(ode(y=start98,
                            times=seq(0,12,1),
                            parms=params.list98,
                            fun=SIS_ode))
  
  ##ode(1999)
  cov99=MCV1_data[4,3] %>% pull("MCV1")
  cov99=cov99
  case99=Case_data[4,3] %>% pull("Cases")
  start99<-as.data.frame(t(model98[13,-1])) %>% pull("13")
  params.list99=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov99,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case99,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model99=as.data.frame(ode(y=start99,
                            times=seq(0,12,1),
                            parms=params.list99,
                            fun=SIS_ode))
  
  ##ode(2000)
  cov00=MCV1_data[5,3] %>% pull("MCV1")
  cov00=cov00
  case00=Case_data[5,3] %>% pull("Cases")
  start00<-as.data.frame(t(model99[13,-1])) %>% pull("13")
  params.list00=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov00,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case00,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model00=as.data.frame(ode(y=start00,
                            times=seq(0,12,1),
                            parms=params.list00,
                            fun=SIS_ode))
  
  ##ode(2001)
  cov01=MCV1_data[6,3] %>% pull("MCV1")
  cov01=cov01
  case01=Case_data[6,3] %>% pull("Cases")
  start01<-as.data.frame(t(model00[13,-1])) %>% pull("13")
  params.list01=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov01,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case01,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model01=as.data.frame(ode(y=start01,
                            times=seq(0,12,1),
                            parms=params.list01,
                            fun=SIS_ode))
  
  
  ##ode(2002)
  ##1st SIA in June 2002 targetting 9m-14yrs with a cov of 94%
  cov02=MCV1_data[7,3] %>% pull("MCV1")
  cov02=cov02
  case02=Case_data[7,3] %>% pull("Cases")
  
  ###Nationmode8
  SIA02<-c(rep(0,9),rep(0.94,28))
  SIA02a<-c(rep(0,9),rep(0.94,2),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA02b<-c(rep(0,11),rep(0.94,7),rep(0,18))##vaccine failure of one dose given at 12-18months
  SIA02c<-c(rep(0,18),rep(0.94,18))##vaccine failure of one dose given at >18months
  
  
  start02b<-as.data.frame(t(model01[13,-1])) %>% pull("13")
  params.list02=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov02,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case02,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model02b=as.data.frame(ode(y=start02b,
                             times=seq(0,5,1), ##run it to end may
                             parms=params.list02,
                             fun=SIS_ode))
  
  myoutput02<-model02b[6,]
  startduring02<-as.data.frame(t(myoutput02[,-1])) %>% pull("6")
  
  ##during SIA
  ##94% SIA,68%mcv1 and 0mcv2
  
  
  S=startduring02[1:noagegps]-((startduring02[1:noagegps]*Beta1hi*SIA02a)+
                                 (startduring02[1:noagegps]*Beta2hi*SIA02b)+
                                 (startduring02[1:noagegps]*Beta3hi*SIA02c))
  NI=startduring02[(1:noagegps)+noagegps]
  MCV1 <-startduring02[(1:noagegps)+2*noagegps]
  MCV2 <-startduring02[(1:noagegps)+3*noagegps]
  SIA <- startduring02[(1:noagegps)+4*noagegps]+
    ((startduring02[1:noagegps]*Beta1hi*SIA02a)+
       (startduring02[1:noagegps]*Beta2hi*SIA02b)+
       (startduring02[1:noagegps]*Beta3hi*SIA02c))
  
  
  startafter02<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model02a=as.data.frame(ode(y=startafter02,
                             times=seq(6,12,1), ##run from july to Aug
                             parms=params.list02,
                             fun=SIS_ode))
  
  
  model02<-rbind.data.frame(model02b,model02a)
  
  ##ode(2003)
  cov03=MCV1_data[8,3] %>% pull("MCV1")
  cov03=cov03
  case03=Case_data[8,3] %>% pull("Cases")
  start03<-as.data.frame(t(model02[13,-1])) %>% pull("13")
  params.list03=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov03,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case03,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model03=as.data.frame(ode(y=start03,
                            times=seq(0,12,1),
                            parms=params.list03,
                            fun=SIS_ode))
  
  ##ode(2004)
  cov04=MCV1_data[9,3] %>% pull("MCV1")
  cov04=cov04
  case04=Case_data[9,3] %>% pull("Cases")
  start04<-as.data.frame(t(model03[13,-1])) %>% pull("13")
  params.list04=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov04,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case04,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model04=as.data.frame(ode(y=start04,
                            times=seq(0,12,1),
                            parms=params.list04,
                            fun=SIS_ode))
  
  ##ode(2005)
  cov05=MCV1_data[10,3] %>% pull("MCV1")
  cov05=cov05
  case05=Case_data[10,3] %>% pull("Cases")
  start05<-as.data.frame(t(model04[13,-1])) %>% pull("13")
  params.list05=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov05,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case05,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model05=as.data.frame(ode(y=start05,
                            times=seq(0,12,1),
                            parms=params.list05,
                            fun=SIS_ode))
  
  
  ##ode(2006)
  ##2nd SIA in JuLY 2002 targetting 9m-5yrs with a cov 90%
  cov06=MCV1_data[11,3] %>% pull("MCV1")
  cov06=cov06
  case06=Case_data[11,3] %>% pull("Cases")
  
  SIA06<-c(rep(0,9),rep(0.9,19),rep(0,9))
  SIA06a<-c(rep(0,9),rep(0.9,2),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA06b<-c(rep(0,11),rep(0.9,7),rep(0,18))##vaccine failure of one dose given at <=18months
  SIA06c<-c(rep(0,18),rep(0.9,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  start06b<-as.data.frame(t(model05[13,-1])) %>% pull("13")
  params.list06=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov06,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case06,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model06b=as.data.frame(ode(y=start06b,
                             times=seq(0,6,1),
                             parms=params.list06,
                             fun=SIS_ode))
  ##during SIA
  myoutput06<-model06b[7,]
  startduring06<-as.data.frame(t(myoutput06[,-1])) %>% pull("7")
  ##implement SIA
  ##90% SIA,67%mcv1 and 0mcv2
  
  S=startduring06[1:noagegps]-((startduring06[1:noagegps]*Beta1hi*SIA06a)+
                                 (startduring06[1:noagegps]*Beta2hi*SIA06b)+
                                 (startduring06[1:noagegps]*Beta3hi*SIA06c))
  NI=startduring06[(1:noagegps)+noagegps]
  MCV1 <-startduring06[(1:noagegps)+2*noagegps]
  MCV2 <-startduring06[(1:noagegps)+3*noagegps]
  SIA <- startduring06[(1:noagegps)+4*noagegps]+
    ((startduring06[1:noagegps]*Beta1hi*SIA06a)+
       (startduring06[1:noagegps]*Beta2hi*SIA06b)+
       (startduring06[1:noagegps]*Beta3hi*SIA06c))
  
  startafter06<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model06a=as.data.frame(ode(y=startafter06,
                             times=seq(7,12,1), ##run from july to Dec
                             parms=params.list06,
                             fun=SIS_ode))
  
  model06<-rbind.data.frame(model06b,model06a)
  
  ##ode(2007)
  cov07=MCV1_data[12,3] %>% pull("MCV1")
  cov07=cov07
  case07=Case_data[12,3] %>% pull("Cases")
  start07<-as.data.frame(t(model06[13,-1])) %>% pull("13")
  params.list07=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov07,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case07,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model07=as.data.frame(ode(y=start07,
                            times=seq(0,12,1),
                            parms=params.list07,
                            fun=SIS_ode))
  
  
  ##ode(2008)
  cov08=MCV1_data[13,3] %>% pull("MCV1")
  cov08=cov08
  case08=Case_data[13,3] %>% pull("Cases")
  case08= case08
  start08<-as.data.frame(t(model07[13,-1])) %>% pull("13")
  params.list08=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov08,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case08,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model08=as.data.frame(ode(y=start08,
                            times=seq(0,12,1),
                            parms=params.list08,
                            fun=SIS_ode))
  
  ##ode(2009)
  ##3rd SIA in September 2009 targeting 9m-5yrs with a cov 82%
  
  cov09=MCV1_data[14,3] %>% pull("MCV1")
  cov09=cov09
  case09=Case_data[14,3] %>% pull("Cases")
  case09= case09
  
  
  ###Nationmode8
  SIA09<-c(rep(0,9),rep(0.82,19),rep(0,9))
  SIA09a<-c(rep(0,9),rep(0.82,2),rep(0,25))##vaccine failure of one dose given at <11mnths
  SIA09b<-c(rep(0,11),rep(0.82,7),rep(0,18))##vaccine failure of one dose given at 10-18mnths
  SIA09c<-c(rep(0,18),rep(0.82,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start09b<-as.data.frame(t(model08[13,-1])) %>% pull("13")
  params.list09=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov09,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case09,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model09b=as.data.frame(ode(y=start09b,
                             times=seq(0,8,1),
                             parms=params.list09,
                             fun=SIS_ode))
  ##during SIA
  myoutput09<-model09b[9,]
  startduring09<-as.data.frame(t(myoutput09[,-1])) %>% pull("9")
  ##implement SIA
  ##82% SIA,77%mcv1 and 0mcv2
  
  S=startduring09[1:noagegps]-((startduring09[1:noagegps]*Beta1hi*SIA09a)+
                                 (startduring09[1:noagegps]*Beta2hi*SIA09b)+
                                 (startduring09[1:noagegps]*Beta3hi*SIA09c))
  NI=startduring09[(1:noagegps)+noagegps]
  MCV1 <-startduring09[(1:noagegps)+2*noagegps]
  MCV2 <-startduring09[(1:noagegps)+3*noagegps]
  SIA <- startduring09[(1:noagegps)+4*noagegps]+
    ((startduring09[1:noagegps]*Beta1hi*SIA09a)+
       (startduring09[1:noagegps]*Beta2hi*SIA09b)+
       (startduring09[1:noagegps]*Beta3hi*SIA09c))
  
  
  startafter09<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model09a=as.data.frame(ode(y=startafter09,
                             times=seq(9,12,1), ##run from july to Aug
                             parms=params.list09,
                             fun=SIS_ode))
  
  model09<-rbind.data.frame(model09b,model09a)
  
  
  ##ode(2010)
  cov10=MCV1_data[15,3] %>% pull("MCV1")
  cov10=cov10
  case10=Case_data[15,3] %>% pull("Cases")
  case10= case10
  start10<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  params.list10=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov10,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case10,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model10=as.data.frame(ode(y=start10,
                            times=seq(0,12,1),
                            parms=params.list10,
                            fun=SIS_ode))
  
  ##ode(2011)
  cov11=MCV1_data[16,3] %>% pull("MCV1")
  cov11=cov11
  case11=Case_data[16,3] %>% pull("Cases")
  case11=case11
  start11<-as.data.frame(t(model10[13,-1])) %>% pull("13")
  params.list11=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov11,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case11,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model11=as.data.frame(ode(y=start11,
                            times=seq(0,12,1),
                            parms=params.list11,
                            fun=SIS_ode))
  
  ##ode(2012)
  ##4th SIA in November 2012 targeting 9m-5yrs with a cov 92%
  cov12=MCV1_data[17,3] %>% pull("MCV1")
  cov12=cov12
  case12=Case_data[17,3] %>% pull("Cases")
  case12=case12
  
  ###Nationmode8
  SIA12<-c(rep(0,9),rep(0.92,19),rep(0,9))
  SIA12a<-c(rep(0,9),rep(0.92,2),rep(0,25))##vaccine failure of one dose given at <=11mnths
  SIA12b<-c(rep(0,11),rep(0.92,7),rep(0,18))##vaccine failure of one dose given at 12-18mnths
  SIA12c<-c(rep(0,18),rep(0.92,9),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start12b<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  params.list12=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov12,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case12,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model12b=as.data.frame(ode(y=start12b,
                             times=seq(0,10,1),
                             parms=params.list12,
                             fun=SIS_ode))
  ##during SIA
  myoutput12<-model12b[11,]
  startduring12<-as.data.frame(t(myoutput12[,-1])) %>% pull("11")
  ##implement SIA
  S=startduring12[1:noagegps]-((startduring12[1:noagegps]*Beta1hi*SIA12a)+
                                 (startduring12[1:noagegps]*Beta2hi*SIA12b)+
                                 (startduring12[1:noagegps]*Beta3hi*SIA12c))
  NI=startduring12[(1:noagegps)+noagegps]
  MCV1 <-startduring12[(1:noagegps)+2*noagegps]
  MCV2 <-startduring12[(1:noagegps)+3*noagegps]
  SIA <- startduring12[(1:noagegps)+4*noagegps]+
    ((startduring12[1:noagegps]*Beta1hi*SIA12a)+
       (startduring12[1:noagegps]*Beta2hi*SIA12b)+
       (startduring12[1:noagegps]*Beta3hi*SIA12c))
  
  
  startafter12<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model12a=as.data.frame(ode(y=startafter12,
                             times=seq(11,12,1), ##run from july to Aug
                             parms=params.list12,
                             fun=SIS_ode))
  
  model12<-rbind.data.frame(model12b,model12a)
  
  ##ode(2013)
  cov13=MCV1_data[18,3] %>% pull("MCV1")
  cov13=cov13
  case13=Case_data[18,3] %>% pull("Cases")
  case13=case13
  start13<-as.data.frame(t(model12[13,-1])) %>% pull("13")
  params.list13=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov13,rep(0,26)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case13,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model13=as.data.frame(ode(y=start13,
                            times=seq(0,12,1),
                            parms=params.list13,
                            fun=SIS_ode))
  
  ##ode(2014)
  ##introduce 2nd dose
  cov14=MCV1_data[19,3] %>% pull("MCV1")
  cov14=cov14
  case14=Case_data[19,3] %>% pull("Cases")
  case14=case14
  cov2_14=MCV2_data[19,3] %>% pull("MCV2")
  cov2_14=cov2_14
  start14<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  params.list14=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov14,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_14,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case14,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model14=as.data.frame(ode(y=start14,
                            times=seq(0,12,1),
                            parms=params.list14,
                            fun=SIS_ode))
  
  
  ##ode(2015)
  cov15=MCV1_data[20,3] %>% pull("MCV1")
  cov15=cov15
  case15=Case_data[20,3] %>% pull("Cases")
  cov2_15=MCV2_data[20,3] %>% pull("MCV2")
  cov2_15=cov2_15
  start15<-as.data.frame(t(model14[13,-1])) %>% pull("13")
  params.list15=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov15,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_15,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case15,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model15=as.data.frame(ode(y=start15,
                            times=seq(0,12,1),
                            parms=params.list15,
                            fun=SIS_ode))
  
  
  ##ode(2016)
  ##5th SIA in May 2016 targeting 9m-14yrs with a cov of 95%
  cov16=MCV1_data[21,3] %>% pull("MCV1")
  cov16=cov16
  case16=Case_data[21,3] %>% pull("Cases")
  cov2_16=MCV2_data[21,3] %>% pull("MCV2")
  cov2_16=cov2_16
  
  ###Nationmode8
  SIA16<-c(rep(0,9),rep(0.95,28))
  SIA16a<-c(rep(0,9),rep(0.95,2),rep(0,25))##vaccine failure of one dose given at <1yr
  SIA16b<-c(rep(0,11),rep(0.95,7),rep(0,18))##vaccine failure of one dose given at 12-18months
  SIA16c<-c(rep(0,18),rep(0.95,18))##vaccine failure of one dose given at >18mnths
  
  start16b<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  params.list16=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov16,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_16,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case16,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model16b=as.data.frame(ode(y=start16b,
                             times=seq(0,4,1),
                             parms=params.list16,
                             fun=SIS_ode))
  ##during SIA
  myoutput16<-model16b[5,]
  startduring16<-as.data.frame(t(myoutput16[,-1])) %>% pull("5")
  ##implement SIA
  ##92% SIA,81%mcv1 and 0mcv2
  S=startduring16[1:noagegps]-((startduring16[1:noagegps]*Beta1hi*SIA16a)+
                                 (startduring16[1:noagegps]*Beta2hi*SIA16b)+
                                 (startduring16[1:noagegps]*Beta3hi*SIA16c))
  NI=startduring16[(1:noagegps)+noagegps]
  MCV1 <-startduring16[(1:noagegps)+2*noagegps]
  MCV2 <-startduring16[(1:noagegps)+3*noagegps]
  SIA <- startduring16[(1:noagegps)+4*noagegps]+
    ((startduring16[1:noagegps]*Beta1hi*SIA16a)+
       (startduring16[1:noagegps]*Beta2hi*SIA16b)+
       (startduring16[1:noagegps]*Beta3hi*SIA16c))
  startafter16<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model16a=as.data.frame(ode(y=startafter16,
                             times=seq(5,12,1), ##run from july to Aug
                             parms=params.list16,
                             fun=SIS_ode))
  
  model16<-rbind.data.frame(model16b,model16a)
  
  ##ode(2017)
  cov17=MCV1_data[22,3] %>% pull("MCV1")
  cov17=cov17
  case17=Case_data[22,3] %>% pull("Cases")
  cov2_17=MCV2_data[22,3] %>% pull("MCV2")
  cov2_17=cov2_17
  start17<-as.data.frame(t(model16[13,-1])) %>% pull("13")
  params.list17=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov17,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_17,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case17,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model17=as.data.frame(ode(y=start17,
                            times=seq(0,12,1),
                            parms=params.list17,
                            fun=SIS_ode))
  
  ##4th plot-pple in S,MCV1,MCV2,SIA,SIA2
  ##ode(2018)
  cov18=MCV1_data[23,3] %>% pull("MCV1")
  cov18=cov18
  case18=Case_data[23,3] %>% pull("Cases")
  cov2_18=MCV2_data[23,3] %>% pull("MCV2")
  cov2_18=cov2_18
  start18<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  params.list18=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov18,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_18,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case18,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model18=as.data.frame(ode(y=start18,
                            times=seq(0,12,1),
                            parms=params.list18,
                            fun=SIS_ode))
  
  ##ode(2019)
  cov19=MCV1_data[24,3] %>% pull("MCV1")
  cov19=cov19
  case19=Case_data[24,3] %>% pull("Cases")
  cov2_19=MCV2_data[24,3] %>% pull("MCV2")
  cov2_19=cov2_19
  start19<-as.data.frame(t(model18[13,-1])) %>% pull("13")
  params.list19=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov19,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_19,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case19,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model19=as.data.frame(ode(y=start19,
                            times=seq(0,12,1),
                            parms=params.list19,
                            fun=SIS_ode))
  
  
  
  ##ode(2020)
  cov20=MCV1_data[25,3] %>% pull("MCV1")
  cov20=cov20
  case20=Case_data[25,3] %>% pull("Cases")
  cov2_20=MCV2_data[25,3] %>% pull("MCV2")
  cov2_20=cov2_20
  start20<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  params.list20=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov20,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_20,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case20,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model20=as.data.frame(ode(y=start20,
                            times=seq(0,12,1),
                            parms=params.list20,
                            fun=SIS_ode))
  
  ##ode(2021)
  cov21=MCV1_data[26,3] %>% pull("MCV1")
  cov21=cov21
  case21=Case_data[26,3] %>% pull("Cases")
  cov2_21=MCV2_data[26,3] %>% pull("MCV2")
  cov2_21=cov2_21
  start21<-as.data.frame(t(model20[13,-1])) %>% pull("13")
  params.list21=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,9),cov21,rep(0,26)),
    alpha_v2=c(rep(0,17),cov2_21,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case21,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model21=as.data.frame(ode(y=start21,
                            times=seq(0,12,1),
                            parms=params.list21,
                            fun=SIS_ode))
  
  
  ###Generate the likelihood
  ##The likelihood is a combination of datasets in years;
  ##2009
  res09<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  myres09=res09 %>% matrix(nrow=noagegps) 
  rownames(myres09)=1:36
  colnames(myres09)=c("S","NI","MCV1","MCV2","SIA")
  ##2011
  res11<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  myres11=res11 %>% matrix(nrow=noagegps) 
  rownames(myres11)=1:36
  colnames(myres11)=c("S","NI","MCV1","MCV2","SIA")
  ##2013
  res13<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  myres13=res13 %>% matrix(nrow=noagegps) 
  rownames(myres13)=1:36
  colnames(myres13)=c("S","NI","MCV1","MCV2","SIA")
  ##2015
  res15<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  myres15=res15 %>% matrix(nrow=noagegps) 
  rownames(myres15)=1:36
  colnames(myres15)=c("S","NI","MCV1","MCV2","SIA")
  ##2017
  res17<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  myres17=res17 %>% matrix(nrow=noagegps) 
  rownames(myres17)=1:36
  colnames(myres17)=c("S","NI","MCV1","MCV2","SIA")
  ##2019
  res19<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  myres19=res19 %>% matrix(nrow=noagegps) 
  rownames(myres19)=1:36
  colnames(myres19)=c("S","NI","MCV1","MCV2","SIA")
  ##2021
  res21<-as.data.frame(t(model21[13,-1])) %>% pull("13")
  myres21=res21 %>% matrix(nrow=noagegps) 
  rownames(myres21)=1:36
  colnames(myres21)=c("S","NI","MCV1","MCV2","SIA")
  
  ##combined pred
  df_pred=rbind(myres09,myres11,myres13,myres15,myres17,myres19,myres21)
  df_pred=df_pred+1e-6
  
  df_RCall=as.data.frame(df_pred)
  rownames(df_RCall)=1:252
  df_RCall=df_RCall %>%
    mutate(Year=rep(seq(2009,2021,2),each=36)) %>%
    mutate(Agegroup=rep(1:36,7)) %>%
    mutate(Age_year=
             case_when(Agegroup<= 6 ~ "0.5",
                       Agegroup<= 36 ~ "0"))
  
  
  df_RCall$Age_year=factor(df_RCall$Age_year,levels=unique(df_RCall$Age_year))
  df_RCall$MI=ifelse(df_RCall$Age_year==0.5,df_RCall$S,0)
  df_RCall$S=ifelse(df_RCall$Age_year==0.5,0,df_RCall$S)
  
  df_RCall2=df_RCall %>%dplyr::select(!c(Agegroup,Age_year,Year)) %>%
    #group_by(Year) %>%
    summarise_all(sum) %>%ungroup()
  
  predc_sero.MCV1=df_RCall2[,"MCV1"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.MCV2=df_RCall2[,"MCV2"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.SIA=df_RCall2[,"SIA"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.NI=df_RCall2[,"NI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.S=df_RCall2[,"S"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.MI=df_RCall2[,"MI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  
  predc_sero.sample=cbind(predc_sero.MCV1,predc_sero.MCV2,predc_sero.SIA,predc_sero.NI)
  predc_sero.sample2=as.vector(predc_sero.sample) %>% unlist()
  
  outDf_RCall[mm,] <-  predc_sero.sample2
}

## for each row in the matrix get quantiles
quantileMatrix <- matrix(NA,nrow=ncol(outDf_RCall), ncol = 3)
for(jj in 1:ncol(outDf_RCall)){
  quantiles <- outDf_RCall[,jj] %>% quantile(probs=c(.5,.025,.975))
  quantileMatrix[jj,] <- quantiles
}

####plots for predicted and observed data
###predicted data
predRC_allgps=quantileMatrix*100
predRC_allgps=as.data.frame(predRC_allgps)
colnames(predRC_allgps)=c("midx","lox","hix")
predRC_allgps=predRC_allgps %>%
  mutate(Immunity_status=c("MCV1","MCV2","SIAs","Natural Infection"))

predRC_allgps$Immunity_status=factor(predRC_allgps$Immunity_status,
                                     levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_MCV1sens10=predRC_allgps
#Save the data
save(predRC_MCV1sens10,file = "Data/predRC_MCV1sens10")

###########generate RC for MCV1=11mnths,mcv2=18 mnths,beta1<=13mnths,13=beta2<=20mnths####
load("Data/mainmod5")
df_chain=Nationmod6$chain[10000:45000,]
df_chain= df_chain%>% as.data.frame()
names(df_chain) =c("Beta1","Beta2","Beta3","omega")
df_chain_new=df_chain
mcmcMatrix=as.matrix(df_chain_new)
numSamples = 1000

outDf_RCall<- matrix(NA,nrow=numSamples, ncol = 4)
for (mm in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  
  Beta1hi <- mcmcMatrix[randomNumber,"Beta1"]
  Beta2hi <- mcmcMatrix[randomNumber, "Beta2"]
  Beta3hi <- mcmcMatrix[randomNumber,"Beta3"]
  omegahi <- mcmcMatrix[randomNumber,"omega"]
  
  cov96=MCV1_data[1,3] %>% pull("MCV1")
  cov96=cov96
  case96=Case_data[1,3] %>% pull("Cases")
  params.list96=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov96,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case96,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  
  model96=as.data.frame(ode(y=state_prop2,
                            times=seq(0,12,1),
                            parms=params.list96,
                            fun=SIS_ode))
  
  
  ##ode(1997)
  cov97=MCV1_data[2,3] %>% pull("MCV1")
  cov97=cov97
  case97=Case_data[2,3] %>% pull("Cases")
  start97<-as.data.frame(t(model96[13,-1])) %>% pull("13")
  params.list97=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov97,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case97,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model97=as.data.frame(ode(y=start97,
                            times=seq(0,12,1),
                            parms=params.list97,
                            fun=SIS_ode))
  
  
  ##ode(1998)
  cov98=MCV1_data[3,3] %>% pull("MCV1")
  case98=Case_data[3,3] %>% pull("Cases")
  cov98=cov98
  case98=Case_data[3,3] %>% pull("Cases")
  start98<-as.data.frame(t(model97[13,-1])) %>% pull("13")
  params.list98=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov98,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case98,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model98=as.data.frame(ode(y=start98,
                            times=seq(0,12,1),
                            parms=params.list98,
                            fun=SIS_ode))
  
  ##ode(1999)
  cov99=MCV1_data[4,3] %>% pull("MCV1")
  cov99=cov99
  case99=Case_data[4,3] %>% pull("Cases")
  start99<-as.data.frame(t(model98[13,-1])) %>% pull("13")
  params.list99=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov99,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case99,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model99=as.data.frame(ode(y=start99,
                            times=seq(0,12,1),
                            parms=params.list99,
                            fun=SIS_ode))
  
  ##ode(2000)
  cov00=MCV1_data[5,3] %>% pull("MCV1")
  cov00=cov00
  case00=Case_data[5,3] %>% pull("Cases")
  start00<-as.data.frame(t(model99[13,-1])) %>% pull("13")
  params.list00=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov00,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case00,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model00=as.data.frame(ode(y=start00,
                            times=seq(0,12,1),
                            parms=params.list00,
                            fun=SIS_ode))
  
  ##ode(2001)
  cov01=MCV1_data[6,3] %>% pull("MCV1")
  cov01=cov01
  case01=Case_data[6,3] %>% pull("Cases")
  start01<-as.data.frame(t(model00[13,-1])) %>% pull("13")
  params.list01=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov01,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case01,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model01=as.data.frame(ode(y=start01,
                            times=seq(0,12,1),
                            parms=params.list01,
                            fun=SIS_ode))
  
  
  ##ode(2002)
  ##1st SIA in June 2002 targetting 9m-14yrs with a cov of 94%
  cov02=MCV1_data[7,3] %>% pull("MCV1")
  cov02=cov02
  case02=Case_data[7,3] %>% pull("Cases")
  
  ##Nationmode6
  SIA02<-c(rep(0,8),rep(0.94,28))
  SIA02a<-c(rep(0,8),rep(0.94,5),rep(0,23))##vaccine failure of one dose given at <=11mnths
  SIA02b<-c(rep(0,13),rep(0.94,7),rep(0,16))##vaccine failure of one dose given at 12-18months
  SIA02c<-c(rep(0,20),rep(0.94,16))##vaccine failure of one dose given at >18months
  
  
  start02b<-as.data.frame(t(model01[13,-1])) %>% pull("13")
  params.list02=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov02,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case02,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model02b=as.data.frame(ode(y=start02b,
                             times=seq(0,5,1), ##run it to end may
                             parms=params.list02,
                             fun=SIS_ode))
  
  myoutput02<-model02b[6,]
  startduring02<-as.data.frame(t(myoutput02[,-1])) %>% pull("6")
  
  ##during SIA
  ##94% SIA,68%mcv1 and 0mcv2
  
  
  S=startduring02[1:noagegps]-((startduring02[1:noagegps]*Beta1hi*SIA02a)+
                                 (startduring02[1:noagegps]*Beta2hi*SIA02b)+
                                 (startduring02[1:noagegps]*Beta3hi*SIA02c))
  NI=startduring02[(1:noagegps)+noagegps]
  MCV1 <-startduring02[(1:noagegps)+2*noagegps]
  MCV2 <-startduring02[(1:noagegps)+3*noagegps]
  SIA <- startduring02[(1:noagegps)+4*noagegps]+
    ((startduring02[1:noagegps]*Beta1hi*SIA02a)+
       (startduring02[1:noagegps]*Beta2hi*SIA02b)+
       (startduring02[1:noagegps]*Beta3hi*SIA02c))
  
  
  startafter02<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model02a=as.data.frame(ode(y=startafter02,
                             times=seq(6,12,1), ##run from july to Aug
                             parms=params.list02,
                             fun=SIS_ode))
  
  
  model02<-rbind.data.frame(model02b,model02a)
  
  ##ode(2003)
  cov03=MCV1_data[8,3] %>% pull("MCV1")
  cov03=cov03
  case03=Case_data[8,3] %>% pull("Cases")
  start03<-as.data.frame(t(model02[13,-1])) %>% pull("13")
  params.list03=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov03,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case03,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model03=as.data.frame(ode(y=start03,
                            times=seq(0,12,1),
                            parms=params.list03,
                            fun=SIS_ode))
  
  ##ode(2004)
  cov04=MCV1_data[9,3] %>% pull("MCV1")
  cov04=cov04
  case04=Case_data[9,3] %>% pull("Cases")
  start04<-as.data.frame(t(model03[13,-1])) %>% pull("13")
  params.list04=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov04,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case04,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model04=as.data.frame(ode(y=start04,
                            times=seq(0,12,1),
                            parms=params.list04,
                            fun=SIS_ode))
  
  ##ode(2005)
  cov05=MCV1_data[10,3] %>% pull("MCV1")
  cov05=cov05
  case05=Case_data[10,3] %>% pull("Cases")
  start05<-as.data.frame(t(model04[13,-1])) %>% pull("13")
  params.list05=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov05,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case05,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model05=as.data.frame(ode(y=start05,
                            times=seq(0,12,1),
                            parms=params.list05,
                            fun=SIS_ode))
  
  
  ##ode(2006)
  ##2nd SIA in JuLY 2002 targetting 9m-5yrs with a cov 90%
  cov06=MCV1_data[11,3] %>% pull("MCV1")
  cov06=cov06
  case06=Case_data[11,3] %>% pull("Cases")
  
  SIA06<-c(rep(0,8),rep(0.9,19),rep(0,9))
  SIA06a<-c(rep(0,8),rep(0.9,5),rep(0,23))##vaccine failure of one dose given at <=11mnths
  SIA06b<-c(rep(0,13),rep(0.9,7),rep(0,16))##vaccine failure of one dose given at <=18months
  SIA06c<-c(rep(0,20),rep(0.9,7),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start06b<-as.data.frame(t(model05[13,-1])) %>% pull("13")
  params.list06=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov06,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case06,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model06b=as.data.frame(ode(y=start06b,
                             times=seq(0,6,1),
                             parms=params.list06,
                             fun=SIS_ode))
  ##during SIA
  myoutput06<-model06b[7,]
  startduring06<-as.data.frame(t(myoutput06[,-1])) %>% pull("7")
  ##implement SIA
  ##90% SIA,67%mcv1 and 0mcv2
  
  S=startduring06[1:noagegps]-((startduring06[1:noagegps]*Beta1hi*SIA06a)+
                                 (startduring06[1:noagegps]*Beta2hi*SIA06b)+
                                 (startduring06[1:noagegps]*Beta3hi*SIA06c))
  NI=startduring06[(1:noagegps)+noagegps]
  MCV1 <-startduring06[(1:noagegps)+2*noagegps]
  MCV2 <-startduring06[(1:noagegps)+3*noagegps]
  SIA <- startduring06[(1:noagegps)+4*noagegps]+
    ((startduring06[1:noagegps]*Beta1hi*SIA06a)+
       (startduring06[1:noagegps]*Beta2hi*SIA06b)+
       (startduring06[1:noagegps]*Beta3hi*SIA06c))
  
  startafter06<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model06a=as.data.frame(ode(y=startafter06,
                             times=seq(7,12,1), ##run from july to Dec
                             parms=params.list06,
                             fun=SIS_ode))
  
  model06<-rbind.data.frame(model06b,model06a)
  
  ##ode(2007)
  cov07=MCV1_data[12,3] %>% pull("MCV1")
  cov07=cov07
  case07=Case_data[12,3] %>% pull("Cases")
  start07<-as.data.frame(t(model06[13,-1])) %>% pull("13")
  params.list07=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov07,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case07,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model07=as.data.frame(ode(y=start07,
                            times=seq(0,12,1),
                            parms=params.list07,
                            fun=SIS_ode))
  
  
  ##ode(2008)
  cov08=MCV1_data[13,3] %>% pull("MCV1")
  cov08=cov08
  case08=Case_data[13,3] %>% pull("Cases")
  case08= case08
  start08<-as.data.frame(t(model07[13,-1])) %>% pull("13")
  params.list08=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov08,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case08,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model08=as.data.frame(ode(y=start08,
                            times=seq(0,12,1),
                            parms=params.list08,
                            fun=SIS_ode))
  
  ##ode(2009)
  ##3rd SIA in September 2009 targeting 9m-5yrs with a cov 82%
  
  cov09=MCV1_data[14,3] %>% pull("MCV1")
  cov09=cov09
  case09=Case_data[14,3] %>% pull("Cases")
  case09= case09
  
  
  ##Nationmode4/5
  SIA09<-c(rep(0,8),rep(0.82,19),rep(0,9))
  SIA09a<-c(rep(0,8),rep(0.82,5),rep(0,23))##vaccine failure of one dose given at <11mnths
  SIA09b<-c(rep(0,13),rep(0.82,7),rep(0,16))##vaccine failure of one dose given at 10-18mnths
  SIA09c<-c(rep(0,20),rep(0.82,7),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  
  start09b<-as.data.frame(t(model08[13,-1])) %>% pull("13")
  params.list09=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov09,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case09,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model09b=as.data.frame(ode(y=start09b,
                             times=seq(0,8,1),
                             parms=params.list09,
                             fun=SIS_ode))
  ##during SIA
  myoutput09<-model09b[9,]
  startduring09<-as.data.frame(t(myoutput09[,-1])) %>% pull("9")
  ##implement SIA
  ##82% SIA,77%mcv1 and 0mcv2
  
  S=startduring09[1:noagegps]-((startduring09[1:noagegps]*Beta1hi*SIA09a)+
                                 (startduring09[1:noagegps]*Beta2hi*SIA09b)+
                                 (startduring09[1:noagegps]*Beta3hi*SIA09c))
  NI=startduring09[(1:noagegps)+noagegps]
  MCV1 <-startduring09[(1:noagegps)+2*noagegps]
  MCV2 <-startduring09[(1:noagegps)+3*noagegps]
  SIA <- startduring09[(1:noagegps)+4*noagegps]+
    ((startduring09[1:noagegps]*Beta1hi*SIA09a)+
       (startduring09[1:noagegps]*Beta2hi*SIA09b)+
       (startduring09[1:noagegps]*Beta3hi*SIA09c))
  
  
  startafter09<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model09a=as.data.frame(ode(y=startafter09,
                             times=seq(9,12,1), ##run from july to Aug
                             parms=params.list09,
                             fun=SIS_ode))
  
  model09<-rbind.data.frame(model09b,model09a)
  
  
  ##ode(2010)
  cov10=MCV1_data[15,3] %>% pull("MCV1")
  cov10=cov10
  case10=Case_data[15,3] %>% pull("Cases")
  case10= case10
  start10<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  params.list10=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov10,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case10,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model10=as.data.frame(ode(y=start10,
                            times=seq(0,12,1),
                            parms=params.list10,
                            fun=SIS_ode))
  
  ##ode(2011)
  cov11=MCV1_data[16,3] %>% pull("MCV1")
  cov11=cov11
  case11=Case_data[16,3] %>% pull("Cases")
  case11=case11
  start11<-as.data.frame(t(model10[13,-1])) %>% pull("13")
  params.list11=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov11,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case11,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model11=as.data.frame(ode(y=start11,
                            times=seq(0,12,1),
                            parms=params.list11,
                            fun=SIS_ode))
  
  ##ode(2012)
  ##4th SIA in November 2012 targeting 9m-5yrs with a cov 92%
  cov12=MCV1_data[17,3] %>% pull("MCV1")
  cov12=cov12
  case12=Case_data[17,3] %>% pull("Cases")
  case12=case12
  
  ##Nationmode4/5
  SIA12<-c(rep(0,8),rep(0.92,19),rep(0,9))
  SIA12a<-c(rep(0,8),rep(0.92,5),rep(0,23))##vaccine failure of one dose given at <=11mnths
  SIA12b<-c(rep(0,13),rep(0.92,7),rep(0,16))##vaccine failure of one dose given at 12-18mnths
  SIA12c<-c(rep(0,20),rep(0.92,7),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start12b<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  params.list12=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov12,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case12,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model12b=as.data.frame(ode(y=start12b,
                             times=seq(0,10,1),
                             parms=params.list12,
                             fun=SIS_ode))
  ##during SIA
  myoutput12<-model12b[11,]
  startduring12<-as.data.frame(t(myoutput12[,-1])) %>% pull("11")
  ##implement SIA
  S=startduring12[1:noagegps]-((startduring12[1:noagegps]*Beta1hi*SIA12a)+
                                 (startduring12[1:noagegps]*Beta2hi*SIA12b)+
                                 (startduring12[1:noagegps]*Beta3hi*SIA12c))
  NI=startduring12[(1:noagegps)+noagegps]
  MCV1 <-startduring12[(1:noagegps)+2*noagegps]
  MCV2 <-startduring12[(1:noagegps)+3*noagegps]
  SIA <- startduring12[(1:noagegps)+4*noagegps]+
    ((startduring12[1:noagegps]*Beta1hi*SIA12a)+
       (startduring12[1:noagegps]*Beta2hi*SIA12b)+
       (startduring12[1:noagegps]*Beta3hi*SIA12c))
  
  
  startafter12<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model12a=as.data.frame(ode(y=startafter12,
                             times=seq(11,12,1), ##run from july to Aug
                             parms=params.list12,
                             fun=SIS_ode))
  
  model12<-rbind.data.frame(model12b,model12a)
  
  ##ode(2013)
  cov13=MCV1_data[18,3] %>% pull("MCV1")
  cov13=cov13
  case13=Case_data[18,3] %>% pull("Cases")
  case13=case13
  start13<-as.data.frame(t(model12[13,-1])) %>% pull("13")
  params.list13=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov13,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case13,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model13=as.data.frame(ode(y=start13,
                            times=seq(0,12,1),
                            parms=params.list13,
                            fun=SIS_ode))
  
  ##ode(2014)
  ##introduce 2nd dose
  cov14=MCV1_data[19,3] %>% pull("MCV1")
  cov14=cov14
  case14=Case_data[19,3] %>% pull("Cases")
  case14=case14
  cov2_14=MCV2_data[19,3] %>% pull("MCV2")
  cov2_14=cov2_14
  start14<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  params.list14=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov14,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_14,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case14,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model14=as.data.frame(ode(y=start14,
                            times=seq(0,12,1),
                            parms=params.list14,
                            fun=SIS_ode))
  
  
  ##ode(2015)
  cov15=MCV1_data[20,3] %>% pull("MCV1")
  cov15=cov15
  case15=Case_data[20,3] %>% pull("Cases")
  cov2_15=MCV2_data[20,3] %>% pull("MCV2")
  cov2_15=cov2_15
  start15<-as.data.frame(t(model14[13,-1])) %>% pull("13")
  params.list15=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov15,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_15,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case15,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model15=as.data.frame(ode(y=start15,
                            times=seq(0,12,1),
                            parms=params.list15,
                            fun=SIS_ode))
  
  
  ##ode(2016)
  ##5th SIA in May 2016 targeting 9m-14yrs with a cov of 95%
  cov16=MCV1_data[21,3] %>% pull("MCV1")
  cov16=cov16
  case16=Case_data[21,3] %>% pull("Cases")
  cov2_16=MCV2_data[21,3] %>% pull("MCV2")
  cov2_16=cov2_16
  
  ###Nationmode4/5
  SIA16<-c(rep(0,8),rep(0.95,28))
  SIA16a<-c(rep(0,8),rep(0.95,5),rep(0,23))##vaccine failure of one dose given at <1yr
  SIA16b<-c(rep(0,13),rep(0.95,7),rep(0,16))##vaccine failure of one dose given at 12-18months
  SIA16c<-c(rep(0,20),rep(0.95,16))##vaccine failure of one dose given at >18mnths
  
  
  start16b<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  params.list16=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov16,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_16,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case16,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model16b=as.data.frame(ode(y=start16b,
                             times=seq(0,4,1),
                             parms=params.list16,
                             fun=SIS_ode))
  ##during SIA
  myoutput16<-model16b[5,]
  startduring16<-as.data.frame(t(myoutput16[,-1])) %>% pull("5")
  ##implement SIA
  ##92% SIA,81%mcv1 and 0mcv2
  S=startduring16[1:noagegps]-((startduring16[1:noagegps]*Beta1hi*SIA16a)+
                                 (startduring16[1:noagegps]*Beta2hi*SIA16b)+
                                 (startduring16[1:noagegps]*Beta3hi*SIA16c))
  NI=startduring16[(1:noagegps)+noagegps]
  MCV1 <-startduring16[(1:noagegps)+2*noagegps]
  MCV2 <-startduring16[(1:noagegps)+3*noagegps]
  SIA <- startduring16[(1:noagegps)+4*noagegps]+
    ((startduring16[1:noagegps]*Beta1hi*SIA16a)+
       (startduring16[1:noagegps]*Beta2hi*SIA16b)+
       (startduring16[1:noagegps]*Beta3hi*SIA16c))
  startafter16<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model16a=as.data.frame(ode(y=startafter16,
                             times=seq(5,12,1), ##run from july to Aug
                             parms=params.list16,
                             fun=SIS_ode))
  
  model16<-rbind.data.frame(model16b,model16a)
  
  ##ode(2017)
  cov17=MCV1_data[22,3] %>% pull("MCV1")
  cov17=cov17
  case17=Case_data[22,3] %>% pull("Cases")
  cov2_17=MCV2_data[22,3] %>% pull("MCV2")
  cov2_17=cov2_17
  start17<-as.data.frame(t(model16[13,-1])) %>% pull("13")
  params.list17=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov17,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_17,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case17,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model17=as.data.frame(ode(y=start17,
                            times=seq(0,12,1),
                            parms=params.list17,
                            fun=SIS_ode))
  
  ##4th plot-pple in S,MCV1,MCV2,SIA,SIA2
  ##ode(2018)
  cov18=MCV1_data[23,3] %>% pull("MCV1")
  cov18=cov18
  case18=Case_data[23,3] %>% pull("Cases")
  cov2_18=MCV2_data[23,3] %>% pull("MCV2")
  cov2_18=cov2_18
  start18<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  params.list18=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov18,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_18,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case18,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model18=as.data.frame(ode(y=start18,
                            times=seq(0,12,1),
                            parms=params.list18,
                            fun=SIS_ode))
  
  ##ode(2019)
  cov19=MCV1_data[24,3] %>% pull("MCV1")
  cov19=cov19
  case19=Case_data[24,3] %>% pull("Cases")
  cov2_19=MCV2_data[24,3] %>% pull("MCV2")
  cov2_19=cov2_19
  start19<-as.data.frame(t(model18[13,-1])) %>% pull("13")
  params.list19=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov19,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_19,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case19,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model19=as.data.frame(ode(y=start19,
                            times=seq(0,12,1),
                            parms=params.list19,
                            fun=SIS_ode))
  
  
  
  ##ode(2020)
  cov20=MCV1_data[25,3] %>% pull("MCV1")
  cov20=cov20
  case20=Case_data[25,3] %>% pull("Cases")
  cov2_20=MCV2_data[25,3] %>% pull("MCV2")
  cov2_20=cov2_20
  start20<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  params.list20=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov20,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_20,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case20,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model20=as.data.frame(ode(y=start20,
                            times=seq(0,12,1),
                            parms=params.list20,
                            fun=SIS_ode))
  
  ##ode(2021)
  cov21=MCV1_data[26,3] %>% pull("MCV1")
  cov21=cov21
  case21=Case_data[26,3] %>% pull("Cases")
  cov2_21=MCV2_data[26,3] %>% pull("MCV2")
  cov2_21=cov2_21
  start21<-as.data.frame(t(model20[13,-1])) %>% pull("13")
  params.list21=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov21,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_21,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case21,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model21=as.data.frame(ode(y=start21,
                            times=seq(0,12,1),
                            parms=params.list21,
                            fun=SIS_ode))
  
  
  ###Generate the likelihood
  ##The likelihood is a combination of datasets in years;
  ##2009
  res09<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  myres09=res09 %>% matrix(nrow=noagegps) 
  rownames(myres09)=1:36
  colnames(myres09)=c("S","NI","MCV1","MCV2","SIA")
  ##2011
  res11<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  myres11=res11 %>% matrix(nrow=noagegps) 
  rownames(myres11)=1:36
  colnames(myres11)=c("S","NI","MCV1","MCV2","SIA")
  ##2013
  res13<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  myres13=res13 %>% matrix(nrow=noagegps) 
  rownames(myres13)=1:36
  colnames(myres13)=c("S","NI","MCV1","MCV2","SIA")
  ##2015
  res15<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  myres15=res15 %>% matrix(nrow=noagegps) 
  rownames(myres15)=1:36
  colnames(myres15)=c("S","NI","MCV1","MCV2","SIA")
  ##2017
  res17<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  myres17=res17 %>% matrix(nrow=noagegps) 
  rownames(myres17)=1:36
  colnames(myres17)=c("S","NI","MCV1","MCV2","SIA")
  ##2019
  res19<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  myres19=res19 %>% matrix(nrow=noagegps) 
  rownames(myres19)=1:36
  colnames(myres19)=c("S","NI","MCV1","MCV2","SIA")
  ##2021
  res21<-as.data.frame(t(model21[13,-1])) %>% pull("13")
  myres21=res21 %>% matrix(nrow=noagegps) 
  rownames(myres21)=1:36
  colnames(myres21)=c("S","NI","MCV1","MCV2","SIA")
  
  ##combined pred
  df_pred=rbind(myres09,myres11,myres13,myres15,myres17,myres19,myres21)
  df_pred=df_pred+1e-6
  
  df_RCall=as.data.frame(df_pred)
  rownames(df_RCall)=1:252
  df_RCall=df_RCall %>%
    mutate(Year=rep(seq(2009,2021,2),each=36)) %>%
    mutate(Agegroup=rep(1:36,7)) %>%
    mutate(Age_year=
             case_when(Agegroup<= 6 ~ "0.5",
                       Agegroup<= 36 ~ "0"))
  
  
  df_RCall$Age_year=factor(df_RCall$Age_year,levels=unique(df_RCall$Age_year))
  df_RCall$MI=ifelse(df_RCall$Age_year==0.5,df_RCall$S,0)
  df_RCall$S=ifelse(df_RCall$Age_year==0.5,0,df_RCall$S)
  
  df_RCall2=df_RCall %>%dplyr::select(!c(Agegroup,Age_year,Year)) %>%
    #group_by(Year) %>%
    summarise_all(sum) %>%ungroup()
  
  predc_sero.MCV1=df_RCall2[,"MCV1"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.MCV2=df_RCall2[,"MCV2"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.SIA=df_RCall2[,"SIA"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.NI=df_RCall2[,"NI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.S=df_RCall2[,"S"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.MI=df_RCall2[,"MI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  
  predc_sero.sample=cbind(predc_sero.MCV1,predc_sero.MCV2,predc_sero.SIA,predc_sero.NI)
  predc_sero.sample2=as.vector(predc_sero.sample) %>% unlist()
  
  outDf_RCall[mm,] <-  predc_sero.sample2
}

## for each row in the matrix get quantiles
quantileMatrix <- matrix(NA,nrow=ncol(outDf_RCall), ncol = 3)
for(jj in 1:ncol(outDf_RCall)){
  quantiles <- outDf_RCall[,jj] %>% quantile(probs=c(.5,.025,.975))
  quantileMatrix[jj,] <- quantiles
}

####plots for predicted and observed data
###predicted data
predRC_allgps=quantileMatrix*100
predRC_allgps=as.data.frame(predRC_allgps)
colnames(predRC_allgps)=c("midx","lox","hix")
predRC_allgps=predRC_allgps %>%
  mutate(Immunity_status=c("MCV1","MCV2","SIAs","Natural Infection"))

predRC_allgps$Immunity_status=factor(predRC_allgps$Immunity_status,
                                     levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_betasens2=predRC_allgps
#Save the data
save(predRC_betasens2,file = "Data/predRC_betasens2")


###########generate RC for MCV1=11mnths,mcv2=18 mnths,beta1<=15mnths,13=beta2<=22mnths####
load("Data/mainmod6")
df_chain=Nationmod7$chain[10000:45000,]
df_chain= df_chain%>% as.data.frame()
names(df_chain) =c("Beta1","Beta2","Beta3","omega")
df_chain_new=df_chain
mcmcMatrix=as.matrix(df_chain_new)
numSamples = 1000

outDf_RCall<- matrix(NA,nrow=numSamples, ncol = 4)
for (mm in 1:numSamples ) {
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  
  Beta1hi <- mcmcMatrix[randomNumber,"Beta1"]
  Beta2hi <- mcmcMatrix[randomNumber, "Beta2"]
  Beta3hi <- mcmcMatrix[randomNumber,"Beta3"]
  omegahi <- mcmcMatrix[randomNumber,"omega"]
  
  cov96=MCV1_data[1,3] %>% pull("MCV1")
  cov96=cov96
  case96=Case_data[1,3] %>% pull("Cases")
  params.list96=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov96,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case96,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  
  model96=as.data.frame(ode(y=state_prop2,
                            times=seq(0,12,1),
                            parms=params.list96,
                            fun=SIS_ode))
  
  
  ##ode(1997)
  cov97=MCV1_data[2,3] %>% pull("MCV1")
  cov97=cov97
  case97=Case_data[2,3] %>% pull("Cases")
  start97<-as.data.frame(t(model96[13,-1])) %>% pull("13")
  params.list97=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov97,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case97,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model97=as.data.frame(ode(y=start97,
                            times=seq(0,12,1),
                            parms=params.list97,
                            fun=SIS_ode))
  
  
  ##ode(1998)
  cov98=MCV1_data[3,3] %>% pull("MCV1")
  case98=Case_data[3,3] %>% pull("Cases")
  cov98=cov98
  case98=Case_data[3,3] %>% pull("Cases")
  start98<-as.data.frame(t(model97[13,-1])) %>% pull("13")
  params.list98=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov98,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case98,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model98=as.data.frame(ode(y=start98,
                            times=seq(0,12,1),
                            parms=params.list98,
                            fun=SIS_ode))
  
  ##ode(1999)
  cov99=MCV1_data[4,3] %>% pull("MCV1")
  cov99=cov99
  case99=Case_data[4,3] %>% pull("Cases")
  start99<-as.data.frame(t(model98[13,-1])) %>% pull("13")
  params.list99=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov99,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case99,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model99=as.data.frame(ode(y=start99,
                            times=seq(0,12,1),
                            parms=params.list99,
                            fun=SIS_ode))
  
  ##ode(2000)
  cov00=MCV1_data[5,3] %>% pull("MCV1")
  cov00=cov00
  case00=Case_data[5,3] %>% pull("Cases")
  start00<-as.data.frame(t(model99[13,-1])) %>% pull("13")
  params.list00=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov00,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case00,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model00=as.data.frame(ode(y=start00,
                            times=seq(0,12,1),
                            parms=params.list00,
                            fun=SIS_ode))
  
  ##ode(2001)
  cov01=MCV1_data[6,3] %>% pull("MCV1")
  cov01=cov01
  case01=Case_data[6,3] %>% pull("Cases")
  start01<-as.data.frame(t(model00[13,-1])) %>% pull("13")
  params.list01=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov01,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case01,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model01=as.data.frame(ode(y=start01,
                            times=seq(0,12,1),
                            parms=params.list01,
                            fun=SIS_ode))
  
  
  ##ode(2002)
  ##1st SIA in June 2002 targetting 9m-14yrs with a cov of 94%
  cov02=MCV1_data[7,3] %>% pull("MCV1")
  cov02=cov02
  case02=Case_data[7,3] %>% pull("Cases")
  
  ##Nationmode6
  SIA02<-c(rep(0,8),rep(0.94,28))
  SIA02a<-c(rep(0,8),rep(0.94,7),rep(0,21))##vaccine failure of one dose given at <=11mnths
  SIA02b<-c(rep(0,15),rep(0.94,7),rep(0,14))##vaccine failure of one dose given at 12-18months
  SIA02c<-c(rep(0,22),rep(0.94,14))##vaccine failure of one dose given at >18months
  
  
  start02b<-as.data.frame(t(model01[13,-1])) %>% pull("13")
  params.list02=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov02,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case02,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  
  model02b=as.data.frame(ode(y=start02b,
                             times=seq(0,5,1), ##run it to end may
                             parms=params.list02,
                             fun=SIS_ode))
  
  myoutput02<-model02b[6,]
  startduring02<-as.data.frame(t(myoutput02[,-1])) %>% pull("6")
  
  ##during SIA
  ##94% SIA,68%mcv1 and 0mcv2
  
  
  S=startduring02[1:noagegps]-((startduring02[1:noagegps]*Beta1hi*SIA02a)+
                                 (startduring02[1:noagegps]*Beta2hi*SIA02b)+
                                 (startduring02[1:noagegps]*Beta3hi*SIA02c))
  NI=startduring02[(1:noagegps)+noagegps]
  MCV1 <-startduring02[(1:noagegps)+2*noagegps]
  MCV2 <-startduring02[(1:noagegps)+3*noagegps]
  SIA <- startduring02[(1:noagegps)+4*noagegps]+
    ((startduring02[1:noagegps]*Beta1hi*SIA02a)+
       (startduring02[1:noagegps]*Beta2hi*SIA02b)+
       (startduring02[1:noagegps]*Beta3hi*SIA02c))
  
  
  startafter02<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model02a=as.data.frame(ode(y=startafter02,
                             times=seq(6,12,1), ##run from july to Aug
                             parms=params.list02,
                             fun=SIS_ode))
  
  
  model02<-rbind.data.frame(model02b,model02a)
  
  ##ode(2003)
  cov03=MCV1_data[8,3] %>% pull("MCV1")
  cov03=cov03
  case03=Case_data[8,3] %>% pull("Cases")
  start03<-as.data.frame(t(model02[13,-1])) %>% pull("13")
  params.list03=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov03,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case03,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model03=as.data.frame(ode(y=start03,
                            times=seq(0,12,1),
                            parms=params.list03,
                            fun=SIS_ode))
  
  ##ode(2004)
  cov04=MCV1_data[9,3] %>% pull("MCV1")
  cov04=cov04
  case04=Case_data[9,3] %>% pull("Cases")
  start04<-as.data.frame(t(model03[13,-1])) %>% pull("13")
  params.list04=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov04,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case04,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model04=as.data.frame(ode(y=start04,
                            times=seq(0,12,1),
                            parms=params.list04,
                            fun=SIS_ode))
  
  ##ode(2005)
  cov05=MCV1_data[10,3] %>% pull("MCV1")
  cov05=cov05
  case05=Case_data[10,3] %>% pull("Cases")
  start05<-as.data.frame(t(model04[13,-1])) %>% pull("13")
  params.list05=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov05,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case05,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model05=as.data.frame(ode(y=start05,
                            times=seq(0,12,1),
                            parms=params.list05,
                            fun=SIS_ode))
  
  
  ##ode(2006)
  ##2nd SIA in JuLY 2002 targetting 9m-5yrs with a cov 90%
  cov06=MCV1_data[11,3] %>% pull("MCV1")
  cov06=cov06
  case06=Case_data[11,3] %>% pull("Cases")
  
  SIA06<-c(rep(0,8),rep(0.9,19),rep(0,9))
  SIA06a<-c(rep(0,8),rep(0.9,7),rep(0,21))##vaccine failure of one dose given at <=11mnths
  SIA06b<-c(rep(0,15),rep(0.9,7),rep(0,14))##vaccine failure of one dose given at <=18months
  SIA06c<-c(rep(0,22),rep(0.9,5),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start06b<-as.data.frame(t(model05[13,-1])) %>% pull("13")
  params.list06=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov06,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case06,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model06b=as.data.frame(ode(y=start06b,
                             times=seq(0,6,1),
                             parms=params.list06,
                             fun=SIS_ode))
  ##during SIA
  myoutput06<-model06b[7,]
  startduring06<-as.data.frame(t(myoutput06[,-1])) %>% pull("7")
  ##implement SIA
  ##90% SIA,67%mcv1 and 0mcv2
  
  S=startduring06[1:noagegps]-((startduring06[1:noagegps]*Beta1hi*SIA06a)+
                                 (startduring06[1:noagegps]*Beta2hi*SIA06b)+
                                 (startduring06[1:noagegps]*Beta3hi*SIA06c))
  NI=startduring06[(1:noagegps)+noagegps]
  MCV1 <-startduring06[(1:noagegps)+2*noagegps]
  MCV2 <-startduring06[(1:noagegps)+3*noagegps]
  SIA <- startduring06[(1:noagegps)+4*noagegps]+
    ((startduring06[1:noagegps]*Beta1hi*SIA06a)+
       (startduring06[1:noagegps]*Beta2hi*SIA06b)+
       (startduring06[1:noagegps]*Beta3hi*SIA06c))
  
  startafter06<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model06a=as.data.frame(ode(y=startafter06,
                             times=seq(7,12,1), ##run from july to Dec
                             parms=params.list06,
                             fun=SIS_ode))
  
  model06<-rbind.data.frame(model06b,model06a)
  
  ##ode(2007)
  cov07=MCV1_data[12,3] %>% pull("MCV1")
  cov07=cov07
  case07=Case_data[12,3] %>% pull("Cases")
  start07<-as.data.frame(t(model06[13,-1])) %>% pull("13")
  params.list07=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov07,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case07,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model07=as.data.frame(ode(y=start07,
                            times=seq(0,12,1),
                            parms=params.list07,
                            fun=SIS_ode))
  
  
  ##ode(2008)
  cov08=MCV1_data[13,3] %>% pull("MCV1")
  cov08=cov08
  case08=Case_data[13,3] %>% pull("Cases")
  case08= case08
  start08<-as.data.frame(t(model07[13,-1])) %>% pull("13")
  params.list08=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov08,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case08,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model08=as.data.frame(ode(y=start08,
                            times=seq(0,12,1),
                            parms=params.list08,
                            fun=SIS_ode))
  
  ##ode(2009)
  ##3rd SIA in September 2009 targeting 9m-5yrs with a cov 82%
  
  cov09=MCV1_data[14,3] %>% pull("MCV1")
  cov09=cov09
  case09=Case_data[14,3] %>% pull("Cases")
  case09= case09
  
  
  ##Nationmode4/5
  SIA09<-c(rep(0,8),rep(0.82,19),rep(0,9))
  SIA09a<-c(rep(0,8),rep(0.82,7),rep(0,21))##vaccine failure of one dose given at <11mnths
  SIA09b<-c(rep(0,15),rep(0.82,7),rep(0,14))##vaccine failure of one dose given at 10-18mnths
  SIA09c<-c(rep(0,22),rep(0.82,5),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start09b<-as.data.frame(t(model08[13,-1])) %>% pull("13")
  params.list09=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov09,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case09,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model09b=as.data.frame(ode(y=start09b,
                             times=seq(0,8,1),
                             parms=params.list09,
                             fun=SIS_ode))
  ##during SIA
  myoutput09<-model09b[9,]
  startduring09<-as.data.frame(t(myoutput09[,-1])) %>% pull("9")
  ##implement SIA
  ##82% SIA,77%mcv1 and 0mcv2
  
  S=startduring09[1:noagegps]-((startduring09[1:noagegps]*Beta1hi*SIA09a)+
                                 (startduring09[1:noagegps]*Beta2hi*SIA09b)+
                                 (startduring09[1:noagegps]*Beta3hi*SIA09c))
  NI=startduring09[(1:noagegps)+noagegps]
  MCV1 <-startduring09[(1:noagegps)+2*noagegps]
  MCV2 <-startduring09[(1:noagegps)+3*noagegps]
  SIA <- startduring09[(1:noagegps)+4*noagegps]+
    ((startduring09[1:noagegps]*Beta1hi*SIA09a)+
       (startduring09[1:noagegps]*Beta2hi*SIA09b)+
       (startduring09[1:noagegps]*Beta3hi*SIA09c))
  
  
  startafter09<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model09a=as.data.frame(ode(y=startafter09,
                             times=seq(9,12,1), ##run from july to Aug
                             parms=params.list09,
                             fun=SIS_ode))
  
  model09<-rbind.data.frame(model09b,model09a)
  
  
  ##ode(2010)
  cov10=MCV1_data[15,3] %>% pull("MCV1")
  cov10=cov10
  case10=Case_data[15,3] %>% pull("Cases")
  case10= case10
  start10<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  params.list10=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov10,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case10,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model10=as.data.frame(ode(y=start10,
                            times=seq(0,12,1),
                            parms=params.list10,
                            fun=SIS_ode))
  
  ##ode(2011)
  cov11=MCV1_data[16,3] %>% pull("MCV1")
  cov11=cov11
  case11=Case_data[16,3] %>% pull("Cases")
  case11=case11
  start11<-as.data.frame(t(model10[13,-1])) %>% pull("13")
  params.list11=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov11,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case11,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model11=as.data.frame(ode(y=start11,
                            times=seq(0,12,1),
                            parms=params.list11,
                            fun=SIS_ode))
  
  ##ode(2012)
  ##4th SIA in November 2012 targeting 9m-5yrs with a cov 92%
  cov12=MCV1_data[17,3] %>% pull("MCV1")
  cov12=cov12
  case12=Case_data[17,3] %>% pull("Cases")
  case12=case12
  
  ##Nationmode4/5
  SIA12<-c(rep(0,8),rep(0.92,19),rep(0,9))
  SIA12a<-c(rep(0,8),rep(0.92,7),rep(0,21))##vaccine failure of one dose given at <=11mnths
  SIA12b<-c(rep(0,15),rep(0.92,7),rep(0,14))##vaccine failure of one dose given at 12-18mnths
  SIA12c<-c(rep(0,22),rep(0.92,5),rep(0,9))##vaccine failure of one dose given at >18mnths
  
  
  start12b<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  params.list12=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov12,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case12,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model12b=as.data.frame(ode(y=start12b,
                             times=seq(0,10,1),
                             parms=params.list12,
                             fun=SIS_ode))
  ##during SIA
  myoutput12<-model12b[11,]
  startduring12<-as.data.frame(t(myoutput12[,-1])) %>% pull("11")
  ##implement SIA
  S=startduring12[1:noagegps]-((startduring12[1:noagegps]*Beta1hi*SIA12a)+
                                 (startduring12[1:noagegps]*Beta2hi*SIA12b)+
                                 (startduring12[1:noagegps]*Beta3hi*SIA12c))
  NI=startduring12[(1:noagegps)+noagegps]
  MCV1 <-startduring12[(1:noagegps)+2*noagegps]
  MCV2 <-startduring12[(1:noagegps)+3*noagegps]
  SIA <- startduring12[(1:noagegps)+4*noagegps]+
    ((startduring12[1:noagegps]*Beta1hi*SIA12a)+
       (startduring12[1:noagegps]*Beta2hi*SIA12b)+
       (startduring12[1:noagegps]*Beta3hi*SIA12c))
  
  
  startafter12<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model12a=as.data.frame(ode(y=startafter12,
                             times=seq(11,12,1), ##run from july to Aug
                             parms=params.list12,
                             fun=SIS_ode))
  
  model12<-rbind.data.frame(model12b,model12a)
  
  ##ode(2013)
  cov13=MCV1_data[18,3] %>% pull("MCV1")
  cov13=cov13
  case13=Case_data[18,3] %>% pull("Cases")
  case13=case13
  start13<-as.data.frame(t(model12[13,-1])) %>% pull("13")
  params.list13=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov13,rep(0,25)),
    alpha_v2=c(rep(0,36)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case13,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  model13=as.data.frame(ode(y=start13,
                            times=seq(0,12,1),
                            parms=params.list13,
                            fun=SIS_ode))
  
  ##ode(2014)
  ##introduce 2nd dose
  cov14=MCV1_data[19,3] %>% pull("MCV1")
  cov14=cov14
  case14=Case_data[19,3] %>% pull("Cases")
  case14=case14
  cov2_14=MCV2_data[19,3] %>% pull("MCV2")
  cov2_14=cov2_14
  start14<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  params.list14=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov14,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_14,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case14,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model14=as.data.frame(ode(y=start14,
                            times=seq(0,12,1),
                            parms=params.list14,
                            fun=SIS_ode))
  
  
  ##ode(2015)
  cov15=MCV1_data[20,3] %>% pull("MCV1")
  cov15=cov15
  case15=Case_data[20,3] %>% pull("Cases")
  cov2_15=MCV2_data[20,3] %>% pull("MCV2")
  cov2_15=cov2_15
  start15<-as.data.frame(t(model14[13,-1])) %>% pull("13")
  params.list15=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov15,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_15,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case15,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model15=as.data.frame(ode(y=start15,
                            times=seq(0,12,1),
                            parms=params.list15,
                            fun=SIS_ode))
  
  
  ##ode(2016)
  ##5th SIA in May 2016 targeting 9m-14yrs with a cov of 95%
  cov16=MCV1_data[21,3] %>% pull("MCV1")
  cov16=cov16
  case16=Case_data[21,3] %>% pull("Cases")
  cov2_16=MCV2_data[21,3] %>% pull("MCV2")
  cov2_16=cov2_16
  
  ###Nationmode4/5
  SIA16<-c(rep(0,8),rep(0.95,28))
  SIA16a<-c(rep(0,8),rep(0.95,7),rep(0,21))##vaccine failure of one dose given at <1yr
  SIA16b<-c(rep(0,15),rep(0.95,7),rep(0,14))##vaccine failure of one dose given at 12-18months
  SIA16c<-c(rep(0,22),rep(0.95,14))##vaccine failure of one dose given at >18mnths
  
  start16b<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  params.list16=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov16,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_16,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case16,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model16b=as.data.frame(ode(y=start16b,
                             times=seq(0,4,1),
                             parms=params.list16,
                             fun=SIS_ode))
  ##during SIA
  myoutput16<-model16b[5,]
  startduring16<-as.data.frame(t(myoutput16[,-1])) %>% pull("5")
  ##implement SIA
  ##92% SIA,81%mcv1 and 0mcv2
  S=startduring16[1:noagegps]-((startduring16[1:noagegps]*Beta1hi*SIA16a)+
                                 (startduring16[1:noagegps]*Beta2hi*SIA16b)+
                                 (startduring16[1:noagegps]*Beta3hi*SIA16c))
  NI=startduring16[(1:noagegps)+noagegps]
  MCV1 <-startduring16[(1:noagegps)+2*noagegps]
  MCV2 <-startduring16[(1:noagegps)+3*noagegps]
  SIA <- startduring16[(1:noagegps)+4*noagegps]+
    ((startduring16[1:noagegps]*Beta1hi*SIA16a)+
       (startduring16[1:noagegps]*Beta2hi*SIA16b)+
       (startduring16[1:noagegps]*Beta3hi*SIA16c))
  startafter16<-c(S,NI,MCV1,MCV2,SIA)
  
  ##continue running till end year
  model16a=as.data.frame(ode(y=startafter16,
                             times=seq(5,12,1), ##run from july to Aug
                             parms=params.list16,
                             fun=SIS_ode))
  
  model16<-rbind.data.frame(model16b,model16a)
  
  ##ode(2017)
  cov17=MCV1_data[22,3] %>% pull("MCV1")
  cov17=cov17
  case17=Case_data[22,3] %>% pull("Cases")
  cov2_17=MCV2_data[22,3] %>% pull("MCV2")
  cov2_17=cov2_17
  start17<-as.data.frame(t(model16[13,-1])) %>% pull("13")
  params.list17=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov17,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_17,rep(0,18)),   
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case17,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model17=as.data.frame(ode(y=start17,
                            times=seq(0,12,1),
                            parms=params.list17,
                            fun=SIS_ode))
  
  ##4th plot-pple in S,MCV1,MCV2,SIA,SIA2
  ##ode(2018)
  cov18=MCV1_data[23,3] %>% pull("MCV1")
  cov18=cov18
  case18=Case_data[23,3] %>% pull("Cases")
  cov2_18=MCV2_data[23,3] %>% pull("MCV2")
  cov2_18=cov2_18
  start18<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  params.list18=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov18,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_18,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case18,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model18=as.data.frame(ode(y=start18,
                            times=seq(0,12,1),
                            parms=params.list18,
                            fun=SIS_ode))
  
  ##ode(2019)
  cov19=MCV1_data[24,3] %>% pull("MCV1")
  cov19=cov19
  case19=Case_data[24,3] %>% pull("Cases")
  cov2_19=MCV2_data[24,3] %>% pull("MCV2")
  cov2_19=cov2_19
  start19<-as.data.frame(t(model18[13,-1])) %>% pull("13")
  params.list19=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov19,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_19,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case19,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model19=as.data.frame(ode(y=start19,
                            times=seq(0,12,1),
                            parms=params.list19,
                            fun=SIS_ode))
  
  
  
  ##ode(2020)
  cov20=MCV1_data[25,3] %>% pull("MCV1")
  cov20=cov20
  case20=Case_data[25,3] %>% pull("Cases")
  cov2_20=MCV2_data[25,3] %>% pull("MCV2")
  cov2_20=cov2_20
  start20<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  params.list20=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov20,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_20,rep(0,18)),    
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case20,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model20=as.data.frame(ode(y=start20,
                            times=seq(0,12,1),
                            parms=params.list20,
                            fun=SIS_ode))
  
  ##ode(2021)
  cov21=MCV1_data[26,3] %>% pull("MCV1")
  cov21=cov21
  case21=Case_data[26,3] %>% pull("Cases")
  cov2_21=MCV2_data[26,3] %>% pull("MCV2")
  cov2_21=cov2_21
  start21<-as.data.frame(t(model20[13,-1])) %>% pull("13")
  params.list21=list(
    noagegps=noagegps,
    alpha_v1=c(rep(0,10),cov21,rep(0,25)),
    alpha_v2=c(rep(0,17),cov2_21,rep(0,18)),  
    alpha_v3=c(rep(0,36)),
    ageout=ageing,
    casedf=c(rep(0,6),rep(case21,30)),
    omega=omegahi,
    beta1=Beta1hi,
    beta2=Beta2hi,
    beta3=Beta3hi,
    agein=c(1/sum(1/ageing) , ageing[-noagegps]))
  
  model21=as.data.frame(ode(y=start21,
                            times=seq(0,12,1),
                            parms=params.list21,
                            fun=SIS_ode))
  
  
  ###Generate the likelihood
  ##The likelihood is a combination of datasets in years;
  ##2009
  res09<-as.data.frame(t(model09[13,-1])) %>% pull("13")
  myres09=res09 %>% matrix(nrow=noagegps) 
  rownames(myres09)=1:36
  colnames(myres09)=c("S","NI","MCV1","MCV2","SIA")
  ##2011
  res11<-as.data.frame(t(model11[13,-1])) %>% pull("13")
  myres11=res11 %>% matrix(nrow=noagegps) 
  rownames(myres11)=1:36
  colnames(myres11)=c("S","NI","MCV1","MCV2","SIA")
  ##2013
  res13<-as.data.frame(t(model13[13,-1])) %>% pull("13")
  myres13=res13 %>% matrix(nrow=noagegps) 
  rownames(myres13)=1:36
  colnames(myres13)=c("S","NI","MCV1","MCV2","SIA")
  ##2015
  res15<-as.data.frame(t(model15[13,-1])) %>% pull("13")
  myres15=res15 %>% matrix(nrow=noagegps) 
  rownames(myres15)=1:36
  colnames(myres15)=c("S","NI","MCV1","MCV2","SIA")
  ##2017
  res17<-as.data.frame(t(model17[13,-1])) %>% pull("13")
  myres17=res17 %>% matrix(nrow=noagegps) 
  rownames(myres17)=1:36
  colnames(myres17)=c("S","NI","MCV1","MCV2","SIA")
  ##2019
  res19<-as.data.frame(t(model19[13,-1])) %>% pull("13")
  myres19=res19 %>% matrix(nrow=noagegps) 
  rownames(myres19)=1:36
  colnames(myres19)=c("S","NI","MCV1","MCV2","SIA")
  ##2021
  res21<-as.data.frame(t(model21[13,-1])) %>% pull("13")
  myres21=res21 %>% matrix(nrow=noagegps) 
  rownames(myres21)=1:36
  colnames(myres21)=c("S","NI","MCV1","MCV2","SIA")
  
  ##combined pred
  df_pred=rbind(myres09,myres11,myres13,myres15,myres17,myres19,myres21)
  df_pred=df_pred+1e-6
  
  df_RCall=as.data.frame(df_pred)
  rownames(df_RCall)=1:252
  df_RCall=df_RCall %>%
    mutate(Year=rep(seq(2009,2021,2),each=36)) %>%
    mutate(Agegroup=rep(1:36,7)) %>%
    mutate(Age_year=
             case_when(Agegroup<= 6 ~ "0.5",
                       Agegroup<= 36 ~ "0"))
  
  
  df_RCall$Age_year=factor(df_RCall$Age_year,levels=unique(df_RCall$Age_year))
  df_RCall$MI=ifelse(df_RCall$Age_year==0.5,df_RCall$S,0)
  df_RCall$S=ifelse(df_RCall$Age_year==0.5,0,df_RCall$S)
  
  df_RCall2=df_RCall %>%dplyr::select(!c(Agegroup,Age_year,Year)) %>%
    #group_by(Year) %>%
    summarise_all(sum) %>%ungroup()
  
  predc_sero.MCV1=df_RCall2[,"MCV1"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.MCV2=df_RCall2[,"MCV2"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.SIA=df_RCall2[,"SIA"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  predc_sero.NI=df_RCall2[,"NI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.S=df_RCall2[,"S"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  #predc_sero.MI=df_RCall2[,"MI"]/rowSums(df_RCall2[,c("MCV1","MCV2","SIA","NI")])
  
  predc_sero.sample=cbind(predc_sero.MCV1,predc_sero.MCV2,predc_sero.SIA,predc_sero.NI)
  predc_sero.sample2=as.vector(predc_sero.sample) %>% unlist()
  
  outDf_RCall[mm,] <-  predc_sero.sample2
}

## for each row in the matrix get quantiles
quantileMatrix <- matrix(NA,nrow=ncol(outDf_RCall), ncol = 3)
for(jj in 1:ncol(outDf_RCall)){
  quantiles <- outDf_RCall[,jj] %>% quantile(probs=c(.5,.025,.975))
  quantileMatrix[jj,] <- quantiles
}

####plots for predicted and observed data
###predicted data
predRC_allgps=quantileMatrix*100
predRC_allgps=as.data.frame(predRC_allgps)
colnames(predRC_allgps)=c("midx","lox","hix")
predRC_allgps=predRC_allgps %>%
  mutate(Immunity_status=c("MCV1","MCV2","SIAs","Natural Infection"))

predRC_allgps$Immunity_status=factor(predRC_allgps$Immunity_status,
                                     levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_betasens4=predRC_allgps
#Save the data
save(predRC_betasens4,file = "Data/predRC_betasens4")




