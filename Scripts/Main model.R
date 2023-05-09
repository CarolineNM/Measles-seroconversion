#####baseline model###
##MCV1 11mnths,MCV2 18months
##beta1<12mnths,12=beta2<=18mnths,beta3>18mnths
rm(list=ls())
library(deSolve)
library(tidyverse)
library(rootSolve)
library(readxl)
library(MASS)
library(MCMCvis)
library(coda)

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
av_cases2<-sum(Case_data$Cases)/nrow(Case_data)
##########Observed seroprevalence###
load("Data/obs_sero")
df_obs<-obs_sero

#####likelihood starts here because we have introduced uncertainty in coverage
startvalue=c(0.9,0.9,0.7,0.04)###beta1,beta2,beta3,Omega
d=length(startvalue)
iterations=50000
chain=array(dim=c(iterations+1,d))
chain[1,]<-startvalue
p_curr <- startvalue
names(p_curr) <- c("beta1","beta2","beta3","omega") 
scaling=1          
scaling_max=10
cov0= diag(p_curr)*c(1e-2,1e-2,1e-2,1e-2) 
cov = cov0
LL_curr=-10000 
acceptance=0
time<-system.time({    
  for(i in 2:iterations){
    p_prop<-0.05*mvrnorm(1,(p_curr),scaling*(2.38^2)/d*cov0)+(1-0.05)*mvrnorm(1,(p_curr),scaling*(2.38^2)/d*cov)##How did this come about?
    
    if(any(p_prop<0)|p_prop["beta1"]>1|p_prop["beta2"]>1|
       p_prop["beta3"]>1|p_prop["omega"]>1){
      LL_prop=-10000
      prior_LL_prop=-10000
      prior_LL_curr=-10000
      
    } else{
      
      ##ode(1996)
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
      
    
      SIA02<-c(rep(0,8),rep(0.94,28))
      SIA02a<-c(rep(0,8),rep(0.94,3),rep(0,25))##vaccine failure of one dose given at <12mnths
      SIA02b<-c(rep(0,11),rep(0.94,7),rep(0,18))##vaccine failure of one dose given at 12-18months
      SIA02c<-c(rep(0,18),rep(0.94,18))##vaccine failure of one dose given at >18months
      
      
      start02b<-as.data.frame(t(model01[13,-1])) %>% pull("13")
      params.list02=list(
        noagegps=noagegps,
        alpha_v1=c(rep(0,10),cov02,rep(0,25)),
        alpha_v2=c(rep(0,36)),  
        alpha_v3=c(rep(0,36)),
        ageout=ageing,
        casedf=c(rep(0,6),rep(case02,30)),
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      
      model02b=as.data.frame(ode(y=start02b,
                                 times=seq(0,5,1), ##run it to end may
                                 parms=params.list02,
                                 fun=SIS_ode))
      
      myoutput02<-model02b[6,]
      startduring02<-as.data.frame(t(myoutput02[,-1])) %>% pull("6")
      
      ##during SIA
      
      
      S=startduring02[1:noagegps]-((startduring02[1:noagegps]*p_prop["beta1"]*SIA02a)+
                                     (startduring02[1:noagegps]*p_prop["beta2"]*SIA02b)+
                                     (startduring02[1:noagegps]*p_prop["beta3"]*SIA02c))
      NI=startduring02[(1:noagegps)+noagegps]
      MCV1 <-startduring02[(1:noagegps)+2*noagegps]
      MCV2 <-startduring02[(1:noagegps)+3*noagegps]
      SIA <- startduring02[(1:noagegps)+4*noagegps]+
        ((startduring02[1:noagegps]*p_prop["beta1"]*SIA02a)+
           (startduring02[1:noagegps]*p_prop["beta2"]*SIA02b)+
           (startduring02[1:noagegps]*p_prop["beta3"]*SIA02c))
      
      
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
      SIA06a<-c(rep(0,8),rep(0.9,3),rep(0,25))##vaccine failure of one dose given at <12mnths
      SIA06b<-c(rep(0,11),rep(0.9,7),rep(0,18))##vaccine failure of one dose given at 12-18months
      SIA06c<-c(rep(0,18),rep(0.9,9),rep(0,9))##vaccine failure of one dose given at >18months
      

      start06b<-as.data.frame(t(model05[13,-1])) %>% pull("13")
      params.list06=list(
        noagegps=noagegps,
        alpha_v1=c(rep(0,10),cov06,rep(0,25)),
        alpha_v2=c(rep(0,36)),  
        alpha_v3=c(rep(0,36)),
        ageout=ageing,
        casedf=c(rep(0,6),rep(case06,30)),
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      model06b=as.data.frame(ode(y=start06b,
                                 times=seq(0,6,1),
                                 parms=params.list06,
                                 fun=SIS_ode))
      ##during SIA
      myoutput06<-model06b[7,]
      startduring06<-as.data.frame(t(myoutput06[,-1])) %>% pull("7")
      ##implement SIA
      
      
      S=startduring06[1:noagegps]-((startduring06[1:noagegps]*p_prop["beta1"]*SIA06a)+
                                     (startduring06[1:noagegps]*p_prop["beta2"]*SIA06b)+
                                     (startduring06[1:noagegps]*p_prop["beta3"]*SIA06c))
      NI=startduring06[(1:noagegps)+noagegps]
      MCV1 <-startduring06[(1:noagegps)+2*noagegps]
      MCV2 <-startduring06[(1:noagegps)+3*noagegps]
      SIA <- startduring06[(1:noagegps)+4*noagegps]+
        ((startduring06[1:noagegps]*p_prop["beta1"]*SIA06a)+
           (startduring06[1:noagegps]*p_prop["beta2"]*SIA06b)+
           (startduring06[1:noagegps]*p_prop["beta3"]*SIA06c))
      
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        
     
      SIA09<-c(rep(0,8),rep(0.82,19),rep(0,9))
      SIA09a<-c(rep(0,8),rep(0.82,3),rep(0,25))##vaccine failure of one dose given at <12mnths
      SIA09b<-c(rep(0,11),rep(0.82,7),rep(0,18))##vaccine failure of one dose given at 12-18months
      SIA09c<-c(rep(0,18),rep(0.82,9),rep(0,9))##vaccine failure of one dose given at >18months
      
       
      start09b<-as.data.frame(t(model08[13,-1])) %>% pull("13")
      params.list09=list(
        noagegps=noagegps,
        alpha_v1=c(rep(0,10),cov09,rep(0,25)),
        alpha_v2=c(rep(0,36)),  
        alpha_v3=c(rep(0,36)),
        ageout=ageing,
        casedf=c(rep(0,6),rep(case09,30)),
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      model09b=as.data.frame(ode(y=start09b,
                                 times=seq(0,8,1),
                                 parms=params.list09,
                                 fun=SIS_ode))
      ##during SIA
      myoutput09<-model09b[9,]
      startduring09<-as.data.frame(t(myoutput09[,-1])) %>% pull("9")
      ##implement SIA
      
      S=startduring09[1:noagegps]-((startduring09[1:noagegps]*p_prop["beta1"]*SIA09a)+
                                     (startduring09[1:noagegps]*p_prop["beta2"]*SIA09b)+
                                     (startduring09[1:noagegps]*p_prop["beta3"]*SIA09c))
      NI=startduring09[(1:noagegps)+noagegps]
      MCV1 <-startduring09[(1:noagegps)+2*noagegps]
      MCV2 <-startduring09[(1:noagegps)+3*noagegps]
      SIA <- startduring09[(1:noagegps)+4*noagegps]+
        ((startduring09[1:noagegps]*p_prop["beta1"]*SIA09a)+
           (startduring09[1:noagegps]*p_prop["beta2"]*SIA09b)+
           (startduring09[1:noagegps]*p_prop["beta3"]*SIA09c))
      
      
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
         
      SIA12<-c(rep(0,8),rep(0.92,19),rep(0,9))
      SIA12a<-c(rep(0,8),rep(0.92,3),rep(0,25))##vaccine failure of one dose given at <12mnths
      SIA12b<-c(rep(0,11),rep(0.92,7),rep(0,18))##vaccine failure of one dose given at 12-18months
      SIA12c<-c(rep(0,18),rep(0.92,9),rep(0,9))##vaccine failure of one dose given at >18months
      
        
      start12b<-as.data.frame(t(model11[13,-1])) %>% pull("13")
      params.list12=list(
        noagegps=noagegps,
        alpha_v1=c(rep(0,10),cov12,rep(0,25)),
        alpha_v2=c(rep(0,36)),  
        alpha_v3=c(rep(0,36)),
        ageout=ageing,
        casedf=c(rep(0,6),rep(case12,30)),
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      model12b=as.data.frame(ode(y=start12b,
                                 times=seq(0,10,1),
                                 parms=params.list12,
                                 fun=SIS_ode))
      ##during SIA
      myoutput12<-model12b[11,]
      startduring12<-as.data.frame(t(myoutput12[,-1])) %>% pull("11")
      ##implement SIA
      S=startduring12[1:noagegps]-((startduring12[1:noagegps]*p_prop["beta1"]*SIA12a)+
                                     (startduring12[1:noagegps]*p_prop["beta2"]*SIA12b)+
                                     (startduring12[1:noagegps]*p_prop["beta3"]*SIA12c))
      NI=startduring12[(1:noagegps)+noagegps]
      MCV1 <-startduring12[(1:noagegps)+2*noagegps]
      MCV2 <-startduring12[(1:noagegps)+3*noagegps]
      SIA <- startduring12[(1:noagegps)+4*noagegps]+
        ((startduring12[1:noagegps]*p_prop["beta1"]*SIA12a)+
           (startduring12[1:noagegps]*p_prop["beta2"]*SIA12b)+
           (startduring12[1:noagegps]*p_prop["beta3"]*SIA12c))
      
      
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
      
      SIA16<-c(rep(0,8),rep(0.95,28))
      SIA16a<-c(rep(0,8),rep(0.95,3),rep(0,25))##vaccine failure of one dose given at <12mnths
      SIA16b<-c(rep(0,11),rep(0.95,7),rep(0,18))##vaccine failure of one dose given at 12-18months
      SIA16c<-c(rep(0,18),rep(0.95,18))##vaccine failure of one dose given at >18months
        
      start16b<-as.data.frame(t(model15[13,-1])) %>% pull("13")
      params.list16=list(
        noagegps=noagegps,
        alpha_v1=c(rep(0,10),cov16,rep(0,25)),
        alpha_v2=c(rep(0,17),cov2_16,rep(0,18)),    
        alpha_v3=c(rep(0,36)),
        ageout=ageing,
        casedf=c(rep(0,6),rep(case16,30)),
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      model16b=as.data.frame(ode(y=start16b,
                                 times=seq(0,4,1),
                                 parms=params.list16,
                                 fun=SIS_ode))
      ##during SIA
      myoutput16<-model16b[5,]
      startduring16<-as.data.frame(t(myoutput16[,-1])) %>% pull("5")
      ##implement SIA
      
      S=startduring16[1:noagegps]-((startduring16[1:noagegps]*p_prop["beta1"]*SIA16a)+
                                     (startduring16[1:noagegps]*p_prop["beta2"]*SIA16b)+
                                     (startduring16[1:noagegps]*p_prop["beta3"]*SIA16c))
      NI=startduring16[(1:noagegps)+noagegps]
      MCV1 <-startduring16[(1:noagegps)+2*noagegps]
      MCV2 <-startduring16[(1:noagegps)+3*noagegps]
      SIA <- startduring16[(1:noagegps)+4*noagegps]+
        ((startduring16[1:noagegps]*p_prop["beta1"]*SIA16a)+
           (startduring16[1:noagegps]*p_prop["beta2"]*SIA16b)+
           (startduring16[1:noagegps]*p_prop["beta3"]*SIA16c))
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
        agein=c(1/sum(1/ageing) , ageing[-noagegps]))
      
      model17=as.data.frame(ode(y=start17,
                                times=seq(0,12,1),
                                parms=params.list17,
                                fun=SIS_ode))
      
      
      
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
        omega=p_prop["omega"],
        beta1=p_prop["beta1"],
        beta2=p_prop["beta2"],
        beta3=p_prop["beta3"],
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
      
      predc_sero=(rowSums(df_pred[,c("NI","MCV1","MCV2","SIA")])/rowSums(df_pred))
      
      LL_prop = likelihood(pred_sero= predc_sero) 
      prior_LL_prop=prior(param=p_prop)
      prior_LL_curr=prior(param=p_curr)
      
      if (i%%100==0){
        print(predc_sero)
      }
    }    
    
    LL_ratio <- LL_prop-LL_curr +  prior_LL_prop - prior_LL_curr
    
    if(log(runif(1)) < LL_ratio){
      chain[i,]= p_prop
      p_curr = p_prop
      LL_curr = LL_prop
      acceptance = acceptance + 1
    }else{
      chain[i,]= p_curr
    }
    
    # adaptive MCMC 
    if (i>200){ # start updating after 200 iterations
      cov <- cov(chain[1:i,])
      scaling <- exp(log(scaling) + (0.9999)^(i)*(acceptance/(i) - 0.234))
      scaling = min(scaling,scaling_max)
    }
    
    if(i%%2==0) print(round(c(i,acceptance/i,LL_curr,LL_prop,scaling),2))
    
    mainmod <- list(chain=chain)
    save(mainmod ,file = "Data/mainmod")
    
    
  }}) # end for loop
time 

###run mainmod2 with different startvalues#######
########plot the chains######
####check chain convergence using gelman rubin taht uses difference in variance between and within the 2 chains
###A scale reduction factor below 1.01 is okay
source("Scripts/Convergencechains.R")
load("Data/mainmod")
chain_conv(burnreps =1000, sfix="chain_1")
chain_conv2(burnreps =1000, sfix="chain_1")
df_chaina=mainmod$chain[10000:50000,]
df_chainb=mainmod2$chain[10000:50000,]
mcmc_chain=as.mcmc(df_chain)
varnames(mcmc_chain)=c("Beta1","Beta2","Beta3","Omega")
mcmc_chainb=as.mcmc(df_chainb)
varnames(mcmc_chain)=c("Beta1","Beta2","Beta3","Omega")
combinedchains = mcmc.list(mcmc_chain,mcmc_chainb)
plot(combinedchains)
data=MCMCsummary(combinedchains) ## Check ESS and Rhat








