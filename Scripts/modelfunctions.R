library(rriskDistributions)
#######vaccine efficacy priors#####
beta1<-0.84[0.72-0.95]
beta2<-0.93[0.85-0.97]
beta3<-0.93[0.85-0.97]
##generate the parameters of the beta distributions
mybeta1<-get.beta.par(p = c(0.025, 0.5, 0.975), q=c(0.72,0.84,0.95))
mybeta2<-get.beta.par(p = c(0.025, 0.5, 0.975), q=c(0.85,0.93,0.97))


######Binomial likelihood####
likelihood<-function(dt=df_obs,pred_sero){
  if(any(pred_sero <=0)|any(pred_sero >=1)){
    lik=-10000
  }else{
    lik=sum(dt$Pos*log(pred_sero)+(dt$Total-dt$Pos)*log(1-pred_sero))
  }
  return(lik)
}

######priors
prior<-function(param){
  beta1prior = dbeta(param["beta1"],shape1=37.36535,shape2=7.37432 )
  beta2prior= dbeta(param["beta2"],shape1=77.991905, shape2=6.189579)
  beta3prior= dbeta(param["beta3"],shape1=77.991905, shape2=6.189579)
  omegaprior=  dunif(param["omega"],min=0,max=1)
  return(log(beta1prior*beta2prior*beta3prior*omegaprior))
}


SIS_ode <- function(times, state, parms) {
  # # parameters
  noagegps=parms$noagegps
  alpha_v1=parms$alpha_v1
  alpha_v2=parms$alpha_v2
  alpha_v3=parms$alpha_v3
  age.out=parms$ageout
  age.in=parms$agein
  omega=parms$omega
  beta1=parms$beta1
  beta2=parms$beta2
  beta3=parms$beta3
  casedf=parms$casedf
  
  # states
  S <- state[1:noagegps]
  NI <- state[(1:noagegps)+noagegps]
  MCV1 <- state[(1:noagegps)+2*noagegps]
  MCV2 <- state[(1:noagegps)+3*noagegps]
  SIA <- state[(1:noagegps)+4*noagegps]
  
  N<-S+NI+MCV1+MCV2+SIA
  
  
  ##equations incorporating aging and vaccination coverage
  lambda=omega*(casedf/av_cases) 
  
  dS<- (1-(lambda+(beta1*alpha_v1)+(beta2*alpha_v2)+(beta1*alpha_v3)+
             (beta2*alpha_v3)+(beta3*alpha_v3)))*age.in*c(sum(N),S[-noagegps])-age.out*S
  
  dNI<-lambda*age.in*c(sum(N),S[-noagegps])+age.in*c(0,NI[-noagegps])-age.out*NI
  dMCV1<- (beta1*alpha_v1)*age.in*c(sum(N),S[-noagegps])+age.in*c(0,MCV1[-noagegps])-age.out*MCV1
  dMCV2<- (beta2*alpha_v2)*age.in*c(sum(N),S[-noagegps])+age.in*c(0,MCV2[-noagegps])-age.out*MCV2
  dSIA<-((beta1*alpha_v3)+(beta2*alpha_v3)+(beta3*alpha_v3))*age.in*c(sum(N),S[-noagegps])+
    age.in*c(0,SIA[-noagegps])-age.out*SIA
  
  return(list(c(dS,dNI,dMCV1,dMCV2,dSIA))) 
}
