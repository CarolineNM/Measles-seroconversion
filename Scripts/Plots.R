rm(list=ls())
library(tidyverse)
library(ggpubr)
#########plot observed and predicted figure 2a main text###########
load("Data/obs_6gps")
df_obs6<-obs_6gps%>% mutate(Data=rep("Observed",42))
names(df_obs6)=c("Year","Age_cat","midxb","loxb","hixb","Data")
load("Data/dfbase2_fit6gps")
dfbase2_fit6gps=dfbase2_fit6gps
df2_pred6b=dfbase2_fit6gps%>% mutate(Data=rep("Predicted",42))
combined_fit=rbind.data.frame(df_obs6,df2_pred6b)
combined_fit$Age_cat=factor(combined_fit$Age_cat,levels=c("<9m","9m-<1yr",
                                                          "1-<2yrs","2-4yrs",
                                                          "5-9yrs","10-14yrs"))

combined_fit$Data=factor(combined_fit$Data,levels=unique(combined_fit$Data))
colors <- c("Observed"="black","Predicted"="royalblue1")
baseline=ggplot(combined_fit, aes(x=Age_cat,y=midxb*100, ymin=loxb*100, ymax=hixb*100,color=Data)) +
  geom_pointrange(size=0.7,position=position_dodge(0.9))+
  facet_grid(~Year)+
  theme_bw() +
  scale_color_manual(values = colors)+
  labs(y = "Measles seroprevalence",
       x = "Age groups",
       color = "Legend")+
  theme(axis.text.x = element_text(size=15,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        strip.text.y = element_text(size=20,face="plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=20,face="plain"),
        legend.position="bottom")

########plot overall RC baseline main figure 2b################
load("Data/predRC_allbase")
dfbaseRC=as.data.frame(predRC_allbase) 
dfbaseRC$Immunity_status=factor(dfbaseRC$Immunity_status,
                                levels=c("MCV1","Natural Infection","SIAs","MCV2"))


RC_year=dfbaseRC %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,fill=Immunity_status))+
  geom_bar(stat="identity",position=position_dodge(0.9),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  
  theme_bw()+
  theme(axis.text.x = element_text(size=20,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        legend.title = element_text(size=20,face="plain"),
        legend.position="none",
        legend.text=element_text(size=20))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution among all seroconversions (%)")+
  coord_flip()
RC_year

##########merge the plots together to  generate figure2###########
tiff("Figures/mainfigure2.tiff", units="in", width=13, height=12, res=300)
ggarrange(baseline,RC_year,
          labels = c("a", "b"),
          ncol = 1, nrow=2)

dev.off()


########plot yearly RC baseline figure s5################
load("Data/predRC_yearlybase")
predRC_yearlybase$Year=factor(predRC_yearlybase$Year,levels=unique(predRC_yearlybase$Year))
predRC_yearlybase$Immunity_status=factor(predRC_yearlybase$Immunity_status,
                                         levels=c("MCV1","Natural Infection","SIAs","MCV2"))

tiff("Figures/figures5.tiff", units="in", width=8, height=4, res=300)
predRC_yearlybase %>% ggplot2::ggplot(aes(x=Year,y=midx,fill=Immunity_status))+
  geom_bar(stat="identity",position=position_stack(0.9),width=0.75)+
  #geom_errorbar(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix,colour=Immunity_status), position=position_dodge(0.75),alpha=1, size=0.4)+
  theme_bw()+
  theme(axis.text.x = element_text(size=16,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=16,face="plain"),
        axis.title.x=element_text(size=16,face="plain"),
        axis.title.y=element_text(size=16,face="plain"),
        strip.text.x = element_text(size=16,face="plain"),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=16))+
  xlab("Survey Year")+
  ylab("Relative contribution (%)")+
  coord_flip()
dev.off()

############figure s6  projections fit##############
#########Combine data and have one graph figure s6######################
load("Data/dfbase2_fit6gps")
dfbase2_fit6gps=dfbase2_fit6gps
df2_pred6b=dfbase2_fit6gps%>% mutate(Data=rep("Baseline",42))
load("Data/dfmcv1_fit6gps")
dfmcv1a=dfmcv1_fit6gps%>% mutate(Data=rep("Increased MCV1 and MCV2",42))
load("Data/dfmcv1b_fit6gps")
dfmcv1b=dfmcv1b_fit6gps%>% mutate(Data=rep("Increased and timely MCV1 and MCV2",42))
load("Data/dfmcv1c_fit6gps")
dfmcv1c=dfmcv1c_fit6gps%>% mutate(Data=rep("Ideal and timely MCV1 and MCV2",42))
combined_fit=rbind.data.frame(dfmcv1a,dfmcv1b,dfmcv1c)
combined_fit=combined_fit %>% mutate(Type=rep("MCV1+MCV2",126))
combined_fit$Age_cat=factor(combined_fit$Age_cat,levels=c("<9m","9m-<1yr",
                                                          "1-<2yrs","2-4yrs",
                                                          "5-9yrs","10-14yrs"))


load("Data/dfmcv2b_fit6gps")
dfmcv2b=dfmcv2b_fit6gps%>% mutate(Data=rep("High MCV2",42))
load("Data/dfmcv2a_fit6gps")
dfmcv2a=dfmcv2a_fit6gps%>% mutate(Data=rep("Higher MCV2",42))
combined_fit2=rbind.data.frame(df2_pred6b,dfmcv2b,dfmcv2a)
combined_fit2=combined_fit2 %>% mutate(Type=rep("MCV2 only",126))


######merge combined_fit2 and combined_fit
combined_dat=rbind.data.frame(combined_fit2,combined_fit)
combined_dat$Data=factor(combined_dat$Data,levels=unique(combined_dat$Data))
colors <- c("Baseline"="black","Higher MCV2"="orange","High MCV2"="blue","Ideal and timely MCV1 and MCV2"="red","Increased and timely MCV1 and MCV2"="green","Increased MCV1 and MCV2"="brown")
tiff("Figures/figures6.tiff", units="in", width=18, height=13, res=300)
ggplot(combined_dat, aes(x=Age_cat,y=midxb*100, ymin=loxb*100, ymax=hixb*100,color=Data)) +
  geom_pointrange(size=0.5,position=position_dodge(0.9))+
  #facet_grid(Type~Year)+
  facet_grid(~Year)+
  theme_bw() +
  scale_color_manual(values = colors)+
  labs(y = "Measles seroprevalence",
       x = "Age groups",
       color = "MCV2 projections")+
  theme(axis.text.x = element_text(size=18,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        strip.text.y = element_text(size=20,face="plain"),
        legend.title = element_blank(),
        legend.text=element_text(size=18,face="plain"),
        legend.position="bottom")

dev.off()
###########figure s8 RC contribution among entire pop############
load("Data/predRC_allbaseold")
dfbaseRC=predRC_allbase %>% mutate(Data=rep("Baseline",6))
load("Data/predRC_mcv2a")
dfmcv2aRC=predRC_mcv2a %>% mutate(Data=rep("Higher MCV2",6))
load("Data/predRC_mcv2b")
dfmcv2bRC=predRC_mcv2b %>% mutate(Data=rep("High MCV2",6)) 
load("Data/predRC_mcv1")
dfmcv1aRC=predRC_mcv1 %>% mutate(Data=rep("Increased MCV1 and MCV2",6)) 
load("Data/predRC_mcv1b")
dfmcv1bRC=predRC_mcv1b %>% mutate(Data=rep("Increased and timely MCV1 and MCV2",6)) 
load("Data/predRC_mcv1c")
dfmcv1cRC=predRC_mcv1c %>% mutate(Data=rep("Ideal and timely MCV1 and MCV2",6)) 

combined_RCb=rbind.data.frame(dfbaseRC,dfmcv2bRC,dfmcv2aRC,dfmcv1aRC,dfmcv1bRC,dfmcv1cRC)
combined_RCb$Data=factor(combined_RCb$Data,levels=unique(combined_RCb$Data))
combined_RCb$Immunity_status=factor(combined_RCb$Immunity_status,levels=c("Immune by MCV1","Immune by Natural Infection",
                                                                          "Immune by SIAs","Susceptible","Protected by Maternal Immunity",
                                                                          "Immune by MCV2"))

combined_RCb=combined_RCb %>% mutate(Immunity_status=recode(Immunity_status,"Immune by MCV1"="MCV1",
                                                            "Immune by Natural Infection"="Natural Infection",
                                                            "Immune by SIAs"="SIAs",
                                                            "Susceptible"="Susceptible",
                                                            "Protected by Maternal Immunity"="Maternal Immunity",
                                                            "Immune by MCV2"="MCV2"))

combined_RCb=combined_RCb %>% filter(Immunity_status!="Maternal Immunity")

colors <- c("Baseline"="black","Higher MCV2"="orange","High MCV2"="blue","Ideal and timely MCV1 and MCV2"="red","Increased and timely MCV1 and MCV2"="green","Increased MCV1 and MCV2"="brown")
tiff("Figures/figures8.tiff", units="in", width=16, height=9, res=300)
combined_RCb %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,fill=Data))+
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  theme_bw()+
  #facet_grid(~Type)+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(size=20,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        legend.title = element_blank(),
        legend.position="bottom",
        #legend.justification = "left",
        #legend.direction = "vertical",
        legend.text=element_text(size=17,face="plain"))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution among all population (%)")+
  coord_flip()
dev.off()


#####################RC of scenarios in the years after MCV2 introduction 2015-2021 main figure 3####
load("Data/RC_newbase3")
dfbaseRC=RC_newbase3 %>% mutate(Data=rep("Baseline",5))
load("Data/RC_mcv2aonly")
dfmcv2aRC=RC_mcv2aonly %>% mutate(Data=rep("Higher MCV2",5))
load("Data/RC_mcv2bonly")
dfmcv2bRC=RC_mcv2bonly %>% mutate(Data=rep("High MCV2",5)) 
load("Data/RC_mcv1mcv2a")
dfmcv1aRC=RC_mcv1mcv2a %>% mutate(Data=rep("Increased MCV1 and MCV2",5)) 
load("Data/RC_mcv1mcv2b")
dfmcv1bRC=RC_mcv1mcv2b %>% mutate(Data=rep("Increased and timely MCV1 and MCV2",5)) 
load("Data/RC_mcv1mcv2c")
dfmcv1cRC=RC_mcv1mcv2c %>% mutate(Data=rep("Ideal and timely MCV1 and MCV2",5)) 

combined_RCb=rbind.data.frame(dfbaseRC,dfmcv2bRC,dfmcv2aRC,dfmcv1aRC,dfmcv1bRC,dfmcv1cRC)
combined_RCb$Data=factor(combined_RCb$Data,levels=unique(combined_RCb$Data))
combined_RCb$Immunity_status=factor(combined_RCb$Immunity_status,levels=c("MCV1","Natural Infection",
                                                                          "SIAs","Susceptible","MCV2"))


colors <- c("Baseline"="black","Higher MCV2"="orange","High MCV2"="blue","Ideal and timely MCV1 and MCV2"="red","Increased and timely MCV1 and MCV2"="green","Increased MCV1 and MCV2"="brown")

tiff("Figures/mainfigure3.tiff", units="in", width=15, height=11, res=300)
combined_RCb %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,fill=Data))+
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  theme_bw()+
  #facet_grid(~Type)+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(size=20,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        legend.title = element_blank(),
        legend.position="bottom",
        #legend.justification = "left",
        #legend.direction = "vertical",
        legend.text=element_text(size=18,face="plain"))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution in the entire population(%)")+
  coord_flip()

dev.off()


###############New RC baseline for seroconversion only in years after MCV2 supp figure s6######
load("Data/RC_newbase2")
dfbaseRC=as.data.frame(RC_newbase2) 
dfbaseRC$Immunity_status=factor(dfbaseRC$Immunity_status,
                                levels=c("MCV1","Natural Infection","SIAs","MCV2"))

tiff("Figures/figures6RCbaseMCV2yrs.tiff", units="in", width=10, height=5, res=300)
dfbaseRC %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,fill=Immunity_status))+
  geom_bar(stat="identity",position=position_dodge(0.9),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  
  theme_bw()+
  theme(axis.text.x = element_text(size=16,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=16,face="plain"),
        axis.title.x=element_text(size=16,face="plain"),
        axis.title.y=element_text(size=16,face="plain"),
        strip.text.x = element_text(size=16,face="plain"),
        legend.title = element_text(size=16,face="plain"),
        legend.position="none",
        legend.text=element_text(size=16))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution among all seroconversions (%)")+
  coord_flip()

dev.off()

################RC in the sensitivity analysis for delayed timeliness of MCV1##########
load("DATA_STEF_REV/predRC_allbase") ####baseline
dfbaseRC=as.data.frame(predRC_allbase) 
dfbaseRC$Immunity_status=factor(dfbaseRC$Immunity_status,
                                levels=c("MCV1","Natural Infection","SIAs","MCV2"))

dfbaseRC=dfbaseRC %>% mutate(Type=rep("Baseline",4))

load("DATA_STEF_REV/predRC_MCV1sens9") 
dfMCV1sens9RC=as.data.frame(predRC_MCV1sens9) 
dfMCV1sens9RC$Immunity_status=factor(dfMCV1sens9RC$Immunity_status,
                                     levels=c("MCV1","Natural Infection","SIAs","MCV2"))

dfMCV1sens9RC=dfMCV1sens9RC %>% mutate(Type=rep("MCV1 at 9months",4))


load("DATA_STEF_REV/predRC_MCV1sens10") 
predRC_MCV1sens10=as.data.frame(predRC_MCV1sens10) 
predRC_MCV1sens10$Immunity_status=factor(predRC_MCV1sens10$Immunity_status,
                                         levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_MCV1sens10=predRC_MCV1sens10 %>% mutate(Type=rep("MCV1 at 10months",4))

newRC=rbind.data.frame(dfbaseRC,dfMCV1sens9RC,predRC_MCV1sens10)
tiff("Final-figures/figures9sensitivity.tiff", units="in", width=13, height=7, res=300)
newRC %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,group=Type,fill=Type))+
  geom_bar(stat="identity",position=position_dodge(0.9),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  #facet_grid(~Type)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        legend.title = element_text(size=20,face="plain"),
        legend.position="bottom",
        legend.text=element_text(size=20))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution among all seroconversions (%)")+
  coord_flip()
dev.off()


################RC in the sensitivity analysis for beta age cutoff##########
load("DATA_STEF_REV/predRC_allbase") ####baseline
dfbaseRC=as.data.frame(predRC_allbase) 
dfbaseRC$Immunity_status=factor(dfbaseRC$Immunity_status,
                                levels=c("MCV1","Natural Infection","SIAs","MCV2"))

dfbaseRC=dfbaseRC %>% mutate(Type=rep("Baseline",4))

load("DATA_STEF_REV/predRC_betasens2") 
predRC_betasens2=as.data.frame(predRC_betasens2) 
predRC_betasens2$Immunity_status=factor(predRC_betasens2$Immunity_status,
                                        levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_betasens2=predRC_betasens2 %>% mutate(Type=rep("2 months increase in Vaccine failure age cutoff",4))


load("DATA_STEF_REV/predRC_betasens4") 
predRC_betasens4=as.data.frame(predRC_betasens4) 
predRC_betasens4$Immunity_status=factor(predRC_betasens4$Immunity_status,
                                        levels=c("MCV1","Natural Infection","SIAs","MCV2"))

predRC_betasens4=predRC_betasens4 %>% mutate(Type=rep("4 months increase in Vaccine failure age cutoff",4))

newRC=rbind.data.frame(dfbaseRC,predRC_betasens2,predRC_betasens4)

tiff("Final-figures/figures10betasensitivity.tiff", units="in", width=13, height=7, res=300)
newRC %>% ggplot2::ggplot(aes(x=Immunity_status,y=midx,group=Type,fill=Type))+
  geom_bar(stat="identity",position=position_dodge(0.9),width=0.75)+
  geom_linerange(mapping=aes(x=Immunity_status,y=midx,ymin=lox,ymax=hix), position=position_dodge(0.75), size=1)+
  #facet_grid(~Type)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20,face="plain",angle=60,hjust=1),
        axis.text.y = element_text(size=20,face="plain"),
        axis.title.x=element_text(size=20,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        strip.text.x = element_text(size=20,face="plain"),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=13))+
  xlab("Modelled seroconversion pathway")+
  ylab("Relative contribution among all seroconversions (%)")+
  coord_flip()
dev.off()





