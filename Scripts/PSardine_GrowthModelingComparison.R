rm(list = ls())
# Growth model comparisons for Pacific Sardine
# Four data sets: 
# 1: original dataset 2012-2021
# 2: Dorval et al. 2015 dataset, 2004-2010
# 3: Dorval et al. 2015 with age 0 from Dataset 1 (add rather than replace)

# load packages ####
library(FSA)
library(nlme)
library(mgcv)
library(car)
library(AICcmodavg)
library(lubridate)
library(minpack.lm)
library(nlstools)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(tidyverse)
library(viridis)
library(ggpubr)

#import data ####
# original dataset 
#original <- read.xlsx("Data/sard.u.xlsx")
original <- read.xlsx("Data/sard_u_raw.xlsx")
original$Model <- "Current Study"

original$sl <- original$standardLength_mm
original$fage <- as.factor(original$age)


# Dorval et al. 2015 dataset 2004-2010
ps.0410 <- read.xlsx("Data/sard0410_u.xlsx")
ps.0410$Model <- "Dorval et al. 2015"

ps.0410$sl <- ps.0410$standardLength_mm
ps.0410$fage <- as.factor(ps.0410$age)

unique(ps.0410$age) # check to make sure all good
table(ps.0410$age) # look at sample size

# Dorval et al. 2015 with adding age zeros from Dataset 1
# need to subset data and choose specific columns for later row bind
age0.data1 <- original %>%
  filter(age==0)

ps.0410.addzero <- rbind(age0.data1,ps.0410)
ps.0410.addzero$Model <- "Dorval et al. 2015, plus age-0"

## comparing growth models #####
# setting up von Bertalanffy parameters for all models
vb1 <- vbFuns()
svvb <- list(Linf=250,K=0.3,t0=-2)

##1 original data set ######
fitvb.org <- nlsLM(sl~vb1(age,Linf,K,t0),data=original,start=svvb,trace=TRUE)
vbvals.org=coef(fitvb.org)
summary(fitvb.org)
bootvb.org <- nlsBoot(fitvb.org)
cbind(coef(fitvb.org),confint(bootvb.org))
vbvals.org <- as.data.frame(t(cbind(coef(fitvb.org),confint(bootvb.org))))
vbvals.org$Dataset <- "Current Study"
rownames(vbvals.org)<-c("coeff","LCI","UCI")
vbvals.org$Parameter <- rownames(vbvals.org)

# predicting another way for plotting
predict2vb.org <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2vb.org <- Boot(fitvb.org,f=predict2vb.org)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.org <- data.frame(ages,
                           predict(fitvb.org,data.frame(age=ages)),
                           confint(f.boot2vb.org))
names(preds1vb.org) <- c("age","fit","LCI","UCI")
headtail(preds1vb.org)
preds1vb.org$Model <- "Current Study"

# only predicting for observed ages 
preds2vb.org <- filter(preds1vb.org,age>=0,age<=9)
headtail(preds2vb.org)

##2 Dorval et al. 2015 data set ######
fitvb.0410 <- nlsLM(sl~vb1(age,Linf,K,t0),data=ps.0410,start=svvb,trace=TRUE)
#vbvals.0410=coef(fitvb.0410)
summary(fitvb.0410)
bootvb.0410 <- nlsBoot(fitvb.0410)
cbind(coef(fitvb.0410),confint(bootvb.0410))
vbvals.0410<-as.data.frame(t(cbind(coef(fitvb.0410),confint(bootvb.0410))))
vbvals.0410$Dataset <- "Dorval et al. 2015"
rownames(vbvals.0410)<-c("coeff","LCI","UCI")
vbvals.0410$Parameter <- rownames(vbvals.0410)

# predicting another way for plotting
predict2vb.0410 <- function(x) predict(x,data.frame(age=ages))
ages <- 0:10
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(-1,11,by=0.2)
f.boot2vb.0410 <- Boot(fitvb.0410,f=predict2vb.0410)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.0410 <- data.frame(ages,
                            predict(fitvb.0410,data.frame(age=ages)),
                            confint(f.boot2vb.0410))
names(preds1vb.0410) <- c("age","fit","LCI","UCI")
headtail(preds1vb.0410)
preds1vb.0410$Model <- "Dorval et al. 2015"
# only predicting for observed ages 
preds2vb.0410 <- filter(preds1vb.0410,age>=0,age<=10)
headtail(preds2vb.0410)

##3: Dorval et al. 2015 with age 0 added from Dataset 1####
fitvb.0410.addzero <- nlsLM(sl~vb1(age,Linf,K,t0),data=ps.0410.addzero,start=svvb,trace=TRUE)
summary(fitvb.0410.addzero)
bootvb.0410.addzero <- nlsBoot(fitvb.0410.addzero)
cbind(coef(fitvb.0410.addzero),confint(bootvb.0410.addzero))
vbvals.0410.addzero<-as.data.frame(t(cbind(coef(fitvb.0410.addzero),confint(bootvb.0410.addzero))))
vbvals.0410.addzero$Dataset<-"Dorval et al. 2015, plus age-0"
rownames(vbvals.0410.addzero)<-c("coeff","LCI","UCI")
vbvals.0410.addzero$Parameter <- rownames(vbvals.0410.addzero)

# predicting another way for plotting
predict2vb.0410.addzero <- function(x) predict(x,data.frame(age=ages))
ages <- 0:10
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(-1,11,by=0.2)
f.boot2vb.0410.addzero <- Boot(fitvb.0410.addzero,f=predict2vb.0410.addzero)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.0410.addzero <- data.frame(ages,
                                    predict(fitvb.0410.addzero,data.frame(age=ages)),
                                    confint(f.boot2vb.0410.addzero))
names(preds1vb.0410.addzero) <- c("age","fit","LCI","UCI")
headtail(preds1vb.0410.addzero)
preds1vb.0410.addzero$Model <- "Dorval et al. 2015, plus age-0"
# only predicting for observed ages 
preds2vb.0410.addzero <- filter(preds1vb.0410.addzero,age>=0,age<=10)
headtail(preds2vb.0410.addzero)

# comparing output ######
# getting coefficients for each model

# only looking at original and dorval FIGURE 6 in MS
# plotting with ggplot
preds2.two <- rbind(preds2vb.org,preds2vb.0410)
preds2.two$Model <- factor(preds2.two$Model,levels=c("Current Study","Dorval et al. 2015"))

original.sub <- original %>% select(age,sl,Model)
ps.0410.sub <- ps.0410 %>% select(age,sl,Model)
data.two <- rbind(original.sub,ps.0410.sub)
data.two$Model <- factor(data.two$Model,levels=c("Current Study","Dorval et al. 2015"))

col.pal2 <- c("#332288","#AB2100")

DorvalCurrent <- ggplot() + 
  geom_point(data=data.two,aes(y=sl,x=age, color=Model),position=position_dodge(width=0.3), size=2,alpha=0.3)+
  geom_line(data=preds2.two,aes(y=fit,x=age, color=Model),linewidth=1)+ #line only for observed ages
  geom_ribbon(data=preds2.two,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal2)+
  scale_fill_manual(values=col.pal2)+
  labs(x="Age (years)", y = "Standard length (mm)",color="Dataset",fill="Dataset")+
  scale_x_continuous(limits=c(-0.5,10),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(20,300),breaks=seq(20,300,20))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

DorvalCurrent
ggsave("Output/GrowthModel_2dataset_comp.jpeg", height=5, width=6)

#Put original and Dorval et al. 2015 VBGM lines on #10 data (Dorval plus age-0) FIGURE 7 in MS
col.pal3 <- c("#332288","#AB2100","#44AA99") 
preds2.three <- rbind(preds2vb.org,preds2vb.0410,preds2vb.0410.addzero)
preds2.three$Model <- factor(preds2.three$Model,levels=c("Current Study","Dorval et al. 2015","Dorval et al. 2015, plus age-0"))

#just plot ps.0410.addzero
Dorvaladdzero <- ggplot() + 
  geom_line(data=preds2.three,aes(y=fit,x=age, color=Model),linewidth=1)+ #line only for observed ages
  geom_ribbon(data=preds2.three,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal3)+
  scale_fill_manual(values=col.pal3)+
  labs(x="Age (years)", y = "Standard length (mm)",color="Dataset",fill="Dataset")+
  scale_x_continuous(limits=c(-0.5,10),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(20,300),breaks=seq(20,300,20))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

Dorvaladdzero
ggsave("Output/GrowthModel_IllustrateSubZero.jpeg", height=5, width=6)
