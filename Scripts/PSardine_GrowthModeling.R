rm(list = ls())
#load packages
library(FSA)
library(nlme)
library(mgcv)
library(car)
library(AICcmodavg)
library(minpack.lm)
library(nlstools)
library(ggplot2)
library(patchwork)
library(plyr)
library(openxlsx)
library(dplyr)
library(tidyverse)

#import data
#sard.u <- read.xlsx("Data/sard.u.xlsx")
sard.u <- read.xlsx("Data/sard_u_raw.xlsx")

unique(sard.u$age) # check to make sure all good
table(sard.u$age) # look at sample size

sard.u$fage <- as.factor(sard.u$age)

######### comparing growth models ##########
#name data
data1=sard.u
sl=data1$standardLength_mm
age=data1$age

#name functions from FSA
l1 <- logisticFuns()
g1 <- GompertzFuns()
vb1 <- vbFuns()

#gompertz
svG1 <- list(Linf=250,gi=0.5,ti=0.2)
fitG1 <- nlsLM(sl~g1(age,Linf,gi,ti),data=data1,start=svG1,trace=TRUE)
gvals=coef(fitG1)
summary(fitG1)
bootG1 <- nlsBoot(fitG1)
cbind(coef(fitG1),confint(bootG1))

# predicting another way for plotting
predict2G1 <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2G1 <- Boot(fitG1,f=predict2G1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1G1 <- data.frame(ages,
                       predict(fitG1,data.frame(age=ages)),
                       confint(f.boot2G1))
names(preds1G1) <- c("age","fit","LCI","UCI")
headtail(preds1G1)
preds1G1$Model <- "Gompertz"
# only predicting for observed ages 
preds2G1 <- filter(preds1G1,age>=0,age<=9)
headtail(preds2G1)


#logistic
svL1 <- list(Linf=250,gninf=0.5,ti=2)
fitL1 <- nlsLM(sl~l1(age,Linf,gninf,ti),data=data1,start=svL1,trace=TRUE)
lvals=coef(fitL1)
summary(fitL1)
bootL1 <- nlsBoot(fitL1)
cbind(coef(fitL1),confint(bootL1))

# predicting another way for plotting
predict2L1 <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2L1 <- Boot(fitL1,f=predict2L1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1L1 <- data.frame(ages,
                       predict(fitL1,data.frame(age=ages)),
                       confint(f.boot2L1))
names(preds1L1) <- c("age","fit","LCI","UCI")
headtail(preds1L1)
preds1L1$Model <- "Logistic"
# only predicting for observed ages 
preds2L1 <- filter(preds1L1,age>=0,age<=9)
headtail(preds2L1)

#von bertalanffy
svvb <- list(Linf=250,K=0.3,t0=-2)
fitvb <- nlsLM(sl~vb1(age,Linf,K,t0),data=data1,start=svvb,trace=TRUE)
vbvals=coef(fitvb)
summary(fitvb)
bootvb <- nlsBoot(fitvb)
cbind(coef(fitvb),confint(bootvb))

# predicting another way for plotting
predict2vb <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2vb <- Boot(fitvb,f=predict2vb)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb <- data.frame(ages,
                       predict(fitvb,data.frame(age=ages)),
                       confint(f.boot2vb))
names(preds1vb) <- c("age","fit","LCI","UCI")
headtail(preds1vb)
preds1vb$Model <- "von Bertalanffy"
# only predicting for observed ages 
preds2vb <- filter(preds1vb,age>=0,age<=9)
headtail(preds2vb)

#AIC values
aictab(list(fitvb,fitL1,fitG1),c("von Bertalanffy","logistic","Gompertz"),second.ord = FALSE)

#obtain BIC values
BIC(fitL1,fitG1,fitvb)

#view parameter ests
vbvals
gvals
lvals

# plotting with ggplot
preds2.all <- rbind(preds2vb,preds2G1,preds2L1)
preds2.all$Model <- factor(preds2.all$Model,levels=c("von Bertalanffy","Gompertz","Logistic"))

col.pal <- c("#332288","#117733","#DDCC77")

ggplot() + 
  geom_point(data=sard.u,aes(y=sl,x=age),size=2,alpha=0.1)+
  #geom_line(data=preds1.all,aes(y=fit,x=age,color=Model),linewidth=1,linetype=2)+ #line over entire range
  geom_line(data=preds2.all,aes(y=fit,x=age, color=Model),linewidth=1)+ #line only for observed ages
  geom_ribbon(data=preds2.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal)+
  scale_fill_manual(values=col.pal)+
  labs(x="Age (years)", y = "Standard length (mm)")+
  scale_x_continuous(limits=c(-1,10),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(20,300),breaks=seq(20,300,20))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

ggsave("Output/GrowthModel_comp.jpeg", height=5, width=6)