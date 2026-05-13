rm(list = ls())
#Simulations to examine potential bias in our sample distribution
#1) Simulating 500 lengths for each age class from a normal distribution centered on the mean length-at-age 
#with a standard deviation equal to the standard deviation of length at that age; 2) Randomly subsample 
#abundant ages to an n (need to spend some time figuring out what a good n is; maybe 100, maybe 250, one 
#that accurately captures the length distribution at age) and simulate length data for age classes with 
#fewer than the n (as in 1). Then these datasets need to be run through Brittany's von Bertalanffy 
#growth model code, which you do have access to.

library(tidyverse)
library(openxlsx)
library(FSA)
library(nlme)
library(mgcv)
library(car)
library(minpack.lm)
library(nlstools)
library(patchwork)

#Bring in data
sardinedata <- read.xlsx("Data/sard_u_raw.xlsx", sheet = 1) #assessment data

mla <- sardinedata %>%
  group_by(age) %>%
  summarise(avg = mean(standardLength_mm),
            std = sd(standardLength_mm),
            n = n()
  )
mla

#function to generate a normal distribution from mean and sd
normv <- function(n, mean, sd){
  out <- rnorm(n*length(mean), mean = mean, sd = sd)
  return(matrix(out, ncol = n, byrow = FALSE))
}

####1 - Simulate 500 lengths for each age class####

#using function above, generate random sample of 500 per age class, turn into dataframe
set.seed(1)
sim1.df <- as.data.frame(normv(500, mean=mla$avg, sd=mla$std))
#rename rows as ages
rownames(sim1.df) <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

#pivot data.frame from wide to long, Vid is not necessary for data, but is for
#gather function (a key), so all the V1s are stacked, then V2s...
sim1 <- sim1.df %>% gather(Vid, standardLength_mm)
#Add back in identifier of age, first V1 is age 0, second is age 1 and so on
#to age 9 for this dataset.
sim1$age <- c(0:9)

#To see how far the simulated dataset falls from the true mean length at age
#and standard deviation
rehash <- sim1 %>%
  group_by(age) %>%
  summarise(avg = mean(standardLength_mm),
            std = sd(standardLength_mm))
old <- options(pillar.sigfig = 7)
rehash #compare with #mla

#plot simulated dataset: age frequency
ggplot(sim1, aes(age))+
  geom_histogram()+
  theme_classic(base_size = 14)

#plot simulated dataset: length frequency
ggplot(sim1, aes(standardLength_mm))+
  geom_histogram(binwidth=20, boundary=0)+
  theme_classic(base_size = 14)

#get rid of Vid column
sim1 <- select(sim1, -Vid)
str(sim1)

#Run Sim1 through VBGF####
sim1$sl <- sim1$standardLength_mm
sim1$fage <- as.factor(sim1$age)

unique(sim1$age) # check to make sure all good
table(sim1$age)

#name functions from FSA
vb1 <- vbFuns()

data2=sim1
sl=data2$sl
age=data2$age

#von bertalanffy
svvb <- list(Linf=250,K=0.3,t0=-2)
fitvb.sim1 <- nlsLM(sl~vb1(age,Linf,K,t0),data=data2,start=svvb,trace=TRUE)
vbvals.sim1=coef(fitvb.sim1)
summary(fitvb.sim1)
bootvb.sim1 <- nlsBoot(fitvb.sim1)
cbind(coef(fitvb.sim1),confint(bootvb.sim1))

# predicting another way for plotting
predict2vb.sim1 <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2vb.sim1 <- Boot(fitvb.sim1,f=predict2vb.sim1)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.sim1 <- data.frame(ages,
                            predict(fitvb.sim1,data.frame(age=ages)),
                            confint(f.boot2vb.sim1))
names(preds1vb.sim1) <- c("age","fit","LCI","UCI")
headtail(preds1vb.sim1)
preds1vb.sim1$Model <- "Simulation #1"
# only predicting for observed ages 
preds2vb.sim1 <- filter(preds1vb.sim1,age>=0,age<=9)
headtail(preds2vb.sim1)

vbvals.sim1
#this will write over something that's already in output, check that it'll be the same
write.xlsx(preds1vb.sim1, "Data/preds1vb.sim1.xlsx")
write.xlsx(preds2vb.sim1, "Data/preds2vb.sim1.xlsx")

#

####2 - Subsample to n 250, then simulate up to the number####
#see ns for different age classes
table(sardinedata$age)

#subset ages 1-4 to 250 
ss14 <- subset(sardinedata, age <= 4 & age >= 1)
table(ss14$age)
ss14.250 <- ss14 %>%
  group_by(age) %>%
  sample_n(250, replace=FALSE)
table(ss14.250$age)

ss0 <- subset(sardinedata, age == 0)
table(ss0$age)
ss59 <- subset(sardinedata, age >= 5)
table(ss59$age)

#bolster 0 & 5-9 to 250 using mla f
#generate data for each age class up to 250 for each age class
set.seed(3)
sim0w <- as.data.frame(normv((250-mla$n[1]), mean=mla$avg[1], sd=mla$std[1]))
sim0l <- sim0w %>% gather(Vid, standardLength_mm)
sim0l$age <- 0
sim0 <- select(sim0l, standardLength_mm, age)

sim5w <- as.data.frame(normv((250-mla$n[6]), mean=mla$avg[6], sd=mla$std[6]))
sim5l <- sim5w %>% gather(Vid, standardLength_mm)
sim5l$age <- 5
sim5 <- select(sim5l, standardLength_mm, age)

sim6w <- as.data.frame(normv((250-mla$n[7]), mean=mla$avg[7], sd=mla$std[7]))
sim6l <- sim6w %>% gather(Vid, standardLength_mm)
sim6l$age <- 6
sim6 <- select(sim6l, standardLength_mm, age)

sim7w <- as.data.frame(normv((250-mla$n[8]), mean=mla$avg[8], sd=mla$std[8]))
sim7l <- sim7w %>% gather(Vid, standardLength_mm)
sim7l$age <- 7
sim7 <- select(sim7l, standardLength_mm, age)

sim8w <- as.data.frame(normv((250-mla$n[9]), mean=mla$avg[9], sd=mla$std[9]))
sim8l <- sim8w %>% gather(Vid, standardLength_mm)
sim8l$age <- 8
sim8 <- select(sim8l, standardLength_mm, age)

sim9w <- as.data.frame(normv((250-mla$n[10]), mean=mla$avg[10], sd=mla$std[10]))
sim9l <- sim9w %>% gather(Vid, standardLength_mm)
sim9l$age <- 9
sim9 <- select(sim9l, standardLength_mm, age)

#combine raw data (ss59) with simulated data (sim5-9) to get 250 for each age class
ss0.s <- select(ss0, standardLength_mm, age)
ss59.s <- select(ss59, standardLength_mm, age)
ss059.250 <- rbind(ss0.s,ss59.s,sim0,sim5,sim6,sim7,sim8,sim9)
table(ss059.250$age)

#combine with 0-4 data
ss14.250.s <- select(ss14.250, standardLength_mm, age)
ss.250 <- rbind(ss14.250.s,ss059.250)
table(ss.250$age)

#Assess if n=250 is good
mla
ss.250.mla <- ss.250 %>% 
  group_by(age) %>%
  summarise(avg=mean(standardLength_mm),
            std=sd(standardLength_mm),
            n=n())
ss.250.mla
#250 is plenty

####2 - Subsample to n 100, then simulate up to the number####
table(sardinedata$age)

#subset ages 0-6 to 100
ss06 <- subset(sardinedata, age <= 6)
table(ss06$age)
ss06.100 <- ss06 %>%
  group_by(age) %>%
  sample_n(100, replace=FALSE)
table(ss06.100$age)

ss79 <- subset(sardinedata, age >= 7)
table(ss79$age)

#for each age class from mla
sim7w <- as.data.frame(normv((100-mla$n[8]), mean=mla$avg[8], sd=mla$std[8]))
sim7l <- sim7w %>% gather(Vid, standardLength_mm)
sim7l$age <- 7
sim7 <- select(sim7l, standardLength_mm, age)

sim8w <- as.data.frame(normv((100-mla$n[9]), mean=mla$avg[9], sd=mla$std[9]))
sim8l <- sim8w %>% gather(Vid, standardLength_mm)
sim8l$age <- 8
sim8 <- select(sim8l, standardLength_mm, age)

sim9w <- as.data.frame(normv((100-mla$n[10]), mean=mla$avg[10], sd=mla$std[10]))
sim9l <- sim9w %>% gather(Vid, standardLength_mm)
sim9l$age <- 9
sim9 <- select(sim9l, standardLength_mm, age)

#combine ss79 and sim data
ss79.s <- select(ss79, standardLength_mm, age)
ss79.100 <- rbind(ss79.s, sim7, sim8, sim9)
table(ss79.100$age)

#combine all age classes
ss06.100.s <- select(ss06.100, standardLength_mm, age)
ss.100 <- rbind(ss79.100, ss06.100.s)
table(ss.100$age)

#Assess if n=100 is good
mla
ss.100.mla <- ss.100 %>% 
  group_by(age) %>%
  summarise(avg=mean(standardLength_mm),
            std=sd(standardLength_mm),
            n=n())
ss.100.mla
#100 is actually plenty


#
####Run Sim2 through VBGF####
ss.100$sl <- ss.100$standardLength_mm
ss.100$fage <- as.factor(ss.100$age)

unique(ss.100$age) # check to make sure all good
table(ss.100$age)

#name functions from FSA
vb1 <- vbFuns()

data3=ss.100
sl=data3$sl
age=data3$age

#von bertalanffy
svvb <- list(Linf=250,K=0.3,t0=-2)
fitvb.sim2 <- nlsLM(sl~vb1(age,Linf,K,t0),data=data3,start=svvb,trace=TRUE)
vbvals.sim2=coef(fitvb.sim2)
summary(fitvb.sim2)
bootvb.sim2 <- nlsBoot(fitvb.sim2)
cbind(coef(fitvb.sim2),confint(bootvb.sim2))

# predicting another way for plotting
predict2vb.sim2 <- function(x) predict(x,data.frame(age=ages))
ages <- 0:9
# constructing mean lengths at ages with bootstrapped confidence intervals 
ages <- seq(0,9,by=0.2)
f.boot2vb.sim2 <- Boot(fitvb.sim2,f=predict2vb.sim2)
# placing ages, predictions and bootstrapped confidence in data frame for later use
preds1vb.sim2 <- data.frame(ages,
                            predict(fitvb.sim2,data.frame(age=ages)),
                            confint(f.boot2vb.sim2))
names(preds1vb.sim2) <- c("age","fit","LCI","UCI")
headtail(preds1vb.sim2)
preds1vb.sim2$Model <- "Simulation #2"
# only predicting for observed ages 
preds2vb.sim2 <- filter(preds1vb.sim2,age>=0,age<=9)
headtail(preds2vb.sim2)

vbvals.sim2
#this will write over something that's already in output, check that it'll be the same
write.xlsx(preds1vb.sim2, "Data/preds1vb.sim2.xlsx")
write.xlsx(preds2vb.sim2, "Data/preds2vb.sim2.xlsx")

#
####Run original data through VBGF for comparisons####
sard.s <- sardinedata %>% select(age, standardLength_mm)
sard.s$sl <- sard.s$standardLength_mm
sard.s$fage <- as.factor(sard.s$age)

unique(sard.s$age) # check to make sure all good
table(sard.s$age)

#von bertalanffy
svvb <- list(Linf=250,K=0.3,t0=-2)
data4=sard.s
sl=data4$sl
age=data4$age

fitvb.org <- nlsLM(sl~vb1(age,Linf,K,t0),data=data4,start=svvb,trace=TRUE)
vbvals.org=coef(fitvb.org)
summary(fitvb.org)
bootvb.org <- nlsBoot(fitvb.org)
cbind(coef(fitvb.org),confint(bootvb.org))

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
preds1vb.org$Model <- "Original"
# only predicting for observed ages 
preds2vb.org <- filter(preds1vb.org,age>=0,age<=9)
headtail(preds2vb.org)

vbvals.org
#this will write over something that's already in output, check that it'll be the same
write.xlsx(preds1vb.org, "Data/preds1vb.org.xlsx")
write.xlsx(preds2vb.org, "Data/preds2vb.org.xlsx")
#

####Graphing True, Sim #1, and Sim #2 ####
#combo true and sims for graphing
#reduce columns for true data (sardinedata)
sard.s$Model <-  "Original"
sard.s.1 <- sard.s %>% select(-"standardLength_mm")
sard.s.1 <- sard.s.1 %>% rename("age" = "age")

#add model name for Sim1 & Sim2
colnames(sim1) <- c("standardLength_mm", "age", "sl", "fage")
sim1$Model <- "Simulation #1"
sim1.1 <- sim1 %>% select(-"standardLength_mm")

ss.100 <- ss.100 %>% rename ("age" = "age")
ss.100$Model <- "Simulation #2"
ss.100.1 <- ss.100 %>% select(-"standardLength_mm")

combo.raw <- rbind(sard.s.1, sim1.1, ss.100.1)

#for graphing length at age comparison (2 sims, 1 true)
col.pal <- c("#332288", "#17AFBE", "#EE6A50")

LAboxcomp <- ggplot()+
  geom_boxplot(data=combo.raw, aes(fage, sl, color = Model), size=0.7,alpha=0.3)+
  #geom_point(data=combo.raw, aes(fage, sl, color = Model), size=1,alpha=0.3, position=position_jitterdodge(),alpha=0.3)+
  scale_color_manual(values=col.pal)+
  labs(x="Age (years)", y = "Standard length (mm)")+
  #scale_x_discrete(limits=c(0,10),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(20,300),breaks=seq(20,300,20))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

LAboxcomp
ggsave("Output/SimLAbox_comp.jpeg", height=5, width=6)

#VBGF graphing
#pull in saved datasets for model graphing, using preds2 for model
#preds1vb.org <- read.xlsx("Data/preds1vb.org.xlsx", sheet = 1)
preds2vb.org <- read.xlsx("Data/preds2vb.org.xlsx", sheet = 1)
#preds1vb.sim1 <- read.xlsx("Data/preds1vb.sim1.xlsx", sheet = 1)
preds2vb.sim1 <- read.xlsx("Data/preds2vb.sim1.xlsx", sheet = 1)
#preds1vb.sim2 <- read.xlsx("Data/preds1vb.sim2.xlsx", sheet = 1)
preds2vb.sim2 <- read.xlsx("Data/preds2vb.sim2.xlsx", sheet = 1)

#preds1.all <- rbind(preds1vb.org,preds1vb.sim1,preds1vb.sim2)
#preds1.all$Model <- factor(preds1.all$Model,levels=c("Original","Simulation #1","Simulation #2"))
preds2.all <- rbind(preds2vb.org,preds2vb.sim1,preds2vb.sim2)
preds2.all$Model <- factor(preds2.all$Model,levels=c("Original","Simulation #1","Simulation #2"))

#merged dataset = combo.raw
combo.raw$Model <- factor(combo.raw$Model, levels=c("Original","Simulation #1","Simulation #2"))

#VBGF graph
GrowthModelVBSim_comp <-
  ggplot() + 
  geom_point(data=combo.raw,aes(y=sl,x=age, color=Model),position=position_dodge(width=0.4), size=2,alpha=0.3)+
  #geom_line(data=preds1.all,aes(y=fit,x=age,color=Model),linewidth=1,linetype=2)+ #line over entire range
  geom_line(data=preds2.all,aes(y=fit,x=age, color=Model),linewidth=3)+ #line only for observed ages
  geom_ribbon(data=preds2.all,aes(x=age,ymin=LCI,ymax=UCI, fill=Model), alpha=0.5)+
  scale_color_manual(values=col.pal)+
  scale_fill_manual(values=col.pal)+
  labs(x="Age (years)", y = "Standard length (mm)")+
  scale_x_continuous(limits=c(-0.3,10),breaks=seq(0,10,1))+
  scale_y_continuous(limits=c(20,300),breaks=seq(20,300,20))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.8,0.2),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))

GrowthModelVBSim_comp
ggsave("Output/GrowthModelVBSim_comp.jpeg", height=5, width=6)
