rm(list = ls())

#load libraries
library(openxlsx)
library(tidyverse)
library(viridis)
library(scales)
library(ggpubr)
library(FSA)

####Data wrangling####
#Dataset of Pacific Sardine ages, one reader per fish
sard.u.raw <- read.xlsx("Data/sard_u_raw.xlsx")

#Length Frequency - All years combined
sl.freq.raw<-ggplot(sard.u.raw, aes(x=standardLength_mm))+
  geom_histogram(binwidth = 20, boundary = 0, fill="black")+
  labs(x="Standard Length (mm)",y="Frequency")+
  scale_x_continuous(breaks=seq(20,300,20))+
  theme_linedraw()+
  theme(panel.grid=element_blank(),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))
sl.freq.raw

age.freq.raw<-ggplot(sard.u.raw, aes(x=as.factor(age)))+
  geom_histogram(stat="count", fill="black")+
  labs(x="Age (years)",y="Frequency")+
  theme_linedraw()+
  theme(panel.grid=element_blank(),
        axis.title.y = element_text(size=12,margin=margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(size=12,margin=margin(t=10,r=0,b=0,l=0)),
        plot.title = element_text(hjust=0.5))
age.freq.raw

ggarrange(sl.freq.raw, age.freq.raw, ncol=2, nrow=1, labels = c("A","B"),
          font.label= list(family="serif",face="bold"))

ggsave("Output/Sardine_SL_Age_FreqRaw.jpeg", width=11, height=5.5)

#Comparison of mean length-at-age in Raw data####
#Normality test of length data
shapiro.test(sard.u.raw$standardLength_mm)

#homogeneity of variance on length data by age
bartlett.test(standardLength_mm~age, data = sard.u.raw)

#Pairwise Comparison of mena lengths at age for each age class
kruskal.test(standardLength_mm~age, data = sard.u.raw)
dunnTest(standardLength_mm~as.factor(age), data = sard.u.raw, method = "hochberg")

#Mean length at age and Outliers
sard.u.raw$sl <- sard.u.raw$standardLength_mm
sard.u.raw$fage <- as.factor(sard.u.raw$age)

ssard <- sard.u.raw %>% 
  group_by(fage) %>% 
  summarize(
    count = n(),
    mean = mean(sl, na.rm = TRUE),
    sd2 = 2*sd(sl, na.rm = TRUE)
  )
old <- options(pillar.sigfig = 7)
ssard

#mean length
msard <- sard.u.raw %>% 
  summarize(
    count = n(),
    mean = mean(sl, na.rm = TRUE),
    sd = sd(sl, na.rm = TRUE)
  )
msard

#mean age
mage <- sard.u.raw %>%
  summarize(
    count = n(),
    mean = mean(age, na.rm = TRUE),
    sd = sd(age, na.rm = TRUE)
  )
mage

####Map####
# once you have the RDS file saved, can just import that instead of connecting to database
#requires the following to be installed: rgdal, raster, mapproj

states <- map_data("state")
can <- map_data("world", region = "Canada")
mex <- map_data("world", region = "Mexico")
north.am.mex <- rbind(can, states, mex)

#Dataset with count of each collection (number of sardine aged in a haul)
small <- read.xlsx("Data/small.xlsx")
table(small$cruise)

#map black and white
ggplot()+
  geom_polygon(data = north.am.mex, aes(x = long, y = lat, group = group), 
               fill = "grey85", color = "black")+
  # add points
  geom_point(data = small, aes(x = startLongDecimal, y = startLatDecimal, size = n), shape = 1)+
 
  scale_x_continuous(breaks=c(-126, -122, -118),  
                     labels=c("\u2212126°W","\u2212122°W","\u2212118°W"), expand = c(0, 0))+ 
  scale_y_continuous(breaks=seq(35,50,5), labels=paste0(seq(35,50,5),'°N'),expand = c(0, 0))+
  coord_map(xlim=c(-127.5, -117.0), ylim=c(30.5, 50.5))+
  labs(x = "Longitude", y = "Latitude")+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill="white"),
        axis.title.x = element_text(margin=margin(t=5,r=0,b=0,l=0)),
        axis.title.y = element_text(margin=margin(t=0,r=5,b=0,l=0)),
        plot.title = element_text(hjust=0.5))
ggsave("Output/Map_BubblePlotBW.jpeg",width=6, height=8)

####Length Age Frequency Dorval####
sard.0410.u <- read.xlsx("Data/sard0410_u.xlsx")
#Mean Length at Age
sard.0410.u$fage <- as.factor(sard.0410.u$age)
sard.0410.u$sl <- sard.0410.u$standardLength_mm
ssard0410 <- sard.0410.u %>% 
  group_by(fage) %>% 
  summarize(
    count = n(),
    mean = mean(sl, na.rm = TRUE),
    sd = sd(sl, na.rm = TRUE)
  )
ssard0410
