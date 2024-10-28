# Load the required packages -----------------------------------------------

{ library(dplyr)
library(ggpubr)
library(moments)
library(car)
library(emmeans)
library(multcomp)
library(plotrix)
library(ggplot2)
library(ccplot)
library(pracma)
library(rcompanion)
library(pscl)
  library(performance)
  library(lme4)
  library(nlme)
  library(viridis)
}

# Load and check the data -------------------------------------------------

# Load the oyster count data
Oyst <- read.csv("./Oyster_count.csv")

# Check variable classifications
str(Oyst)
# Correct
Oyst$Frame_ID <- as.factor(Oyst$Frame_ID)
Oyst$Configuration <- as.factor(Oyst$Configuration)
Oyst$Panel <- as.factor(Oyst$Panel)
Oyst$Species <- as.factor(Oyst$Species)
# Check again
str(Oyst)

# Separate the species
Mgigas <- filter(Oyst, Species=="M_gigas")
Aephi <- filter(Oyst, Species=="A_ephi")

# Pacific oyster count ----------------------------------------------------

sum(Mgigas$Count)
# 28

# BY CONFIGURATION
Mgigas_sum <- Mgigas %>% group_by(Configuration) %>% 
  summarise(sum=sum(Count))
Mgigas_sum
# Configuration   sum
# Grouped          19
# Random            9

ggplot(data=Mgigas_sum, aes(x=Configuration, y=sum, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Pacific oysters")+
  xlab("Configuration")+
  ylim(c(0, 20))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# BY PANEL
Mgigas_sum <- Mgigas %>% group_by(Panel) %>% 
  summarise(sum=sum(Count))
Mgigas_sum
# Panel     sum
# Flat       18
# Ripple      3
# 2.5 cm      1
# 5 cm        6

Mgigas_sum$Panel = factor(Mgigas_sum$Panel, levels = c("Flat", "Ripple", "25_cm", "5_cm"), ordered=TRUE)

ggplot(data=Mgigas_sum, aes(x=Panel, y=sum, fill=Panel))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Pacific oysters")+
  xlab("Panel")+
  ylim(c(0, 20))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Analysis
mod1 <- aov(Count~Configuration+Panel, data=Mgigas)
summary(mod1)
# No effect of configuration
# F(1, 211) = 1.44, p = 0.23
# No effect of panel
# F(3, 211) = 2.10, p = 0.10

# Saddle oyster count ----------------------------------------------------

sum(Aephi$Count)
# 44

# BY CONFIGURATION
Aephi_sum <- Aephi %>% group_by(Configuration) %>% 
  summarise(sum=sum(Count))
Aephi_sum
# Configuration   sum
# Grouped          19
# Random           25

ggplot(data=Aephi_sum, aes(x=Configuration, y=sum, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Saddle oysters")+
  xlab("Configuration")+
  ylim(c(0, 25))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# BY PANEL
Aephi_sum <- Aephi %>% group_by(Panel) %>% 
  summarise(sum=sum(Count))
Aephi_sum
# Panel      sum
# Flat         6
# Ripple      12
# 2.5 cm      17
# 5 cm         9

Aephi_sum$Panel = factor(Aephi_sum$Panel, levels = c("Flat", "Ripple", "25_cm", "5_cm"), ordered=TRUE)

ggplot(data=Aephi_sum, aes(x=Panel, y=sum, fill=Panel))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Saddle oysters")+
  xlab("Panel")+
  ylim(c(0, 25))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Load the data -----------------------------------------------------------

bio <- read.csv("./Biodiversity_data.csv")
head(bio)

# Check and edit variables classes
str(bio)
bio$Frame_ID <- as.factor(bio$Frame_ID)
bio$Configuration <- as.factor(bio$Configuration)
bio$Panel <- as.factor(bio$Panel)
bio$Panel_ID <- as.factor(bio$Panel_ID)
bio$Type <- as.factor(bio$Type)
bio$Measure <- as.factor(bio$Measure)
str(bio)

# Separate INNS species
Ceum <- filter(bio, Species=="C_eumyota")
Wsub <- filter(bio, Species=="W_subatra")

# C eumyota cover ----------------------------------------------------

max(Ceum$Count_perc)
# 7

# BY CONFIGURATION
Ceum_avg <- Ceum %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg
# Configuration    mean      se
# Grouped         0.120  0.0538
# Random          0.278  0.0947

ggplot(data=Ceum_avg, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of Corella eumyota")+
  xlab("Configuration")+
  ylim(c(0, 0.4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# BY PANEL
Ceum_avg <- Ceum %>% group_by(Panel) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg 
# Panel       mean       se
# cm_25        0.5    0.176 
# cm_5       0.667    0.255 
# Flat      0.0139   0.0139
# Ripple         0        0        

Ceum_avg$Panel = factor(Ceum_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=Ceum_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of Corella eumyota")+
  xlab("Panel")+
  ylim(c(0, 1))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Analysis
mod1 <- aov(Count_perc~Configuration+Panel, data=Ceum)
summary(mod1)
# No effect of configuration
# F(1, 211) = 2.34, p = 0.13
# Significant effect of tile complexity
# F(3, 211) = 9.59, p < 0.0001

TukeyHSD(mod1, which="Panel")
# Flat-Ripple  0.9995198
# Flat-cm_25   0.0100477 *
# Flat-cm_5    0.0002037 *
# Ripple-cm_25 0.0075486 *
# Ripple-cm_5  0.0001409 *
# cm_5-cm_25   0.7860050

# W subatra cover ----------------------------------------------------

max(Wsub$Count_perc)
# 20

# BY CONFIGURATION
Wsub_avg <- Wsub %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg
# Configuration    mean      se
# Grouped          1.81   0.371
# Random          0.889   0.207

ggplot(data=Wsub_avg, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average cover of Watersipora subatra")+
  xlab("Configuration")+
  ylim(c(0, 3))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# BY PANEL
Wsub_avg <- Wsub %>% group_by(Panel) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg 
# Panel       mean       se
# cm_25       2.28    0.484 
# cm_5        5.78    0.817 
# Flat      0.0278   0.0195
# Ripple         0        0     

Wsub_avg$Panel = factor(Wsub_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=Wsub_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average cover of Watersipora subatra")+
  xlab("Panel")+
  ylim(c(0, 8))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

mod1 <- aov(Count_perc~Configuration+Panel, data=Wsub)
summary(mod1)
# Significant effect of configuration
# F(1, 211) = 8.95, p = 0.003
# Significant effect of panel
# F(3, 211) = 64.07, p < 0.0001

TukeyHSD(mod1, which="Panel")
# Flat-Ripple   0.9998584
# Flat-cm_25    0.0000144 *
# Flat-cm_5     0.0000000 *
# Ripple-cm_25  0.0000109 *
# Ripple-cm_5   0.0000000 *
# cm_25-cm_5    0.0000000 *
