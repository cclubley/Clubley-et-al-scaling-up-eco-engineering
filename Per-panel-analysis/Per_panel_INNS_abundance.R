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
  library(glmmTMB)
  library(scales)
}

# Load and check the data -------------------------------------------------

# Load the oyster count data
Oyst <- read.csv("./Oyster_count.csv")

# Check variable classifications
str(Oyst)
# Correct
Oyst$Frame_ID <- as.factor(Oyst$Frame_ID)
Oyst$Configuration <- as.factor(Oyst$Configuration)
Oyst$Panel_type <- as.factor(Oyst$Panel_type)
Oyst$Species <- as.factor(Oyst$Species)
# Check again
str(Oyst)

# Separate the species
Mgigas <- filter(Oyst, Species=="M_gigas")
Aephi <- filter(Oyst, Species=="A_ephi")

# Pacific oyster count ----------------------------------------------------

sum(Mgigas$Count)
# 28

# Summary table - configuration
Mgigas_sum <- Mgigas %>% group_by(Configuration) %>% 
  summarise(sum=sum(Count))
Mgigas_sum
# Configuration   sum
# Grouped          19
# Random            9

ggplot(data=Mgigas_sum, aes(x=Configuration, y=sum))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Pacific oysters")+
  xlab("Configuration")+
  ylim(c(0, 20))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Summary table - panel
Mgigas_sum <- Mgigas %>% group_by(Panel_type) %>% 
  summarise(sum=sum(Count))
Mgigas_sum
# Panel_type     sum
# Flat            18
# Ripple           3
# 2.5 cm           1
# 5 cm             6

Mgigas_sum$Panel_type = factor(Mgigas_sum$Panel_type, levels = c("Flat", "Ripple", "25_cm", "5_cm"), ordered=TRUE)
ggplot(data=Mgigas_sum, aes(x=Panel_type, y=sum))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Pacific oysters")+
  xlab("Panel")+
  ylim(c(0, 20))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Check for normality
ggdensity(Mgigas$Count)
ggqqplot(Mgigas$Count)

# Shapiro-Wilk test
shapiro.test(Mgigas$Count)

# Check for zero-inflation
prop_zeros <- mean(Mgigas$Count == 0)
prop_zeros # 0.92

# Model
set.seed(1234)
mod1 <- glmmTMB(Count ~ Configuration+Panel_type + (1|Configuration/Frame_ID), ziformula=~1, data=Mgigas, family=poisson)
check_model(mod1)

summary(mod1)
#                     Estimate  Std. Error  z value  Pr(>|z|)  
# (Intercept)          -1.8153      1.1756   -1.544    0.1226  
# ConfigurationRandom  -0.9128      0.5838   -1.564    0.1179  
# Panel5_cm             1.8956      1.2332    1.537    0.1243  
# PanelFlat             2.4063      1.1322    2.125    0.0336 *
# PanelRipple           0.4383      1.2570    0.349    0.7273   

# Saddle oyster count ----------------------------------------------------

sum(Aephi$Count)
# 44

# Summary table - configuration
Aephi_sum <- Aephi %>% group_by(Configuration) %>% 
  summarise(sum=sum(Count))
Aephi_sum
# Configuration   sum
# Grouped          19
# Random           25

ggplot(data=Aephi_sum, aes(x=Configuration, y=sum))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Total recruited Saddle oysters")+
  xlab("Configuration")+
  ylim(c(0, 25))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Summary table - by panel
Aephi_sum <- Aephi %>% group_by(Panel_type) %>% 
  summarise(sum=sum(Count))
Aephi_sum
# Panel_type      sum
# Flat              6
# Ripple           12
# 2.5 cm           17
# 5 cm              9

Aephi_sum$Panel_type = factor(Aephi_sum$Panel_type, levels = c("Flat", "Ripple", "25_cm", "5_cm"), ordered=TRUE)
ggplot(data=Aephi_sum, aes(x=Panel_type, y=sum))+
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
bio$Panel_type <- as.factor(bio$Panel_type)
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

# Summary table - configuration
Ceum_avg <- Ceum %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg
# Configuration    mean      se
# Grouped         0.120  0.0538
# Random          0.278  0.0947

ggplot(data=Ceum_avg, aes(x=Configuration, y=mean))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of Corella eumyota")+
  xlab("Configuration")+
  ylim(c(0, 0.4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Summary table - panel
Ceum_avg <- Ceum %>% group_by(Panel_type) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg 
# Panel       mean       se
# cm_25        0.5    0.176 
# cm_5       0.667    0.255 
# Flat      0.0139   0.0139
# Ripple         0        0        

Ceum_avg$Panel_type = factor(Ceum_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=Ceum_avg, aes(x=Panel_type, y=mean))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of Corella eumyota")+
  xlab("Panel")+
  ylim(c(0, 1))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Check for normality
ggdensity(Ceum$Count_perc)
ggqqplot(Ceum$Count_perc)

# Shapiro-Wilk test
shapiro.test(Ceum$Count_perc)

Ceum$n2 <- rescale(Ceum$Count_perc, to=c(0.01, 0.99))
set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration+Panel_type + (1|Configuration/Frame_ID), data=Ceum, family=beta_family())
check_model(mod1)

summary(mod1)
#                     Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)          -2.6662      0.1818  -14.668    <2e-16 ***
# ConfigurationRandom   0.1006      0.1302    0.773    0.4396    
# Panelcm_5             0.1078      0.2194    0.491    0.6234    
# PanelFlat            -0.3440      0.1937   -1.776    0.0757 .  
# PanelRipple          -0.3585      0.1938   -1.850    0.0643 .    
 
# W subatra cover ----------------------------------------------------

max(Wsub$Count_perc)
# 20

# Summary table - configuration
Wsub_avg <- Wsub %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg
# Configuration    mean      se
# Grouped          1.81   0.371
# Random          0.889   0.207

ggplot(data=Wsub_avg, aes(x=Configuration, y=mean))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average cover of Watersipora subatra")+
  xlab("Configuration")+
  ylim(c(0, 3))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Summary table - panel
Wsub_avg <- Wsub %>% group_by(Panel_type) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg 
# Panel       mean       se
# cm_25       2.28    0.484 
# cm_5        5.78    0.817 
# Flat      0.0278   0.0195
# Ripple         0        0     

Wsub_avg$Panel_type = factor(Wsub_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=Wsub_avg, aes(x=Panel_type, y=mean))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average cover of Watersipora subatra")+
  xlab("Panel")+
  ylim(c(0, 8))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Check for normality
ggdensity(Wsub$Count_perc)
ggqqplot(Wsub$Count_perc)

# Shapiro-Wilk test
shapiro.test(Wsub$Count_perc)

Wsub$n2 <- rescale(Wsub$Count_perc, to=c(0.01, 0.99))
set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration+Panel_type + (1|Configuration/Frame_ID), data=Wsub, family=beta_family())
check_model(mod1)

summary(mod1)
#                     Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)          -2.0574      0.2002  -10.276   < 2e-16 ***
# ConfigurationRandom  -0.4165      0.1887   -2.207    0.0273 *  
# Panelcm_5             1.1860      0.2278    5.207  1.92e-07 ***
# PanelFlat            -0.8259      0.1983   -4.164  3.13e-05 ***
# PanelRipple          -0.8573      0.1999   -4.289  1.80e-05 ***    
 