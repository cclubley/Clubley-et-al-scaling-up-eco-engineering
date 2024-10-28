# Load the required packages ----------------------------------------------

{ library(ggpubr)
library(moments)
library(pracma)
library(rcompanion)
library(car)
  library(FSA)
  library(dplyr)
library(plotrix)
  library(performance)
  library(lme4)
  library(ggplot2)
  library(ccplot)
  library(multcomp)
  library(emmeans)
  library(viridis)
}

# Load data ---------------------------------------------------------------

# Sediment cover
bio <- read.csv("./Biodiversity_data.csv")
head(bio)

# Check and edit variables classes
str(bio)
bio$Frame_ID <- as.factor(bio$Frame_ID)
bio$Configuration <- as.factor(bio$Configuration)
bio$Panel_ID <- as.factor(bio$Panel_ID)
bio$Panel <- as.factor(bio$Panel)
bio$Type <- as.factor(bio$Type)
bio$Measure <- as.factor(bio$Measure)

# Filter just the sediment data
Sedc <- filter(bio, Species=="Sediment")

# Sediment depth
Sedd <- read.csv("./Sediment_depth_data.csv")

# Check variable classifications
str(Sedd)
# Correct
Sedd$Frame_ID <- as.factor(Sedd$Frame_ID)
Sedd$Configuration <- as.factor(Sedd$Configuration)
Sedd$Panel <- as.factor(Sedd$Panel)
# Check again
str(Sedd)

# Calculate summary tables ------------------------------------------------

# Overall average cover
mean(Sedc$Count_perc)
std.error(Sedc$Count_perc)
# 62.06019 +- 2.191206

# The effect of configuration on cover
Sedc %>% group_by(Configuration) %>% 
  summarise(mean = mean(Count_perc), se = std.error(Count_perc))
# Configuration    Mean     SE
# Grouped          63.3   3.05
# Random           60.8   3.16

# The effect of panel on cover
Sedc %>% group_by(Panel) %>% 
  summarise(mean = mean(Count_perc), se = std.error(Count_perc))
# Panel       Mean     SE
# Flat        65.5   4.07
# Rippl       60.8   3.78
# 25_cm       60.7   5.06
# 5_cm        56.4   5.01

# Overall average depth
mean(Sedd$Sed_depth)
std.error(Sedd$Sed_depth)
# 0.195 +- 0.02

# The effect of configuration on depth
Sedd %>% 
  group_by(Configuration) %>% 
  summarise(mean = mean(Sed_depth),
            std = std.error(Sed_depth))
# Configuration   Mean     SE
# Grouped        0.142  0.016
# Random         0.231  0.036

# The effect of panel on depth
Sedd %>% 
  group_by(Panel) %>% 
  summarise(mean = mean(Sed_depth),
            std = std.error(Sed_depth))
# Panel      Mean    SE
# Flat      0.133  0.02
# Ripple    0.280  0.06
# 25_cm     0.127  0.02
# 5_cm      0.240  0.05

# Plot cover --------------------------------------------------------------

# Plot configuration on its own
sedc_avg <- Sedc %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))

ggplot(data=sedc_avg, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sediment cover (%)")+
  xlab("Configuration")+
  ylim(c(0, 100))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Plot Panel on its own
sedc_avg <- Sedc %>% group_by(Panel) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))

# Re-order the tile types
sedc_avg$Panel = factor(sedc_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=sedc_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sediment cover (%)")+
  xlab("Panel")+
  ylim(c(0, 100))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Plot depth --------------------------------------------------------------

# Plot arrangement on its own
sedd_avg <- Sedd %>% group_by(Configuration) %>% 
  summarise(mean=mean(Sed_depth), se=std.error(Sed_depth))

ggplot(data=sedd_avg, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sediment depth (cm)")+
  xlab("Configuration")+
  ylim(c(0, 0.4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Plot panel on its own
sedd_avg <- Sedd %>% group_by(Panel) %>% 
  summarise(mean=mean(Sed_depth), se=std.error(Sed_depth))

# Re-order the tile types
sedd_avg$Panel = factor(sedd_avg$Panel, levels = c("Flat", "Ripple", "25_cm", "5_cm"), ordered=TRUE)

ggplot(data=sedd_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sediment depth (cm)")+
  xlab("Panel")+
  ylim(c(0, 0.4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")
