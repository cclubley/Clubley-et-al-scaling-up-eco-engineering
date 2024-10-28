# Load the required packages ----------------------------------------------

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
Oyst_c <- read.csv("./Oyster_count.csv")

# Check variable classifications
str(Oyst_c)
# Correct
Oyst_c$Frame_ID <- as.factor(Oyst_c$Frame_ID)
Oyst_c$Configuration <- as.factor(Oyst_c$Configuration)
Oyst_c$Panel <- as.factor(Oyst_c$Panel)
Oyst_c$Species <- as.factor(Oyst_c$Species)
# Check again
str(Oyst_c)

# Calculate the biodiversity at the frame level
Oyst_c <- Oyst_c %>% group_by(Frame_ID, Configuration, Species) %>% 
  summarise(Count=sum(Count))

# Separate the species
C_Mgigas <- filter(Oyst_c, Species=="M_gigas")
C_Aephi <- filter(Oyst_c, Species=="A_ephi")

# Pacific oyster count ----------------------------------------------------

mod1 <- aov(Count ~ Configuration, data=C_Mgigas)
summary(mod1) 
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1    4.17    4.167    1.662   0.211
# Residuals          22   55.17    2.508               

# Saddle oyster count ----------------------------------------------------

mod1 <- aov(Count ~ Configuration, data=C_Aephi)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1     1.5      1.5    0.126   0.726
# Residuals          22   261.8     11.9               

# Load the data -----------------------------------------------------------

bio <- read.csv("./Biodiversity_data.csv")
head(bio)

# Check and edit variables classes
str(bio)
bio$Frame_ID <- as.factor(bio$Frame_ID)
bio$Configuration <- as.factor(bio$Configuration)
bio$Panel <- as.factor(bio$Panel)
bio$Type <- as.factor(bio$Type)
bio$Measure <- as.factor(bio$Measure)
str(bio)

# Calculate the biodiversity at the frame level
bio <- bio %>% group_by(Frame_ID, Configuration, Species, Func_group, Type, Measure) %>% 
  summarise(Count_perc=sum(Count_perc))

# Correct for the percentages
for(i in 1:nrow(bio)){
  if(bio$Measure[i] == "Percent"){
    bio$Count_perc[i] <- bio$Count_perc[i]*(9/100)
  }
}

# Separate non-native species
Ceum <- filter(bio, Species=="C_eumyota")
Wsub <- filter(bio, Species=="W_subatra")

# C eumyota cover ----------------------------------------------------

max(Ceum$Count_perc)

# BY CONFIGURATION
Ceum_avg <- Ceum %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg
# Configuration     mean      se
# Grouped         0.0975  0.0548
# Random           0.225   0.104
0.225/0.0975 # There was ~ 2.3 times more C eumyota in the random treatment 
# than in the grouped treatment

mod1 <- aov(Count_perc~Configuration, data=Ceum)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1  0.0975  0.09754    1.173   0.291
# Residuals          22  1.8299  0.08318                

# W subatra cover ----------------------------------------------------

max(Wsub$Count_perc)

# BY CONFIGURATION
Wsub_avg <- Wsub %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg
# Configuration    mean      se
# Grouped          1.47   0.355
# Random           0.72   0.127
1.47/0.72 # There was ~ 2 times more W subatra cover in the grouped treatment 
# than in the random treatment

mod1 <- aov(Count_perc~Configuration, data=Wsub)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)  
# Configuration       1   3.375    3.375    3.957  0.0593 .
# Residuals          22  18.765    0.853                 