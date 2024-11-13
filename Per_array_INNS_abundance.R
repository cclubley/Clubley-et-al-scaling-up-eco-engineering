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
  library(glmmTMB)
  library(scales)
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

#Summary table
Mg_avg <- C_Mgigas %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count), se=std.error(Count))
Mg_avg
# Configuration   mean     se
# Grouped         1.58  0.543
# Random          0.75  0.351

# Check for normality
ggdensity(C_Mgigas$Count)
ggqqplot(C_Mgigas$Count)

# Shapiro-Wilk test
shapiro.test(C_Mgigas$Count)
# p = 6.774e-05

# Check for zero-inflation
prop_zeros <- mean(C_Mgigas$Count == 0)
prop_zeros 

set.seed(1234)
mod1 <- glmmTMB(Count ~ Configuration, ziformula=~1, data=C_Mgigas, family=poisson)
check_model(mod1)

summary(mod1)
# Fixed effects:
#                         Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)               0.7931      0.2496    3.177   0.00149 **
# ConfigurationRandom      -0.4149      0.5322   -0.780   0.43559  

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

# Summary table
Ceum_avg <- Ceum %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Ceum_avg
# Configuration    mean    se
# Grouped        0.0975  0.0548
# Random          0.225   0.104

# Check for normality
ggdensity(Ceum$Count_perc)
ggqqplot(Ceum$Count_perc)

# Shapiro-Wilk test
shapiro.test(Ceum$Count_perc)
# p = 1.221e-06

# Check for zero-inflation
prop_zeros <- mean(Ceum$Count_perc == 0)
prop_zeros

Ceum$n2 <- rescale(Ceum$Count_perc, to=c(0.01, 0.99))
set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration, ziformula=~1, data=Ceum, family=beta_family())
check_model(mod1)

summary(mod1)
# Fixed effects:
#                       Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)            -1.5696      0.3765   -4.169  3.07e-05 ***
# ConfigurationRandom     0.4297      0.4577    0.939     0.348     

# W subatra cover ----------------------------------------------------

# Summary table
Wsub_avg <- Wsub %>% group_by(Configuration) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
Wsub_avg
# Configuration    mean      se
# Grouped          1.47   0.355
# Random           0.72   0.127

# Check for normality
ggdensity(Wsub$Count_perc)
ggqqplot(Wsub$Count_perc)

# Shapiro-Wilk test
shapiro.test(Wsub$Count_perc)
# p = 0.004214

Wsub$n2 <- rescale(Wsub$Count_perc, to=c(0.01, 0.99))
set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration, data=Wsub, family=beta_family())
check_model(mod1)

summary(mod1)
# Fixed effects:
#                       Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)            -0.2926      0.3132   -0.934     0.350
# ConfigurationRandom    -0.5864      0.4459   -1.315     0.188