# Load required packages --------------------------------------------------

{ library(dplyr)
  library(ggplot2)
  library(ccplot)
  library(plotrix)
  library(ggpubr)
  library(moments)
  library(pracma)
  library(rcompanion)
  library(car)
  library(multcomp)
  library(emmeans)
  library(performance)
  library(viridis)
  library(blmeco)
  library(lme4)
  library(scales)
  library(glmmTMB)
}

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

# Remove sediment
bio <- filter(bio, !Species=="Sediment")

# Calculate the biodiversity at the frame level
bio_f <- bio %>% group_by(Frame_ID, Configuration, Species, Func_group, Type, Measure) %>% 
  summarise(Count_perc=sum(Count_perc))

# Correct for the percentages
for(i in 1:nrow(bio_f)){
  if(bio_f$Measure[i] == "Percent"){
    bio_f$Count_perc[i] <- bio_f$Count_perc[i]*(9/100)
  }
}

rm(bio)

# Calculate the total number of species detected -------------------------

# Which taxa were unique?
# Remove rows where species weren't present
unique <- filter(bio_f, Count_perc > 0)
unique <- unique %>% group_by(Mosaic) %>% distinct(Species, .keep_all=TRUE)
table(unique$Species)
# Calliostoma zizyphinum - grouped
# Mimachlamys varia - random
# Unknown polychaete - random
# Ocenebra erinaceus - random
# Onchidoris bilamellata - random
# Polycera quadrilineata - random
# Tritia reticulata - random

rich <- filter(bio_f, Count_perc > 0)

# Separate mobile taxa ---------------------------------------------------

mobile <- filter(bio_f, Type=="Mobile")
sedentary <- filter(bio_f, Type=="Sedentary")
mobile <- rbind(mobile, sedentary)

# Separate sessile taxa ---------------------------------------------------

sessile <- filter(bio_f, Type=="Sessile")
sessile <- filter(sessile, !Species=="A_ephippium")
sessile <- filter(sessile, !Species=="M_gigas")

# Calculate total taxon richness ---------------------------------------------

# Remove rows where species weren't present
trich <- filter(bio_f, Count_perc > 0)

# Calculate taxon richness by tile
trich <- trich %>% group_by(Frame_ID) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
nrow(trich)

# Summary table
tavg <- trich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
tavg
# Configuration   mean     se
# Grouped         13.5  0.744
# Random          14.8  0.494

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(trich$n, main="Density plot of taxon richness", xlab="Taxon richness")
ggqqplot(trich$n, main="Q-Q plot of taxon richness")

# Shapiro-Wilk test
shapiro.test(trich$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=trich)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1    9.38    9.375     1.96   0.175
# Residuals          22  105.25    4.784                          

# Calculate mobile taxon richness -----------------------------------------

# Remove rows where species weren't present
mrich <- filter(mobile, Count_perc > 0)

# Calculate taxon richness by frame
mrich <- mrich %>% group_by(Frame_ID) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)

# Summary table
mrich_a_t_avg <- mrich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mrich_a_t_avg
# Configuration   mean     se
# Grouped         4.42  0.452
# Random          5.17  0.297

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mrich$n, main="Density plot of mobile taxon richness", xlab="Mobile taxon richness")
ggqqplot(mrich$n, main="Q-Q plot of mobile taxon richness")

# Shapiro-Wilk test
shapiro.test(mrich$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=mrich)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1    3.38    3.375    1.924   0.179
# Residuals          22   38.58    1.754  

# Calculate sessile taxon richness -----------------------------------------

# Remove rows where species weren't present
srich <- filter(sessile, Count_perc > 0)

# Calculate taxon richness by tile
srich <- srich %>% group_by(Frame_ID) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)

# Summary table
srich_a_t_avg <- srich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
srich_a_t_avg
# Configuration   mean     se
# Grouped         8.00  0.477
# Random          8.75  0.351

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(srich$n, main="Density plot of sessile taxon richness", xlab="Sessile taxon richness")
ggqqplot(srich$n, main="Q-Q plot of sessile taxon richness")

# Shapiro-Wilk test
shapiro.test(srich$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=srich)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1    3.37    3.375    1.605   0.218
# Residuals          22   46.25    2.102               

# Calculate mobile taxon abundance ---------------------------------

mobile_a <- mobile %>% group_by(Configuration, Frame_ID) %>% 
  summarise(n=sum(Count_perc))
mobile_a

# Summary table
mab_a_t_avg <- mobile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mab_a_t_avg
# Configuration   mean     se
# Grouped         14.6   2.19
# Random          21.8   3.59

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mobile_a$n, main="Density plot of mobile abundance", xlab="Mobile taxon abundance")
ggqqplot(mobile_a$n, main="Q-Q plot of mobile taxon abundance")

# Shapiro-Wilk test
shapiro.test(mobile_a$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=mobile_a)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)  
# Configuration       1   315.4    315.4    2.974  0.0986 .
# Residuals          22  2332.6    106.0                 

# Calculate sessile taxon abundance ---------------------------------

sessile_a <- sessile %>% group_by(Frame_ID, Configuration) %>% 
  summarise(n=sum(Count_perc))
sessile_a

# Summary table
sab_a_t_avg <- sessile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
sab_a_t_avg
# Configuration   mean    se
# Grouped         68.0  3.09
# Random          59.9  3.61

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(sessile_a$n, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggqqplot(sessile_a$n, main="Q-Q plot of sessile taxon abundance")

# Shapiro-Wilk test
shapiro.test(sessile_a$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=sessile_a)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1   395.1    395.1    2.919   0.102
# Residuals          22  2977.7    135.4               
