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
bio$Panel_type <- as.factor(bio$Panel_type)
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

# Separate taxa ----------------------------------------------------------

# Mobile macroinvertebrates
mobile <- filter(bio_f, Species=="C_zizyphinum"|Species=="C_maenas"|Species=="Limpet_recruit"|Species=="L_littorea"|Species=="L_obtusata"|Species=="M_varia"|Species=="N_lapillus"|Species=="O_erinaceus"|Species=="O_bilamellata"|Species=="P_nothus"|Species=="P_quadrilineata"|Species=="Polychaete"|Species=="S_umbilicalis"|Species=="T_reticulata")

# Sessile macroinvertebrates
sessile <- filter(bio_f, Species=="Cirripedia"|Species=="C_eumyota"|Species=="O_littoralis"|Species=="S_spirorbis"|Species=="Spirobranchus"|Species=="W_subatra"|Species=="A_ephippium"|Species=="M_edulis"|Species=="M_gigas")

# Macrophytes
macrophyte <- filter(bio_f, Species=="F_serratus"|Species=="F_spiralis"|Species=="F_vesiculosus"|Species=="Juv_fucus"|Species=="U_intestinalis"|Species=="U_lactuca")

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
# Grouped         4.17  0.423
# Random          4.75  0.250

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
# Configuration       1    2.04    2.042    1.407   0.248
# Residuals          22   31.92    1.451  

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
# Grouped         4.83  0.474
# Random          5.33  0.482

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
# Configuration       1    1.50    1.500    0.547   0.467
# Residuals          22   60.33    2.742               

# Calculate mobile taxon abundance ---------------------------------

mobile_a <- mobile %>% group_by(Configuration, Frame_ID) %>% 
  summarise(n=sum(Count_perc))
mobile_a

# Summary table
mab_a_t_avg <- mobile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mab_a_t_avg
# Configuration   mean     se
# Grouped         14.3   2.15
# Random          21.3   3.45

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
# Configuration       1    294    294.00    2.962  0.0993 .
# Residuals          22   2183     99.24                 

# Calculate sessile taxon abundance ---------------------------------

sessile <- filter(sessile, Species=="Cirripedia"|Species=="C_eumyota"|Species=="O_littoralis"|Species=="S_spirorbis"|Species=="Spirobranchus"|Species=="W_subatra")

sessile_a <- sessile %>% group_by(Frame_ID, Configuration) %>% 
  summarise(n=sum(Count_perc))
sessile_a

# Summary table
sab_a_t_avg <- sessile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
sab_a_t_avg
# Configuration   mean    se
# Grouped         6.75  1.08
# Random          8.03  1.05

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
# Configuration       1    9.87    9.869    0.725   0.404
# Residuals          22  299.67   13.621               

# Calculate macrophyte taxon abundance ---------------------------------

macro_a <- macrophyte %>% group_by(Frame_ID, Configuration) %>% 
  summarise(n=sum(Count_perc))
macro_a

# Summary table
aab_a_t_avg <- macro_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
aab_a_t_avg
# Configuration   mean    se
# Grouped         61.2  3.69
# Random          51.8  3.81

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(macro_a$n, main="Density plot of macrophyte abundance", xlab="Macrophyte taxon abundance")
ggqqplot(macro_a$n, main="Q-Q plot of macrophyte taxon abundance")

# Shapiro-Wilk test
shapiro.test(macro_a$n)

# Analysis ----------------------------------------------------------------

mod1 <- aov(n ~ Configuration, data=macro_a)
summary(mod1)
#                    Df  Sum Sq  Mean Sq  F value  Pr(>F)
# Configuration       1     530    529.9    3.132   0.091
# Residuals          22    3722    169.2               
