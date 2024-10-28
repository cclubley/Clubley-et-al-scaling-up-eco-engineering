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
bio$Panel_ID <- as.factor(bio$Panel_ID)
bio$Panel <- as.factor(bio$Panel)
bio$Type <- as.factor(bio$Type)
bio$Measure <- as.factor(bio$Measure)
str(bio)

# Remove sediment
bio <- filter(bio, !Species=="Sediment")

# Calculate the total number of species detected -------------------------

# Calculate total taxon richness
total <- bio %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count()
total
# 29 taxa in total

# Which taxa were unique?
# Remove rows where species weren't present
unique <- filter(bio, Count_perc > 0)
unique <- unique %>% group_by(Configuration, Panel) %>% distinct(Species, .keep_all=TRUE)
table(unique$Species)
# M varia, O erinaceus, Polychaete, O bilamellata & P quadrilineata 
# were unique species:
# Mimachlamys varia = Random, 5 cm
# Unknown polychaete = Random, 5 cm
# Opisthobranch sp = Random, 5 cm
# Ocenebra erinaceus = Random, 2.5 cm
# Polycera quadrilineata = Random, Flat

rich <- filter(bio, Count_perc > 0)

# FOR CONFIGURATION ONLY
mosrich <- rich %>% group_by(Configuration) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
mosrich
# Configuration   n
# Grouped        23
# Random         28

# FOR PANEL ONLY
prich <- rich %>% group_by(Panel) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Panel)
prich
# Panel       n
# Flat       18
# Ripple     18
# 2.5 cm     23
# 5 cm       25

# Plot --------------------------------------------------------------------

# FOR CONFIGURATION ONLY
ggplot(data=mosrich, aes(x=Configuration, y=n, fill=Configuration))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of taxa detected")+
  xlab("Configuration")+
  ylim(c(0, 30))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# FOR PANEL ONLY
# Re-order the panel types
prich$Panel = factor(prich$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=prich, aes(x=Panel, y=n, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of taxa detected")+
  xlab("Panel")+
  ylim(c(0, 30))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Separate mobile taxa ---------------------------------------------------

mobile <- filter(bio, Type=="Mobile")
sedentary <- filter(bio, Type=="Sedentary")
table(sedentary$Species)
mobile <- rbind(mobile, sedentary)

# Separate sessile taxa ---------------------------------------------------

sessile <- filter(bio, Type=="Sessile")
table(sessile$Species)
sessile <- filter(sessile, !Species=="A_ephippium")
sessile <- filter(sessile, !Species=="M_gigas")

# Calculate total taxon richness ---------------------------------------------

# Remove rows where species weren't present
trich <- filter(bio, Count_perc > 0)

trich$Tile_ID <- 1:nrow(trich)

# Calculate taxon richness by tile
trich <- trich %>% group_by(Frame_ID, Panel_ID, Panel) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
nrow(trich) 

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(trich$n, main="Density plot of taxon richness", xlab="Taxon richness")
ggqqplot(trich$n, main="Q-Q plot of taxon richness")

# Shapiro-Wilk test
shapiro.test(trich$n)
# p = 1.935e-07
# Check the skewness
skewness(trich$n)
# 0.6350097

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(trich$n, plotit=FALSE))
# p = 7.774e-05

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Interaction between configuration and panel
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel + (1|Configuration/Frame_ID), data=trich, family=poisson)
summary(mod1)
# AIC = 860.4

# The interaction term is not significant (P = 0.09), so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.5
summary(mod2)
# AIC = 856.8 

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.3
summary(mod3)
# AIC = 856.0
check_model(mod3) 
plot(mod3) 
qqnorm(resid(mod3))
qqline(resid(mod3)) 
hist(residuals(mod3))

# plot residuals vs nominal variables
E2 <- resid(mod3)
boxplot(E2 ~ Configuration, data=trich)
abline(0,0)
boxplot(E2 ~ Panel, data=trich)
abline(0,0)
boxplot(E2 ~ Frame_ID, data=trich)
abline(0,0)

summary(mod3)
#                  Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)      1.899763    0.068981   27.541   < 2e-16 ***
# Panelcm_5       -0.004521    0.096956   -0.047     0.963    
# PanelFlat       -0.809451    0.095147   -8.507   < 2e-16 ***
# PanelRipple     -0.612016    0.090665   -6.750  1.48e-11 ***

# Post-hoc test
em <- emmeans(mod3, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate       SE   df   z.ratio   p.value
# Flat - Ripple   -0.19743   0.0920  Inf    -2.146    0.1386
# Flat - cm_25     0.80945   0.0951  Inf     8.507    <.0001 ***
# Flat - cm_5      0.80493   0.0953  Inf     8.443    <.0001 ***
# Ripple - cm_25   0.61202   0.0907  Inf     6.750    <.0001 ***
# Ripple - cm_5    0.60750   0.0909  Inf     6.686    <.0001 ***
# cm_25 - cm_5     0.00452   0.0970  Inf     0.047    1.0000

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Average taxon richness for panel
rich_a_t_avg <- trich %>% group_by(Panel) %>% 
  summarise(mean=mean(n), se=std.error(n))
rich_a_t_avg
# Tile_type   mean     se
# cm_25       6.72  0.308
# cm_5        6.67  0.378
# Flat        2.99  0.179
# Ripple      3.64  0.184

# Re-order the panel types
rich_a_t_avg$Panel = factor(rich_a_t_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=rich_a_t_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average taxon richness")+
  xlab("Panel")+
  ylim(c(0, 10))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate mobile taxon richness -----------------------------------------

# Remove rows where species weren't present
mrich <- filter(mobile, Count_perc > 0)

# Calculate taxon richness by panel
mrich <- mrich %>% group_by(Frame_ID, Panel_ID, Panel) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)

# Get the information for every tile
mrich2 <- mobile %>% group_by(Frame_ID, Panel_ID, Panel) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
# Set the richness to 0
mrich2$n <- 0
# Combine the two and remove duplicates
mrich <- rbind(mrich, mrich2)
mrich <- mrich[order(mrich$Panel_ID, -abs(mrich$n) ), ]
mrich <- mrich[ !duplicated(mrich$Panel_ID), ] 
nrow(mrich) 

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mrich$n, main="Density plot of mobile taxon richness", xlab="Mobile taxon richness")
ggqqplot(mrich$n, main="Q-Q plot of mobile taxon richness")

# Shapiro-Wilk test
shapiro.test(mrich$n)
# p = 1.547e-14
# Check the skewness
skewness(mrich$n)
# 0.97

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mrich$n, plotit=FALSE))
# p = 1.554e-13

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Interaction between configuraiton and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel + (1|Configuration/Frame_ID), data=mrich, family=poisson)
summary(mod1)
# AIC = 561.5

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.2
summary(mod2)
# AIC = 560.4

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.1
summary(mod3)
check_model(mod3) 
plot(mod3) 
qqnorm(resid(mod3))
qqline(resid(mod3)) 
hist(residuals(mod3)) 

# plot residuals vs nominal variables
E2 <- resid(mod3)
boxplot(E2 ~ Panel, data=mrich) 
abline(0,0)
boxplot(E2 ~ Frame_ID, data=mrich)
abline(0,0)

summary(mod3)
#                Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)      0.7094      0.1439    4.930  8.23e-07 ***
# Panelcm_5       -0.0981      0.2031   -0.483     0.629   
# PanelFlat       -1.3083      0.1990   -6.574  4.90e-11 ***
# PanelRipple     -0.9966      0.1820   -5.475  4.38e-08 ***

# Post-hoc test
em <- emmeans(mod3, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate      SE    df  z.ratio   p.value
# Flat - Ripple    -0.3118   0.205   Inf   -1.524    0.4233
# Flat - cm_25      1.3083   0.199   Inf    6.574    <.0001 ***
# Flat - cm_5       1.2102   0.208   Inf    5.822    <.0001 ***
# Ripple - cm_25    0.9966   0.182   Inf    5.475    <.0001 ***
# Ripple - cm_5     0.8984   0.192   Inf    4.687    <.0001 ***
# cm_25 - cm_5      0.0981   0.203   Inf    0.483    0.9629
 
# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Average taxon richness for panel type
mrich_a_t_avg <- mrich %>% group_by(Panel) %>% 
  summarise(mean=mean(n), se=std.error(n))
mrich_a_t_avg
# Panel     mean     se
# cm_25     2.22   0.20
# cm_5      1.81   0.22
# Flat      0.57   0.09
# Ripple    0.78   0.11

# Re-order the tile types
mrich_a_t_avg$Panel = factor(mrich_a_t_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=mrich_a_t_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average taxon richness")+
  xlab("Panel")+
  ggtitle("Average richness of mobile and sedentary taxa")+
  ylim(c(0, 3))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate sessile taxon richness -----------------------------------------

# Remove rows where species weren't present
srich <- filter(sessile, Count_perc > 0)

# Calculate taxon richness by tile
srich <- srich %>% group_by(Frame_ID, Panel_ID, Panel) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(srich$n, main="Density plot of sessile taxon richness", xlab="Sessile taxon richness")
ggqqplot(srich$n, main="Q-Q plot of sessile taxon richness")

# Shapiro-Wilk test
shapiro.test(srich$n)
# p = 1.967e-09
# Check the skewness
skewness(srich$n)
# 0.6481309

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(srich$n, plotit=FALSE))
# p = 7.056e-08

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel + (1|Configuration/Frame_ID), data=srich, family=poisson)
summary(mod1)
# AIC = 728.6

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.9
summary(mod2)
# AIC = 723

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.6
summary(mod3)
# AIC = 721.3
check_model(mod3)
plot(mod3) 
qqnorm(resid(mod3))
qqline(resid(mod3))
hist(residuals(mod3)) 

# plot residuals vs nominal variables
E2 <- resid(mod3)
boxplot(E2 ~ Panel, data=mrich) 
abline(0,0)
boxplot(E2 ~ Frame_ID, data=mrich)
abline(0,0)

summary(mod3)
#               Estimate   Std. Error   z value   Pr(>|z|)    
# (Intercept)    1.47908      0.07956    18.592    < 2e-16 ***
# Panelcm_5      0.04335      0.11131     0.389      0.697    
# PanelFlat     -0.68057      0.11216    -6.068   1.30e-09 ***
# PanelRipple   -0.46748      0.10668    -4.382   1.17e-05 *** 


# Post-hoc test
em <- emmeans(mod3, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple  -0.2131    0.106   Inf  -2.005   0.1862
# Flat - cm_25    0.6806    0.112   Inf   6.068   <.0001 ***
# Flat - cm_5     0.7239    0.111   Inf   6.525   <.0001 ***
# Ripple - cm_25  0.4675    0.107   Inf   4.382   0.0001 **
# Ripple - cm_5   0.5108    0.105   Inf   4.846   <.0001 ***
# cm_25 - cm_5    -0.0434   0.111   Inf  -0.389   0.9799

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Average taxon richness for panel type
srich_a_t_avg <- srich %>% group_by(Panel) %>% 
  summarise(mean=mean(n), se=std.error(n))
srich_a_t_avg
# Panel     mean     se
# cm_25     4.39  0.200
# cm_5      4.58  0.212
# Flat      2.22  0.111
# Ripple    2.75  0.106

# Re-order the panel types
srich_a_t_avg$Panel = factor(srich_a_t_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=srich_a_t_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average taxon richness")+
  xlab("Panel")+
  ggtitle("Average richness of sessile taxa")+
  ylim(c(0, 5))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate mobile taxon abundance ---------------------------------

mobile_a <- mobile %>% group_by(Frame_ID, Configuration, Panel_ID, Panel) %>% 
  summarise(n=sum(Count_perc))
mobile_a

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mobile_a$n, main="Density plot of mobile abundance", xlab="Mobile taxon abundance")
ggqqplot(mobile_a$n, main="Q-Q plot of mobile taxon abundance")

# Shapiro-Wilk test
shapiro.test(mobile_a$n)
# p = < 2.2e-16
# Check the skewness
skewness(mobile_a$n)
# 2.165933

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mobile_a$n, plotit=FALSE))
# p = 9.289e-13

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel + (1|Configuration/Frame_ID), data=mobile_a, family=poisson)
summary(mod1)
# AIC = 834

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.2
summary(mod2)
# AIC = 833
# Configuration is marginally not significant (p = 0.05)
check_model(mod2) 
plot(mod2) 
qqnorm(resid(mod2))
qqline(resid(mod2)) 
hist(residuals(mod2))

# plot residuals vs nominal variables
E2 <- resid(mod2)
boxplot(E2 ~ Panel, data=mobile_a) 
abline(0,0)
boxplot(E2 ~ Frame_ID, data=mobile_a)
abline(0,0)

summary(mod2)
#                      Estimate  Std. Error  z value   Pr(>|z|)    
# (Intercept)           1.02863     0.19964    5.152   2.57e-07 ***
# ConfigurationRandom   0.39212     0.23823    1.646     0.0998 .  
# Panelcm_5             0.05302     0.18626    0.285     0.7759 
# PanelFlat            -1.45879     0.15470   -9.430    < 2e-16 ***
# PanelRipple          -1.25400     0.14552   -8.617    < 2e-16 ***

# Post-hoc test
em <- emmeans(mod2, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple  -0.205     0.164   Inf  -1.246    0.5975
# Flat - cm_25    1.459     0.155   Inf   9.430    <.0001 ***
# Flat - cm_5     1.512     0.175   Inf   8.642    <.0001 ***
# Ripple - cm_25  1.254     0.146   Inf   8.617    <.0001 ***
# Ripple - cm_5   1.307     0.167   Inf   7.832    <.0001 ***
# cm_25 - cm_5   -0.053     0.186   Inf  -0.285    0.9920

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot ------------------------------------------------------------------------

mobile_a_avg <- mobile_a %>% group_by(Panel) %>% 
  summarise(mean=mean(n), se=std.error(n))
mobile_a_avg
# Panel    mean     se
# cm_25    4.75  0.628
# cm_5     3.31  0.547
# Flat    0.917  0.219
# Ripple   1.12  0.194

# Re-order the tile types
mobile_a_avg$Panel = factor(mobile_a_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=mobile_a_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average mobile taxon abundance")+
  xlab("Panel")+
  ylim(c(0, 6))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

mobile_a_avg_m <- mobile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mobile_a_avg_m
# Mosaic   mean     se
# Grouped  1.62  0.233
# Random   2.43  0.311

ggplot(data=mobile_a_avg_m, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average mobile taxon abundance")+
  xlab("Configuration")+
  ylim(c(0, 6))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Calculate sessile taxon abundance ---------------------------------

sessile_a <- sessile %>% group_by(Frame_ID, Configuration, Panel_ID, Panel) %>% 
  summarise(n=sum(Count_perc))
sessile_a

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(sessile_a$n, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggqqplot(sessile_a$n, main="Q-Q plot of sessile taxon abundance")

# Shapiro-Wilk test
shapiro.test(sessile_a$n)
# p = 0.00231
# Check the skewness
skewness(sessile_a$n)
# 0.2462215

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(sessile_a$n, plotit=FALSE))
# p = 0.01046

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Use a beta family because the values are proportions/percentages

# Normalise the data between 0 and 1 (but actuallt 0.01 and 0.99 because glmmTMB
# with a beta family won't take absolute 0 and 1 values)
sessile_a$n2 <- rescale(sessile_a$n, to=c(0.01, 0.99))

# All the transformation does is change the scale
ggdensity(sessile_a$n, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggdensity(sessile_a$n2, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggqqplot(sessile_a$n, main="Q-Q plot of sessile taxon abundance")
ggqqplot(sessile_a$n2, main="Q-Q plot of sessile taxon abundance")

set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration*Panel + (1|Configuration/Frame_ID), data=sessile_a, family=beta_family())
summary(mod1)
# AIC = -153                 

set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
anova(mod1, mod2) # P = 0.6
summary(mod2)
check_model(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2)) 
hist(residuals(mod2))

# plot residuals vs nominal variables
E2 <- resid(mod2)
boxplot(E2 ~ Panel, data=sessile_a) 
abline(0,0)
boxplot(E2 ~ Configuration, data=sessile_a) 
abline(0,0)

summary(mod2)
#                     Estimate  Std. Error  z value   Pr(>|z|)    
# (Intercept)           0.1279      0.1500    0.853     0.3936   
# ConfigurationRandom  -0.2746      0.1381   -1.989     0.0467 *  
# Panelcm_5             0.1955      0.1854    1.054     0.2917      
# PanelFlat            -0.9580      0.1524   -6.285   3.27e-10 ***
# PanelRipple          -0.5869      0.1498   -3.918   8.95e-05 ***

# Post-hoc test
em <- emmeans(mod3, ~Panel)
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate     SE   df z.ratio p.value
# Flat - Ripple    -0.2131  0.106  Inf  -2.005  0.1862
# Flat - cm_25      0.6806  0.112  Inf   6.068  <.0001 ***
# Flat - cm_5       0.7239  0.111  Inf   6.525  <.0001 ***
# Ripple - cm_25    0.4675  0.107  Inf   4.382  0.0001 ***
# Ripple - cm_5     0.5108  0.105  Inf   4.846  <.0001 ***
# cm_25 - cm_5     -0.0434  0.111  Inf  -0.389  0.9799

cld(em)
# Flat    a
# Ripple  b
# 2.5 cm  c
# 5 cm    c

# Plot ------------------------------------------------------------------------

sessile_a_avg <- sessile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
sessile_a_avg
# Configuration   mean     se
# Grouped         83.9   3.40
# Random          73.9   3.97

ggplot(data=sessile_a_avg, aes(x=Configuration, y=mean, fill=Configuration))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sessile taxon cover (%)")+
  xlab("Configuration")+
  scale_y_continuous(limits=c(0, 125), breaks=seq(0, 125, 25))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

sessile_a_avg <- sessile_a %>% group_by(Panel) %>% 
  summarise(mean=mean(n), se=std.error(n))
sessile_a_avg
# Panel      mean     se
# Flat       58.3   3.68
# Ripple     72.6   3.45
# 2.5 cm    101.0   7.20
# 5 cm      110.0   5.54

# Re-order the panel types
sessile_a_avg$Panel = factor(sessile_a_avg$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=sessile_a_avg, aes(x=Panel, y=mean, fill=Panel))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sessile taxon abundance")+
  xlab("Panel")+
  scale_y_continuous(limits=c(0, 125), breaks=seq(0, 125, 25))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Plot the average count of mobile species groups -------------------------

# Calculate average abundance by tile
mobileavgt <- mobile %>% group_by(Configuration, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))

mobileavgt <- filter(mobileavgt, !Func_group=="NA")

ggplot(data=mobileavgt, aes(x=Configuration, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (count)")+
  xlab("Configuration")+
  ylim(c(0, 1.5))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate average abundance by panel
mobileavgt <- mobile %>% group_by(Panel, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
mobileavgt <- filter(mobileavgt, !Func_group=="NA")

# Re-order the panel types
mobileavgt$Panel = factor(mobileavgt$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=mobileavgt, aes(x=Panel, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (count)")+
  xlab("Panel")+
  ylim(c(0, 1.5))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Plot the average cover of sessile species groups -------------------------

# Calculate average abundance by tile
sessileavgm <- sessile %>% group_by(Configuration, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))

sessileavgm <- filter(sessileavgm, !Func_group=="NA")

ggplot(data=sessileavgm, aes(x=Configuration, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.5, end=0.9)+
  ylab("Average cover (%)")+
  xlab("Configuration")+
  ylim(c(0, 40))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.12, 0.9))

# Calculate average abundance by tile
sessileavgt <- sessile %>% group_by(Panel, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
sessileavgt <- filter(sessileavgt, !Func_group=="NA")
# Re-order the panel types
sessileavgt$Panel = factor(sessileavgt$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

ggplot(data=sessileavgt, aes(x=Panel, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.5, end=0.9)+
  ylab("Average cover (%)")+
  xlab("Panel")+
  ylim(c(0, 40))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.12, 0.9))
