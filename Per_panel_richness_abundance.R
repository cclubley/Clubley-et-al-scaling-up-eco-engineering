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
# Mimachlamys varia = Random, 5 cm
# Unknown polychaete = Random, 5 cm
# Onchidoris bilamellata = Random, 5 cm
# Ocenebra erinaceus = Random, 2.5 cm
# Polycera quadrilineata = Random, Flat

rich <- filter(bio, Count_perc > 0)

# Summary table - configuration
mosrich <- rich %>% group_by(Configuration) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
mosrich
# Configuration   n
# Grouped        23
# Random         28

# Summary table - panel
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

# Configuration
ggplot(data=mosrich, aes(x=Configuration, y=n))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of taxa detected")+
  xlab("Configuration")+
  ylim(c(0, 30))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Panel
prich$Panel = factor(prich$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=prich, aes(x=Panel, y=n))+
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
mobile <- rbind(mobile, sedentary)

# Separate sessile taxa ---------------------------------------------------

sessile <- filter(bio, Type=="Sessile")
sessile <- filter(sessile, !Species=="A_ephippium")
sessile <- filter(sessile, !Species=="M_gigas")

# Calculate total taxon richness ---------------------------------------------

# Remove rows where species weren't present
trich <- filter(bio, Count_perc > 0)
trich$Tile_ID <- 1:nrow(trich)

# Calculate taxon richness by panel
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

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(trich$n, plotit=FALSE))

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)

# Interaction between configuration and panel
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel + (1|Configuration/Frame_ID), data=trich, family=poisson)
summary(mod1)
# AIC = 860.4

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.5
summary(mod2)
# AIC = 856.8 
check_model(mod2) 
plot(mod2) 
qqnorm(resid(mod2))
qqline(resid(mod2)) 
hist(residuals(mod2))

# plot residuals vs nominal variables
E2 <- resid(mod2)
boxplot(E2 ~ Configuration, data=trich)
abline(0,0)
boxplot(E2 ~ Panel, data=trich)
abline(0,0)
boxplot(E2 ~ Frame_ID, data=trich)
abline(0,0)

summary(mod2)
#                      Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)          1.861179    0.077773   23.931   < 2e-16 ***
# ConfigurationRandom  0.077303    0.071944    1.074     0.283    
# Panelcm_5           -0.004572    0.095978   -0.048     0.962    
# PanelFlat           -0.809591    0.094909   -8.530   < 2e-16 ***
# PanelRipple         -0.611883    0.090408   -6.768  1.31e-11 ***

# Post-hoc test
em <- emmeans(mod2, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate       SE   df   z.ratio   p.value
# Flat - Ripple   -0.19771   0.0920  Inf    -2.149    0.1378
# Flat - 2.5 cm    0.80959   0.0949  Inf     8.530    <.0001
# Flat - 5 cm      0.80502   0.0951  Inf     8.465    <.0001
# Ripple - 2.5 cm  0.61188   0.0904  Inf     6.768    <.0001
# Ripple - 5 cm    0.60731   0.0906  Inf     6.703    <.0001
# 2.5 cm - 5 cm    0.00457   0.0960  Inf     0.048    1.0000

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
ggplot(data=rich_a_t_avg, aes(x=Panel, y=mean))+
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

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mrich$n, plotit=FALSE))

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
check_model(mod2) 
plot(mod2) 
qqnorm(resid(mod2))
qqline(resid(mod2)) 
hist(residuals(mod2)) 

# plot residuals vs nominal variables
E2 <- resid(mod2)
boxplot(E2 ~ Panel, data=mrich) 
abline(0,0)
boxplot(E2 ~ Frame_ID, data=mrich)
abline(0,0)

summary(mod2)
#                     Estimate  Std. Error  z value  Pr(>|z|)  
# (Intercept)           0.5941      0.1661    3.576  0.000349 ***
# ConfigurationRandom   0.2455      0.1615    1.521  0.128351    
# Panelcm_5            -0.1111      0.1990   -0.558  0.576747    
# PanelFlat            -1.3145      0.1984   -6.626  3.45e-11 ***
# PanelRipple          -1.0026      0.1813   -5.531  3.19e-08 ***

# Post-hoc test
em <- emmeans(mod2, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate      SE    df  z.ratio   p.value
# Flat - Ripple     -0.312   0.205   Inf   -1.522    0.4243
# Flat - 2.5 cm      1.314   0.198   Inf    6.626    <.0001
# Flat - 5 cm        1.203   0.207   Inf    5.817    <.0001
# Ripple - 2.5 cm    1.003   0.181   Inf    5.531    <.0001
# Ripple - 5 cm      0.892   0.191   Inf    4.679    <.0001
# 2.5 cm - 5 cm      0.111   0.199   Inf    0.558    0.9444

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
ggplot(data=mrich_a_t_avg, aes(x=Panel, y=mean))+
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

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(srich$n, plotit=FALSE))

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
check_model(mod2)
plot(mod2) 
qqnorm(resid(mod2))
qqline(resid(mod2))
hist(residuals(mod2)) 

# plot residuals vs nominal variables
E2 <- resid(mod2)
boxplot(E2 ~ Panel, data=mrich) 
abline(0,0)
boxplot(E2 ~ Frame_ID, data=mrich)
abline(0,0)

summary(mod2)
#                     Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept)          1.45980     0.08863   16.471   < 2e-16 ***
# ConfigurationRandom  0.03818     0.07665    0.498     0.618    
# Panelcm_5            0.04335     0.11131    0.389     0.697    
# PanelFlat           -0.68057     0.11216   -6.068  1.29e-09 ***
# PanelRipple         -0.46748     0.10668   -4.382  1.17e-05 ***

# Post-hoc test
em <- emmeans(mod2, "Panel")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate     SE    df  z.ratio  p.value
# Flat - Ripple    -0.2131  0.106   Inf   -2.005   0.1862
# Flat - cm_25      0.6806  0.112   Inf    6.068   <.0001
# Flat - cm_5       0.7239  0.111   Inf    6.525   <.0001
# Ripple - cm_25    0.4675  0.107   Inf    4.382   0.0001
# Ripple - cm_5     0.5108  0.105   Inf    4.846   <.0001
# cm_25 - cm_5     -0.0434  0.111   Inf   -0.389   0.9799

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
ggplot(data=srich_a_t_avg, aes(x=Panel, y=mean))+
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

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mobile_a$n, plotit=FALSE))

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
ggplot(data=mobile_a_avg, aes(x=Panel, y=mean))+
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

ggplot(data=mobile_a_avg_m, aes(x=Configuration, y=mean))+
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

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(sessile_a$n, plotit=FALSE))

# Analysis ----------------------------------------------------------------

# Generalised linear mixed models (GLMMs)
# Use a beta family because the values are proportions/percentages

# Normalise the data between 0 and 1 (but actuallt 0.01 and 0.99 because glmmTMB
# with a beta family won't take absolute 0 and 1 values)
sessile_a$n2 <- rescale(sessile_a$n, to=c(0.01, 0.99))

set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration*Panel + (1|Configuration/Frame_ID), data=sessile_a, family=beta_family())
summary(mod1)
# AIC = -153                 

# Interaction not significant so try removing
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
em <- emmeans(mod2, ~Panel)
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate     SE   df z.ratio p.value
#  Flat - Ripple    -0.371  0.122  Inf  -3.031  0.0130
# Flat - cm_25       0.958  0.152  Inf   6.285  <.0001
# Flat - cm_5        1.154  0.153  Inf   7.536  <.0001
# Ripple - cm_25     0.587  0.150  Inf   3.918  0.0005
# Ripple - cm_5      0.782  0.152  Inf   5.159  <.0001
# cm_25 - cm_5      -0.196  0.185  Inf  -1.054  0.7172

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

ggplot(data=sessile_a_avg, aes(x=Configuration, y=mean))+
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
ggplot(data=sessile_a_avg, aes(x=Panel, y=mean))+
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
