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
bio$Panel_position <- as.factor(bio$Panel_position)
bio$Panel_type <- as.factor(bio$Panel_type)
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
rm(total)

# Which taxa were unique?
# Remove rows where species weren't present
unique <- filter(bio, Count_perc > 0)
unique <- unique %>% group_by(Configuration, Panel_type) %>% distinct(Species, .keep_all=TRUE)
table(unique$Species)
rm(unique)
# M varia, O erinaceus, Unknown polychaete, O bilamellata & P quadrilineata 
# were unique species:
# Mimachlamys varia = Random, 5 cm
# Unknown polychaete = Random, 5 cm
# Onchidoris bilamellata = Random, 5 cm
# Ocenebra erinaceus = Random, 2.5 cm
# Polycera quadrilineata = Random, Flat

# Species richness
rich <- filter(bio, Count_perc > 0)

# Configuration
conrich <- rich %>% group_by(Configuration) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
conrich
# Configuration   n
# Grouped        23
# Random         28

# Panel type
panrich <- rich %>% group_by(Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Panel_type)
panrich
# Panel_type   n
# Flat        18
# Ripple      18
# 2.5 cm      23
# 5 cm        25

rm(rich)

# Plot --------------------------------------------------------------------

# Configuration
ggplot(data=conrich, aes(x=Configuration, y=n))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of taxa detected")+
  xlab("Configuration")+
  ylim(c(0, 30))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")
rm(conrich)

# Panel type
# Re-order the panel types
panrich$Panel_type = factor(panrich$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=panrich, aes(x=Panel_type, y=n))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Number of taxa detected")+
  xlab("Panel type")+
  ylim(c(0, 30))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")
rm(panrich)

# Separate taxa ----------------------------------------------------------

# Mobile macroinvertebrates
mobile <- filter(bio, Species=="C_zizyphinum"|Species=="C_maenas"|Species=="Limpet_recruit"|Species=="L_littorea"|Species=="L_obtusata"|Species=="M_varia"|Species=="N_lapillus"|Species=="O_erinaceus"|Species=="O_bilamellata"|Species=="P_nothus"|Species=="P_quadrilineata"|Species=="Polychaete"|Species=="S_umbilicalis"|Species=="T_reticulata")

# Sessile macroinvertebrates
sessile <- filter(bio, Species=="Cirripedia"|Species=="C_eumyota"|Species=="O_littoralis"|Species=="S_spirorbis"|Species=="Spirobranchus"|Species=="W_subatra"|Species=="A_ephippium"|Species=="M_edulis"|Species=="M_gigas")

# Macrophytes
macrophyte <- filter(bio, Species=="F_serratus"|Species=="F_spiralis"|Species=="F_vesiculosus"|Species=="Juv_fucus"|Species=="U_intestinalis"|Species=="U_lactuca")

# Calculate total taxon richness ---------------------------------------------

# Remove rows where species weren't present
prich <- filter(bio, Count_perc > 0)

# Calculate taxon richness by tile
prich <- prich %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
nrow(prich)

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(prich$n, main="Density plot of taxon richness", xlab="Taxon richness")
ggqqplot(prich$n, main="Q-Q plot of taxon richness")

# Shapiro-Wilk test
shapiro.test(prich$n)
# p = 1.935e-07

# Try data transformation
# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(prich$n, plotit=FALSE))
# p = 7.774e-05

# Analysis ----------------------------------------------------------------

# As the data is non-normal, use a generalised linear mixed effects model. Configuration and Panel type are fixed factors, whilst frame ID is a random factor nested in Configuration

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=prich, family=poisson)
summary(mod1)
# AIC = 860.4

# The interaction term is not significant, so try model reduction 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.5
summary(mod2)
# AIC = 856.8 

# Configuration is not significant, so simplify further
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.3
summary(mod3)
# AIC = 856.0
check_model(mod3)
# plot(mod3) 
# qqnorm(resid(mod3))
# qqline(resid(mod3)) 
# hist(residuals(mod3))
# E2 <- resid(mod3)
# boxplot(E2 ~ Configuration, data=prich)
# abline(0,0)
# boxplot(E2 ~ Panel_type, data=prich)
# abline(0,0)
# boxplot(E2 ~ Frame_ID, data=prich)
# abline(0,0)

summary(mod3)
# AIC      BIC     logLik    deviance   df.resid 
# 856.0    876.2   -422.0    844.0      210

# Fixed effects:
#                  Estimate    Std. Error  z value  Pr(>|z|)    
# (Intercept)       1.899763   0.068981    27.541   < 2e-16 ***
# Tile_typecm_5    -0.004521   0.096956    -0.047   0.963    
# Tile_typeFlat    -0.809451   0.095147    -8.507   < 2e-16 ***
# Tile_typeRipple  -0.612016   0.090665    -6.750   1.48e-11 ***

# Post-hoc test
em <- emmeans(mod3, "Panel_type")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple  -0.19743   0.0920  Inf  -2.146    0.1386
# Flat - cm_25    0.80945   0.0951  Inf   8.507    <.0001 ***
# Flat - cm_5     0.80493   0.0953  Inf   8.443    <.0001 ***
# Ripple - cm_25  0.61202   0.0907  Inf   6.750    <.0001 ***
# Ripple - cm_5   0.60750   0.0909  Inf   6.686    <.0001 ***
# cm_25 - cm_5    0.00452   0.0970  Inf   0.047    1.0000

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Averge taxon richness for configuration
rich_a_c_avg <- prich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
rich_a_c_avg
# Configuration  mean    se
# Grouped        4.27 0.221
# Random         4.61 0.233

# Create a plot using only the panel type, as the configuration was non-significant
# Average taxon richness for tile type
rich_a_p_avg <- prich %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
rich_a_p_avg
# Panel_type   mean     se
# cm_25       6.72  0.308
# cm_5        6.67  0.378
# Flat        2.99  0.179
# Ripple      3.64  0.184

# Re-order the panel types
rich_a_p_avg$Panel_type = factor(rich_a_p_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=rich_a_p_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average taxon richness")+
  xlab("Panel type")+
  ylim(c(0, 10))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

rm(mod1, mod2, mod3, rich_a_p_avg, E2, em)

# Calculate mobile taxon richness -----------------------------------------

# Remove rows where species weren't present
mrich <- filter(mobile, Count_perc > 0)

# Calculate taxon richness by tile
mrich <- mrich %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
# There should be 216 rows - one for each tile (24 x 9), even if the taxon 
# richness is 0
# mrich is only 130 rows, so I need to add back in the missing rows
# Get the information for every tile
mrich2 <- mobile %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
# Set the richness to 0
mrich2$n <- 0
# Combine the two and remove duplicates
mrich <- rbind(mrich, mrich2)
mrich <- mrich[order(mrich$Panel_ID, -abs(mrich$n) ), ]
mrich <- mrich[ !duplicated(mrich$Panel_ID), ] 
nrow(mrich) # 216 - perfect

rm(mrich2)

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mrich$n, main="Density plot of mobile taxon richness", xlab="Mobile taxon richness")
ggqqplot(mrich$n, main="Q-Q plot of mobile taxon richness")

# Shapiro-Wilk test
shapiro.test(mrich$n)
# p = 1.169e-14

# Try data transformation
# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mrich$n, plotit=FALSE))
# p = 6.827e-14
# The data can't be transformed

# Analysis ----------------------------------------------------------------

# Interaction between mosaic and tile type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=mrich, family=poisson)
summary(mod1)
# AIC = 552.9    

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.15
summary(mod2)
# AIC = 552.3     

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.16
summary(mod3)
# AIC = 552.2     
check_model(mod3)
# plot(mod3)
# qqnorm(resid(mod3))
# qqline(resid(mod3)) 
# hist(residuals(mod3)) 
# E2 <- resid(mod3)
# boxplot(E2 ~ Panel_type, data=mrich) 
# abline(0,0)
# boxplot(E2 ~ Frame_ID, data=mrich)
# abline(0,0)

summary(mod3)
# AIC      BIC     logLik    deviance   df.resid 
# 552.2    572.5   -270.1    540.2      210 

# Fixed effects:
#                  Estimate    Std. Error  z value  Pr(>|z|)    
# (Intercept)       0.63799    0.14595     4.371    1.24e-05 ***
# Panel_typecm_5   -0.06047    0.20509    -0.295       0.768    
# Panel_typeFlat   -1.25936    0.20277    -6.211    5.27e-10 ***
# Panel_typeRipple -0.92290    0.18447    -5.003    5.65e-07 ***

# Post-hoc test
em <- emmeans(mod3, "Panel_type")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple    -0.3365  0.206   Inf    -1.632   0.3606
# Flat - cm_25      1.2594  0.203   Inf     6.211   <.0001
# Flat - cm_5       1.1989  0.210   Inf     5.708   <.0001
# Ripple - cm_25    0.9229  0.184   Inf     5.003   <.0001
# Ripple - cm_5     0.8624  0.192   Inf     4.482   <.0001
# cm_25 - cm_5      0.0605  0.205   Inf     0.295   0.9911
 
# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Averge taxon richness for configuration
mrich_a_c_avg <- mrich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mrich_a_c_avg
# Configuration  mean     se
# Grouped        0.954 0.108
# Random         1.20  0.111

# Create a plot using only the panel type, as the configuration was non-significant
# Average taxon richness for tile type
mrich_a_p_avg <- mrich %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
mrich_a_p_avg
# Tile_type   mean     se
# cm_25       2.06  0.169 
# cm_5        1.75  0.216 
# Flat        0.56  0.095
# Ripple      0.78  0.112 

# Re-order the panel types
mrich_a_p_avg$Panel_type = factor(mrich_a_p_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=mrich_a_p_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average mobile macroinvertebrate richness (count)")+
  xlab("Configuration")+
  ylim(c(0, 2.5))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

rm(mod1, mod2, mod3, mrich_a_p_avg, E2, em)

# Calculate sessile taxon richness -----------------------------------------

# Remove rows where species weren't present
srich <- filter(sessile, Count_perc > 0)

# Calculate taxon richness by tile
srich <- srich %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
nrow(srich)
# Get the information for every tile
srich2 <- sessile %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
# Set the richness to 0
srich2$n <- 0
# Combine the two and remove duplicates
srich <- rbind(srich, srich2)
srich <- srich[order(srich$Panel_ID, -abs(srich$n) ), ]
srich <- srich[ !duplicated(srich$Panel_ID), ] 
nrow(srich) # 216 

rm(srich2)

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(srich$n, main="Density plot of sessile taxon richness", xlab="Sessile taxon richness")
ggqqplot(srich$n, main="Q-Q plot of sessile taxon richness")

# Shapiro-Wilk test
shapiro.test(srich$n)
# p = 1.233e-14

# Try data transformation
# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(srich$n, plotit=FALSE))
# p = 1.014e-11
# The data can't be transformed

# Analysis ----------------------------------------------------------------

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=srich, family=poisson)
summary(mod1)
# AIC = 602.6    

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.69
summary(mod2)
# AIC = 598.1            

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 1.0
summary(mod3)
# AIC = 596.1        
check_model(mod3)
# plot(mod3) 
# qqnorm(resid(mod3))
# qqline(resid(mod3)) 
# hist(residuals(mod3))
# E2 <- resid(mod3)
# boxplot(E2 ~ Panel_type, data=srich) 
# abline(0,0)
# boxplot(E2 ~ Frame_ID, data=srich)
# abline(0,0)

summary(mod3)
# AIC      BIC     logLik    deviance   df.resid 
# 596.1    616.3   -292.0    584.1      210 

# Fixed effects:
#                  Estimate    Std. Error  z value  Pr(>|z|)    
# (Intercept)        0.8629     0.1132     7.625    2.44e-14 ***
# Panel_typecm_5     0.2288     0.1557     1.469       0.142    
# Panel_typeFlat    -0.9263     0.1652    -5.606    2.07e-08 ***
# Panel_typeRipple  -0.7891     0.1594    -4.949    7.44e-07 ***

# Post-hoc test
em <- emmeans(mod3, "Panel_type")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple   -0.137    0.166   Inf  -0.827   0.8415
# Flat - cm_25     0.926    0.165   Inf   5.606   <.0001
# Flat - cm_5      1.155    0.157   Inf   7.353   <.0001
# Ripple - cm_25   0.789    0.159   Inf   4.949   <.0001
# Ripple - cm_5    1.018    0.151   Inf   6.742   <.0001
# cm_25 - cm_5    -0.229    0.156   Inf  -1.469   0.4561

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot --------------------------------------------------------------------

# Averge taxon richness for configuration
srich_a_c_avg <- srich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
srich_a_c_avg
# Configuration  mean     se
# Grouped        1.57  0.123
# Random         1.57  0.123

# Create a plot using only the panel type, as the configuration was non-significant
# Average taxon richness for panel type
srich_a_p_avg <- srich %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
srich_a_p_avg
# Tile_type   mean     se
# cm_25       2.36  0.243 
# cm_5        3.03  0.205 
# Flat        0.94  0.095
# Ripple      1.08  0.078

# Re-order the tile types
srich_a_p_avg$Panel_type = factor(srich_a_p_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=srich_a_p_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sessile macroinvertebrate richness")+
  xlab("Configuration")+
  ylim(c(0, 4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate macrophyte taxon richness -----------------------------------------

# Remove rows where species weren't present
arich <- filter(macrophyte, Count_perc > 0)

# Calculate taxon richness by tile
arich <- arich %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
nrow(arich)
# Get the information for every tile
arich2 <- macrophyte %>% group_by(Frame_ID, Panel_ID, Panel_type) %>% 
  distinct(Species, .keep_all=TRUE) %>% 
  count(Configuration)
# Set the richness to 0
arich2$n <- 0
# Combine the two and remove duplicates
arich <- rbind(arich, arich2)
arich <- arich[order(arich$Panel_ID, -abs(arich$n) ), ]
arich <- arich[ !duplicated(arich$Panel_ID), ] 
nrow(arich) # 216 

rm(arich2)

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(arich$n, main="Density plot of macrophyte richness", xlab="Macrophyte taxon richness")
ggqqplot(arich$n, main="Q-Q plot of macrophyte taxon richness")

# Shapiro-Wilk test
shapiro.test(arich$n)
# p = 4.133e-15

# Try data transformation
# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(arich$n, plotit=FALSE))
# p = 9.915e-15
# The data can't be transformed

# Analysis ----------------------------------------------------------------

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=arich, family=poisson)
summary(mod1)
# AIC = 605.6        

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.9
summary(mod2)
# AIC = 599.9                

# Configuration is not significant, so simplify further by removing
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
anova(mod2, mod3) # p = 0.6
summary(mod3)
# AIC = 598.2            
check_model(mod3)
# plot(mod3) 
# qqnorm(resid(mod3))
# qqline(resid(mod3)) 
# hist(residuals(mod3))
# E2 <- resid(mod3)
# boxplot(E2 ~ Panel_type, data=arich) 
# abline(0,0)
# boxplot(E2 ~ Frame_ID, data=arich)
# abline(0,0)

summary(mod3)
# AIC      BIC     logLik    deviance   df.resid 
# 598.2    618.5   -293.1    586.2      210 

# Fixed effects:
#                  Estimate    Std. Error  z value  Pr(>|z|)    
# (Intercept)        0.8353     0.1098     7.610    2.74e-14 ***
# Panel_typecm_5    -0.1993     0.1636    -1.219     0.22297    
# Panel_typeFlat    -0.4392     0.1463    -3.002     0.00268 ** 
# Panel_typeRipple  -0.2600     0.1409    -1.845     0.06509 .  

# Post-hoc test
em <- emmeans(mod3, "Panel_type")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple    -0.1792  0.131   Inf  -1.368    0.5195
# Flat - cm_25      0.4392  0.146   Inf   3.002    0.0142 *
# Flat - cm_5       0.2398  0.155   Inf   1.546    0.4097
# Ripple - cm_25    0.2600  0.141   Inf   1.845    0.2524
# Ripple - cm_5     0.0606  0.150   Inf   0.404    0.9777
# cm_25 - cm_5      0.1993  0.164   Inf   1.219    0.6149

# Compact letter display
cld(em)
# Flat    a
# Ripple  ab
# 2.5 cm  b
# 5 cm    ab

# Plot --------------------------------------------------------------------

# Averge taxon richness for configuration
arich_a_c_avg <- arich %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
arich_a_c_avg
# Configuration  mean     se
# Grouped        1.74 0.0748
# Random         1.83 0.0727

# Create a plot using only the panel type, as the configuration was non-significant

# Average taxon richness for panel type
arich_a_p_avg <- arich %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
arich_a_p_avg
# Tile_type   mean     se
# cm_25       2.31  0.111 
# cm_5        1.89  0.153 
# Flat        1.49  0.081
# Ripple      1.78  0.077

# Re-order the tile types
arich_a_p_avg$Panel_type = factor(arich_a_p_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=arich_a_p_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average macrophyte richness")+
  xlab("Configuration")+
  ylim(c(0, 3))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate mobile taxon abundance ---------------------------------

mobile_a <- mobile %>% group_by(Frame_ID, Configuration, Panel_position, Panel_type) %>% 
  summarise(n=sum(Count_perc))
mobile_a

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(mobile_a$n, main="Density plot of mobile abundance", xlab="Mobile taxon abundance")
ggqqplot(mobile_a$n, main="Q-Q plot of mobile taxon abundance")

# Shapiro-Wilk test
shapiro.test(mobile_a$n)
# p = < 2.2e-16

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(mobile_a$n, plotit=FALSE))
# p = 7.204e-13

# Analysis ----------------------------------------------------------------

# Interaction between configuration and panel type
set.seed(1234)
mod1 <- glmer(n ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=mobile_a, family=poisson)
summary(mod1)
# AIC = 830.7        

# The interaction term is not significant, so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.12
summary(mod2)
# AIC = 830.6   

# Configuration is not significant
set.seed(1234)
mod3 <- update(mod2, .~. -Configuration)
# Compare models with a likelihood ratio test
anova(mod2, mod3) # p = 0.1
summary(mod3)
# AIC = 831.1        
check_model(mod3) 
# plot(mod3) 
# qqnorm(resid(mod3))
# qqline(resid(mod3)) 
# hist(residuals(mod3))
# E2 <- resid(mod3)
# boxplot(E2 ~ Panel_type, data=mobile_a) 
# abline(0,0)
# boxplot(E2 ~ Configuration, data=mobile_a) 
# abline(0,0)
# boxplot(E2 ~ Frame_ID, data=mobile_a)
# abline(0,0)

summary(mod3)
# AIC      BIC     logLik    deviance   df.resid 
# 831.1    851.3   -409.5    819.1      210 

# Fixed effects:
#                  Estimate   Std. Error  z value  Pr(>|z|)    
# (Intercept)        1.19007    0.16941   7.025    2.14e-12 ***
# Panel_typecm_5     0.07443    0.18893   0.394       0.694    
# Panel_typeFlat    -1.43732    0.15632  -9.195     < 2e-16 ***
# Panel_typeRipple  -1.21725    0.14647  -8.310     < 2e-16 ***

# Post-hoc test
em <- emmeans(mod3, "Panel_type")
contrast(em, "pairwise", adjust="Tukey")
# contrast        estimate  SE      df    z.ratio  p.value
# Flat - Ripple    -0.2201  0.165   Inf  -1.334    0.5409
# Flat - cm_25      1.4373  0.156   Inf   9.195    <.0001
# Flat - cm_5       1.5117  0.177   Inf   8.560    <.0001
# Ripple - cm_25    1.2173  0.146   Inf   8.310    <.0001
# Ripple - cm_5     1.2917  0.168   Inf   7.690    <.0001
# cm_25 - cm_5     -0.0744  0.189   Inf  -0.394    0.9793

# Compact letter display
cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot ------------------------------------------------------------------------

mobile_a_avg <- mobile_a %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
mobile_a_avg
# Panel_type  mean    se
# cm_25       4.58  0.602
# cm_5        3.25  0.550
# Flat        0.90  0.220
# Ripple      1.12  0.194

# Re-order the tile types
mobile_a_avg$Panel_type = factor(mobile_a_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=mobile_a_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average mobile taxon abundance")+
  xlab("Panel type")+
  ylim(c(0, 6))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

mobile_a_avg_m <- mobile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
mobile_a_avg_m
# Configuration   mean     se
# Grouped         1.59  0.232
# Random          2.37  0.302

ggplot(data=mobile_a_avg_m, aes(x=Configuration, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average mobile taxon abundance")+
  xlab("Configuration")+
  ylim(c(0, 4))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Calculate sessile taxon abundance ---------------------------------

sessile <- filter(sessile, Species=="Cirripedia"|Species=="C_eumyota"|Species=="O_littoralis"|Species=="S_spirorbis"|Species=="Spirobranchus"|Species=="W_subatra")

sessile_a <- sessile %>% group_by(Frame_ID, Configuration, Panel_position, Panel_type) %>% 
  summarise(n=sum(Count_perc))
sessile_a

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(sessile_a$n, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggqqplot(sessile_a$n, main="Q-Q plot of sessile taxon abundance")

# Shapiro-Wilk test
shapiro.test(sessile_a$n)
# p = 2.548e-16

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(sessile_a$n, plotit=FALSE))
# p = 2.36e-05

# Analysis ----------------------------------------------------------------

# Use a beta family because the values are proportions/percentages
# Normalise the data between 0 and 1 (but actually 0.01 and 0.99 because glmmTMB
# with a beta family won't take absolute 0 and 1 values)
sessile_a$n2 <- rescale(sessile_a$n, to=c(0.01, 0.99))

ggdensity(sessile_a$n2, main="Density plot of sessile abundance", xlab="Sessile taxon abundance")
ggqqplot(sessile_a$n2, main="Q-Q plot of sessile taxon abundance")

set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=sessile_a, family=beta_family())
summary(mod1)
# AIC = -487.9
check_model(mod1)
# qqnorm(resid(mod1))
# qqline(resid(mod1)) 
# hist(residuals(mod1))
# E2 <- resid(mod1)
# boxplot(E2 ~ Panel_type, data=sessile_a) 
# abline(0,0)
# boxplot(E2 ~ Configuration, data=sessile_a) 
# abline(0,0)

summary(mod1)
#    AIC      BIC   logLik  deviance  df.resid 
# -487.9   -450.8    255.0    -509.9       205

# Conditional model:
#                                      Estimate  Std. Error  z value   Pr(>|z|)    
# (Intercept)                          -1.26032     0.19124   -6.590   4.39e-11 ***
# ConfigurationRandom                   0.03204     0.26603    0.120   0.90414    
# Panel_typecm_5                        0.31176     0.26746    1.166   0.24375    
# Panel_typeFlat                       -1.41575     0.23636   -5.990   2.10e-09 ***
# Panel_typeRipple                     -0.81124     0.22214   -3.652   0.00026 ***
# ConfigurationRandom:Panel_typecm_5    0.76848     0.35867    2.143   0.03215 *  
# ConfigurationRandom:Panel_typeFlat    0.03989     0.32303    0.123   0.90173    
# ConfigurationRandom:Panel_typeRipple -0.05220     0.30919   -0.169   0.86593  

# Post-hoc test
em <- emmeans(mod1, ~Configuration:Panel_type)
contrast(em, "pairwise", adjust="Tukey")

cld(em)
# Grouped Flat    a
# Grouped Ripple  a
# Grouped 2.5 cm  b
# Grouped 5 cm    b
# Random  Flat    a
# Random  Ripple  a
# Random  2.5 cm  b
# Random  5 cm    c

# Plot ------------------------------------------------------------------------

# Configuration
sessile_a_c_avg <- sessile_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
sessile_a_c_avg
# Configuration  mean     se
# Grouped        8.33  0.867
# Random         9.92   1.17

# Interaction
sessile_a_avg <- sessile_a %>% group_by(Configuration, Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
sessile_a_avg
sessile_a_avg$Panel_type = factor(sessile_a_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=sessile_a_avg, aes(x=Configuration, y=mean, fill=Panel_type))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average sessile taxon cover (%)")+
  xlab("Configuration")+
  scale_y_continuous(limits=c(0, 35), breaks=seq(0, 35, 5))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Calculate macrophyte abundance ---------------------------------

macro_a <- macrophyte %>% group_by(Frame_ID, Configuration, Panel_position, Panel_type) %>% 
  summarise(n=sum(Count_perc))
macro_a

# Check for normality -----------------------------------------------------

# Check for normality
ggdensity(macro_a$n, main="Density plot of macrophyte abundance", xlab="macrophyte abundance")
ggqqplot(macro_a$n, main="Q-Q plot of macrophyte abundance")

# Shapiro-Wilk test
shapiro.test(macro_a$n)
# p = 3.829e-06

# Tukey's Ladder of Powers transformation
ggdensity(transformTukey(macro_a$n, plotit=FALSE))
# p = 1.094e-05

# Analysis ----------------------------------------------------------------

macro_a$n2 <- rescale(macro_a$n, to=c(0.01, 0.99))

ggdensity(macro_a$n2, main="Density plot of macrophyte abundance", xlab="macrophyte abundance")
ggqqplot(macro_a$n2, main="Q-Q plot of macrophyte abundance")

set.seed(1234)
mod1 <- glmmTMB(n2 ~ Configuration*Panel_type + (1|Configuration/Frame_ID), data=macro_a, family=beta_family())
summary(mod1)
# AIC = -129.4                       
# The interaction term is not significant (P = 0.09), so try removing 
set.seed(1234)
mod2 <- update(mod1, .~. -Configuration:Panel_type)
# Compare models with a likelihood ratio test
anova(mod1, mod2) # p = 0.9
summary(mod2)
# AIC = -134.7     
# Configuration is marginally significant
check_model(mod2) 
# qqnorm(resid(mod2))
# qqline(resid(mod2)) 
# hist(residuals(mod2))
# E2 <- resid(mod2)
# boxplot(E2 ~ Panel_type, data=macro_a) 
# abline(0,0)
# boxplot(E2 ~ Configuration, data=macro_a) 
# abline(0,0)

summary(mod2)
#    AIC      BIC   logLik  deviance  df.resid 
# -134.7   -107.7     75.4    -150.7       208 

# Conditional model:
#                       Estimate  Std. Error  z value   Pr(>|z|)    
# (Intercept)           0.04346   0.16474     0.264     0.79191    
# ConfigurationRandom  -0.31024   0.15822    -1.961     0.04990 *  
# Panel_typecm_5        0.06451   0.19857     0.325     0.74527    
# Panel_typeFlat       -0.74932   0.15948    -4.698     2.62e-06 ***
# Panel_typeRipple     -0.47946   0.15815    -3.032     0.00243 **  

# Post-hoc test
em <- emmeans(mod2, ~Panel_type)
contrast(em, "pairwise", adjust="Tukey")
# Flat - Ripple   -0.2699 0.126 Inf  -2.136  0.1416
# cm_25 - Flat     0.7493 0.159 Inf   4.698  <.0001
# cm_5 - Flat      0.8138 0.160 Inf   5.079  <.0001
# cm_25 - Ripple   0.4795 0.158 Inf   3.032  0.0130
# cm_5 - Ripple    0.5440 0.159 Inf   3.415  0.0036
# cm_25 - cm_5    -0.0645 0.199 Inf  -0.325  0.9882

cld(em)
# Flat    a
# Ripple  a
# 2.5 cm  b
# 5 cm    b

# Plot ------------------------------------------------------------------------

macro_a_avg <- macro_a %>% group_by(Panel_type) %>% 
  summarise(mean=mean(n), se=std.error(n))
macro_a_avg
# Panel_type  mean    se
# cm_25       89.1  7.18
# cm_5        86.8  4.94
# Flat        55.9  3.83
# Ripple      65.6  3.80

# Re-order the tile types
macro_a_avg$Panel_type = factor(macro_a_avg$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=macro_a_avg, aes(x=Panel_type, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average macrophyte abundance")+
  xlab("Panel type")+
  ylim(c(0, 100))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

macro_a_avg_m <- macro_a %>% group_by(Configuration) %>% 
  summarise(mean=mean(n), se=std.error(n))
macro_a_avg_m
# Mosaic   mean     se
# Grouped  75.6  3.35
# Random   64.0  3.56

ggplot(data=macro_a_avg_m, aes(x=Configuration, y=mean))+
  geom_bar(width=.7, stat="identity", position=position_dodge(.8))+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.2, end=0.8)+
  ylab("Average macrophyte abundance")+
  xlab("Configuration")+
  ylim(c(0, 100))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.8))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="none")

# Plot the average count of mobile species groups -------------------------

# Configuration
# Calculate average abundance by tile
mobileavgt <- mobile %>% group_by(Configuration, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
mobileavgt <- filter(mobileavgt, !Func_group=="Unknown")
ggplot(data=mobileavgt, aes(x=Configuration, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (count)")+
  xlab("Configuration")+
  ylim(c(0, 1))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Panel type
# Calculate average abundance by panel
mobileavgt <- mobile %>% group_by(Panel_type, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
mobileavgt <- filter(mobileavgt, !Func_group=="Unknown")
# Re-order the tile types
mobileavgt$Panel_type = factor(mobileavgt$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=mobileavgt, aes(x=Panel_type, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (count)")+
  xlab("Panel type")+
  ylim(c(0, 1.5))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Plot the average cover of sessile species groups -------------------------

# Don't bother, they're all suspension feeders!

# Plot the average count of macrophyte groups -------------------------

# Configuration
# Calculate average abundance by tile
macroavgt <- macrophyte %>% group_by(Configuration, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
macroavgt <- filter(macroavgt, !Func_group=="Unknown")
ggplot(data=macroavgt, aes(x=Configuration, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (%)")+
  xlab("Configuration")+
  ylim(c(0, 40))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))

# Panel type
# Calculate average abundance by panel
macroavgt <- macrophyte %>% group_by(Panel_type, Func_group) %>% 
  summarise(mean=mean(Count_perc), se=std.error(Count_perc))
macroavgt <- filter(macroavgt, !Func_group=="Unknown")
# Re-order the tile types
macroavgt$Panel_type = factor(macroavgt$Panel_type, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)
ggplot(data=macroavgt, aes(x=Panel_type, y=mean, fill=Func_group))+
  geom_bar(width=.7, stat="identity", position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma", begin=0.1, end=0.9)+
  ylab("Average abundance (%)")+
  xlab("Panel type")+
  ylim(c(0, 40))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position=c(0.1, 0.9))
