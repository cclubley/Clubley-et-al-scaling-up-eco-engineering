# Load required packages --------------------------------------------------

{library(reshape2)
  library(dplyr)
  library(vegan)
  library(ggplot2)
  library(ccplot)
  library(pracma)
}

# Load the data -----------------------------------------------------------

bio <- read.csv("./Biodiversity_data.csv")
head(bio)

# Change data classifications
str(bio)
bio$Frame_ID <- as.factor(bio$Frame_ID)
bio$Configuration <- as.factor(bio$Configuration)
bio$Panel_ID <- as.factor(bio$Panel_ID)
bio$Panel <- as.factor(bio$Panel)
bio$Type <- as.factor(bio$Type)
str(bio)

# Remove the unneeded data
bio <- filter(bio, !Species=="Sediment")

# Format the data ---------------------------------------------------------

# Convert to wide format
bio2 <- dcast(bio, Frame_ID+Configuration+Panel_ID+Panel ~ Species, value.var = c("Count_perc"))
colnames(bio2)

## Write out the .csv file and do the final formatting for PRIMER in Excel
#write.csv(bio2, "bio_for_primer.csv", row.names=FALSE)

# nMDS plot ---------------------------------------------------------------

# Convert species data to a matrix
bio_mat <- as.matrix(bio2[, 5:33])

# Transform/standardise the data
# Use a fourth-root transformation to minimise the influence of the most 
# abundant groups.
bio_mat <- nthroot(bio_mat, 4)
bio_mat[is.na(bio_mat)] <- 0 # Deal with any NA entries

# Configure samples in 2-dimensional space
bioMDS <- metaMDS(bio_mat, distance="bray", maxit=999) 
# Outputs the iteration of the nMDS until a solution is reached 
bioMDS
# Stress value = 0.16 - matches stress value from PRIMER

# Export the data from the nmds object
NMDS1 <- bioMDS$points[,1] 
NMDS2 <- bioMDS$points[,2]
bio_plot <- cbind(bio2, NMDS1, NMDS2)

# Re-order panel types
bio_plot$Panel = factor(bio_plot$Panel, levels = c("Flat", "Ripple", "cm_25", "cm_5"), ordered=TRUE)

# Plot ordination (Arrangement*Tile_type)
ggplot(bio_plot, aes(NMDS1, NMDS2, color=Configuration, shape=Panel))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Configuration:Panel), position=position_jitter(.1))+ # separates overlapping points
  theme_minimal()

# Fit vectors to ordination - which variables (species) are correlated with the ordination?
fit <- envfit(bioMDS, bio_mat)
arrow <- data.frame(fit$vectors$arrows, R=fit$vectors$r, P=fit$vectors$pvals)
arrow$sp <- rownames(arrow)
arrow_p <- filter(arrow, sp=="F_vesiculosus"| sp=="W_subatra"| sp=="F_spiralis"| sp=="U_intestinalis")
# Species chosen based on output in PRIMER - these are the species for which Pearson's 
# correlation was > 0.05

ggplot(bio_plot, aes(NMDS1, NMDS2))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Panel, shape=Configuration), size=2, position=position_jitter(.1))+ # separates overlapping points
  geom_segment(arrow_p, mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length=unit(.2, "cm")*arrow_p$R))+ # add arrows (scaled by R-squared value)
  geom_text(arrow_p, mapping=aes(x=NMDS1, y=NMDS2, label=sp))+
  xlim(c(-1.4, 1.2))+
  ylim(c(-1.2, 1.2))+
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="right")

# Just configuration
ggplot(bio_plot, aes(NMDS1, NMDS2))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Configuration), position=position_jitter(.1))+ # separates overlapping points
  stat_ellipse(aes(fill=Configuration), alpha=.2, type='t', size =1, geom="polygon")+ # changes shading on ellipses
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="right")

# Just panel
ggplot(bio_plot, aes(NMDS1, NMDS2))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Panel), position=position_jitter(.1))+ # separates overlapping points
  stat_ellipse(aes(fill=Panel), alpha=.2, type='t', size =1, geom="polygon")+ # changes shading on ellipses
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="right")
