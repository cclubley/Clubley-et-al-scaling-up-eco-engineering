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
bio$Panel_type <- as.factor(bio$Panel_type)
bio$Type <- as.factor(bio$Type)
bio$Measure <- as.factor(bio$Measure)
str(bio)

# Remove the unneeded data
bio <- filter(bio, !Species=="Sediment")

# Calculate the biodiversity at the frame level
bio <- bio %>% group_by(Frame_ID, Configuration, Species, Func_group, Type, Measure) %>% 
  summarise(Count_perc=sum(Count_perc))

# Correct for the percentages
for(i in 1:nrow(bio)){
  if(bio$Measure[i] == "Percent"){
    bio$Count_perc[i] <- bio$Count_perc[i]*(9/100)
  }
}

# Format the data ---------------------------------------------------------

# Convert to wide format
bio2 <- dcast(bio, Frame_ID+Configuration ~ Species, value.var=c("Count_perc"))
colnames(bio2)

# Write out the .csv file and do the final formatting for PRIMER in Excel
#write.csv(bio2, "bio_for_primer_frame.csv", row.names=FALSE)

# nMDS plot ---------------------------------------------------------------

# Convert species data to a matrix
bio_mat <- as.matrix(bio2[, 3:31])

# Transform/standardise the data
# Use a fourth-root transformation to minimise the influence of the most 
# abundant groups.
bio_mat <- nthroot(bio_mat, 4)
bio_mat[is.na(bio_mat)] <- 0 # Deal with any NA entries

# Configure samples in 2-dimensional space
bioMDS <- metaMDS(bio_mat, distance="bray", maxit=999) 
# Outputs the iteration of the nMDS until a solution is reached 
bioMDS
# Stress value = 0.19 - matches stress value from PRIMER

# Export the data from the nmds object
NMDS1 <- bioMDS$points[,1] 
NMDS2 <- bioMDS$points[,2]
bio_plot <- cbind(bio2, NMDS1, NMDS2)

# Just mosaic
ggplot(bio_plot, aes(NMDS1, NMDS2))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Configuration), position=position_jitter(.1))+ # separates overlapping points
  stat_ellipse(aes(fill=Configuration), alpha=.2, type='t', size =1, geom="polygon")+ # changes shading on ellipses
  cc_theme()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.position="right")

