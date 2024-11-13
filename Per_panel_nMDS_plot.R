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

# Separate the species
spc <- bio_plot[, 5:35]

# Create interaction term
int <- interaction(bio_plot$Configuration, bio_plot$Panel)

# Run SIMPER for interaction
si <- simper(spc, int)
summary(si)

# Create an empty data frame to store the results
summary_table <- data.frame(Contrast = character(),
                            Species = character(),
                            Contribution = numeric(),
                            stringsAsFactors = FALSE)

# Loop over each contrast to extract and store the relevant data
for (i in seq_along(si)) {
  contrast_name <- names(si)[i]
  species_contrib <- as.data.frame(si[[i]]$average)
  colnames(species_contrib) <- "avg"
  species_contrib$sp <- row.names(species_contrib)
  significant_species <- filter(species_contrib, avg >= 0.02)
  
  if (nrow(significant_species) == 0) next
  
  # Add the data to the summary table
  for (j in 1:nrow(significant_species)) {
    summary_table <- rbind(summary_table, 
                           data.frame(Contrast = contrast_name,
                                      Species = rownames(significant_species)[j],
                                      Contribution = significant_species[j, "avg"],
                                      stringsAsFactors = FALSE))
  }
}

table(summary_table$Species)

# Extract species scores from metaMDS result
species_scores <- as.data.frame(scores(bioMDS, "species"))
# Filter species scores to include only significant species from SIMPER
significant_species_scores <- species_scores[rownames(species_scores) %in% summary_table$Species, ]

# Plot
ggplot(bio_plot, aes(NMDS1, NMDS2))+
  geom_point(bio_plot, mapping=aes(NMDS1, NMDS2, color=Panel, shape=Configuration), size=2, position=position_jitter(.1))+ 
  geom_segment(data=significant_species_scores, aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length=unit(0.2, "cm")), 
               color="red", size=0.5)+  # Add vectors for significant species
  geom_text(data=significant_species_scores, aes(x=NMDS1, y=NMDS2, label=rownames(significant_species_scores)),
            color="darkred", size=3, vjust=-0.5) +  # Label species at the end of their vectors
  cc_theme()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.position="right")

