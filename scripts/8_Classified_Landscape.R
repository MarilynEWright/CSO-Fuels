# Script name: 8_Classified_Landscape.R
# Author: M. E. Wright
# Date created: April 26, 2023
#

# Notes: 
#      
# ---

# ---
# 1. Loading packages and data ----
# ---
library(tidyverse)
library(sf)
library(sp)
library(terra)
library(rgdal)
library(ggplot2)


# load in the probability raster 
prob <- rast("data/probability_rasters/composite_probability.tif")

# load in the ACCEL proportion of maxSDI raster
SDI <- rast("data/SDI/ACCEL_SDI/proportion_of_SDI_83_Max_30m.tif")

# reproject the SDI layer to match WGS84 UTM10
SDI <- terra::project(SDI, crs(prob))

# ---
# 2. Create composite SDI/probability raster ----
# ---

# reclassify the SDI layer 
# <= 0.35 --> 1
# >  0.35 --> 2

SDI[SDI > 0.35] <- 2
SDI[SDI <= 0.35] <- 1


# because SDI is 30-m pixels and prob is 10-m, we need to resample the probability
# raster to match SDI
prob.resample <- resample(prob, SDI)


# reclassify resampled probability layer 
# <= 0.0326 --> 4
# >  0.0326 --> 5

prob.resample[prob.resample > 0.0326] <- 5
prob.resample[prob.resample <= 0.0326] <- 4


# divide the SDI layer by the probability layer to get 4 unique categories
# 0.25 = no priority actions (low proportion max SDI and low probability of use)
# 0.50 = fuels treatment (high proportion max SDI and low probability of use)
# 0.40 = conflict zones (high proportion max SDI and high probability of use)
# 0.20 = conservation (low proportion max SDI and high probability of use)

composite <- SDI / prob.resample

## SAVE COMPOSITE RASTER ##
writeRaster(composite, "data/results/treatment_categories.tif")


# ---
# 3. Consider proportions of landscape in each category ----
# ---

# load in composite raster 
comp <- rast("data/results/treatment_categories.tif")


# load in the CA land ownership shapefile
owner <- st_read("data/CA_land_ownership/ownership22_1.shp")

# load in CSO range shapefile
CSO_range <- st_read("data/CSO_range/Spotted_Owl_Range_-_CWHR_B270_[ds897].shp")
  
# load in Sierra Nevada shapefile (for clipping ownership layer)
SNC <- st_read("data/shapefiles/SNC_Boundary.shp")

# reproject the ownership layer to match WGS84 UTM10
owner <- st_transform(owner, crs = crs(comp))
SNC <- st_transform(SNC, crs = crs(comp))
CSO_range <- st_transform(CSO_range, crs = crs(comp))

# ---
# Amend the full Sierra Nevada extent to include only CSO habitat area (west slope)
# clip CSO_range to Sierra Nevada extent
CSO.Sierra <- st_intersection(CSO_range, SNC)

# create spatial vector of CSO.Sierra 
CSO.vect <- vect(CSO.Sierra)

# mask composite by CSO habitat
CSO.comp <- mask(comp, mask = CSO.vect)
# ---


# ---
# SIERRA NEVADA FULL LANDSCAPE #
# number of cells in each category 
Sierra <- freq(comp, digits = 2)

# add a classification column 
Sierra <- Sierra %>%
  mutate(Classification = ifelse(value == 0.50, "Fuel Treatment", 
                                 ifelse(value == 0.40, "Conflict Zone",
                                        ifelse(value == 0.20, "Conservation", "Low Priority Area"))))

# add a proportion column 
total.Sierra <- sum(Sierra$count) # total count of all cells

Sierra <- Sierra %>%
  mutate(proportion = count/total.Sierra)

# add proportion of Sierra column 
Sierra <- Sierra %>%
  mutate(proportion.total = count/total.Sierra)



# ---
# CSO HABITAT WITHIN SIERRA NEVADA FULL LANDSCAPE #
# number of cells in each category 
CSO.Sierra <- freq(CSO.comp, digits = 2)

# add a classification column 
CSO.Sierra <- CSO.Sierra %>%
  mutate(Classification = ifelse(value == 0.50, "Fuel Treatment", 
                                 ifelse(value == 0.40, "Conflict Zone",
                                        ifelse(value == 0.20, "Conservation", "Low Priority Area"))))

# add a proportion column 
total.CSO.Sierra <- sum(CSO.Sierra$count) # total count of all cells

CSO.Sierra <- CSO.Sierra %>%
  mutate(proportion = count/total.CSO.Sierra)

# add proportion of Sierra column 
CSO.Sierra <- CSO.Sierra %>%
  mutate(proportion.total = count/total.CSO.Sierra)




# ---

# clip the ownership layer to the CSO habitat extent 
owner.CSO <- st_intersection(owner, CSO_SNC)

# separate the ownership layer into distinct groups
USFS <- owner.CSO %>%
 filter(OWN_GROUP == "USDA Forest Service")

NP <- owner.CSO %>%
  filter(OWN_GROUP == "National Park Service")

# for private, do a reverse clip (difference) with the CSO habitat boundary
private <- st_difference(CSO_SNC, st_union(st_geometry(owner.CSO)))

# ---
# create SpatVect for each layer 
USFS.vect <- vect(USFS)
NP.vect <- vect(NP)
Private.vect <- vect(private)

# ---
# mask composite by each ownership layer
USFS.comp <- mask(CSO.comp, mask = USFS.vect)
NP.comp <- mask(CSO.comp, mask = NP.vect)
Private.comp <- mask(CSO.comp, mask = Private.vect)

# --- 
# proportion of cells in each category 

# US FOREST SERVICE LAND #
USFS.cat <- freq(USFS.comp, digits = 2)

# add a classification column 
USFS.cat <- USFS.cat %>%
  mutate(Classification = ifelse(value == 0.50, "Fuel Treatment", 
                                 ifelse(value == 0.40, "Conflict Zone",
                                        ifelse(value == 0.20, "Conservation", "Low Priority Area"))))


# add a proportion column 
total <- sum(USFS.cat$count) # total count of all cells

USFS.cat <- USFS.cat %>%
  mutate(proportion = count/total)

# add proportion of Sierra column 
USFS.cat <- USFS.cat %>%
  mutate(proportion.total = count/total.CSO.Sierra)

# ---
# NATIONAL PARK LAND #
NP.cat <- freq(NP.comp, digits = 2)

# add a classification column 
NP.cat <- NP.cat %>%
  mutate(Classification = ifelse(value == 0.50, "Fuel Treatment", 
                                 ifelse(value == 0.40, "Conflict Zone",
                                        ifelse(value == 0.20, "Conservation", "Low Priority Area"))))


# add a proportion column 
total <- sum(NP.cat$count) # total count of all cells

NP.cat <- NP.cat %>%
  mutate(proportion = count/total)

# add proportion of Sierra column 
NP.cat <- NP.cat %>%
  mutate(proportion.total = count/total.CSO.Sierra)

# ---
# PRIVATE LAND #
Private.cat <- freq(Private.comp, digits = 2)

# add a classification column 
Private.cat <- Private.cat %>%
  mutate(Classification = ifelse(value == 0.50, "Fuel Treatment", 
                                 ifelse(value == 0.40, "Conflict Zone",
                                        ifelse(value == 0.20, "Conservation", "Low Priority Area"))))


# add a proportion column 
total <- sum(Private.cat$count) # total count of all cells

Private.cat <- Private.cat %>%
  mutate(proportion = count/total)

# add proportion of Sierra column 
Private.cat <- Private.cat %>%
  mutate(proportion.total = count/total.CSO.Sierra)

# ---
# create one dataframe with all proportions #
# add ownership variable to each dataframe 
Sierra <- Sierra %>%
  mutate(level = "Sierra bioregion")
CSO.Sierra <- CSO.Sierra %>%
  mutate(level = "CSO Habitat w/in Sierra bioregion")
USFS.cat <- USFS.cat %>%
  mutate(level = "National Forest")
NP.cat <- NP.cat %>%
  mutate(level = "National Park")
Private.cat <- Private.cat %>%
  mutate(level = "Private land")

# create one dataframe
proportions <- rbind(Sierra, CSO.Sierra, USFS.cat, NP.cat, Private.cat)

# try ggplot
# filter out Sierra bioregion for this **
plot <- proportions %>%
  filter(level != "Sierra bioregion" & level != "CSO Habitat w/in Sierra bioregion")

plot$level <- factor(plot$level, levels = c("National Park",
                                            "National Forest",
                                            "Private land"))
plot$Classification <- factor(plot$Classification, levels = c("Low Priority Area",
                                                              "Fuel Treatment",
                                                              "Conflict Zone",
                                                              "Conservation"))
## SAVE plot FOR FIGURE CREATION ##
saveRDS(plot, "figures/figure_data/CSO_habitat_proportion_by_ownership")

# FIGURE 5: CLASSIFICATION BY OWNERSHIP TYPE ----
library(patchwork)
library(grid)


A <- ggplot(plot, aes(x= Classification, y= proportion.total, fill= level)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("coral3", "darkolivegreen4", "blue4")) +
  ylab("Proportion of Sierra Bioregion Landscape") +
  xlab("") +
  ylim(ymin= 0, ymax = 0.5) +
  theme_classic() +
  guides(fill= guide_legend(title = "Ownership")) +
  theme(plot.title = element_text(hjust = 0.5)) 

# load images for quadrants
stop <- readPNG("figures/images/stop.png")
caution <- readPNG("figures/images/caution.png")
axe <- readPNG("figures/images/axe.png")
tree <- readPNG("figures/images/tree.png")

# add to graph 
A <-  A +
  annotation_custom(rasterGrob(stop, interpolate = TRUE), xmin = 0.5, xmax = 1.5,
                    ymin = 0.4, ymax = 0.5) +
  annotation_custom(rasterGrob(axe, interpolate = TRUE), xmin = 1.5, xmax = 2.5,
                    ymin = 0.3, ymax = 0.4) +
  annotation_custom(rasterGrob(caution, interpolate = TRUE), xmin = 2.5, xmax = 3.5,
                    ymin = 0.2, ymax = 0.3) +
  annotation_custom(rasterGrob(tree, interpolate = TRUE), xmin = 3.5, xmax = 4.5,
                    ymin = 0.05, ymax = 0.15)

tiff("figures/Classification_by_Owner.tiff", units="in", width=8, height=5, res=300)
A
dev.off()
