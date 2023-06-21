# Script name: 7_SDI.R
# Author: M. E. Wright
# Date created: April 21, 2023
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
library(scales)
library(rgdal)
library(ggplot2)
library(png)
library(grid)
library(RCurl)

# load in the probability raster 
prob <- rast("data/probability_rasters/composite_probability.tif")

# load in the ACCEL proportion of maxSDI raster
SDI <- rast("data/SDI/ACCEL_SDI/proportion_of_SDI_83_Max_30m.tif")

# load in the CA land ownership shapefile
owner <- st_read("data/CA_land_ownership/ownership22_1.shp")

# check the crs of each layer (project to WGS84 UTM10 as needed)

# reproject the SDI layer to match WGS84 UTM10
SDI <- terra::project(SDI, crs(prob))
owner <- st_transform(owner, crs = crs(prob))


# ---
# 2. Point cloud ----
# ---

# take a regular sample across the extent 
point_sample <- spatSample(prob, size = 50000, method = "regular", replace = FALSE,
                           na.rm = TRUE, as.points = TRUE)

# pull out all data from sample (probability of use, SDI, and ownership)
point_SDI <- terra::extract(SDI, point_sample)
point_prob <- terra::extract(prob, point_sample)

# for ownership data, need to create an sf
point_sf <- st_as_sf(point_sample) # sf object for point layer

#point_sf <- st_transform(point_sf, crs = crs(owner)) # transform and match coordinates

point_sf <- point_sf %>%
  mutate(ID = row_number()) # add an id column 

point_owner <- st_intersection(owner, point_sf) # intersect point layer and ownership


# create a dataframe with all information 
point_sample_info <- data.frame(point_prob, point_SDI$proportion_of_SDI_83_Max_30m)

point_sample_info <- left_join(point_sample_info, point_owner, by = "ID")


# rename variables in dataframe
point_sample_info <- point_sample_info %>%
  rename(Prob_of_use = CFO.California.CanopyBaseHeight.Summer2020.00010m.x,
         Prop_max_SDI = point_SDI.proportion_of_SDI_83_Max_30m)

# select only columns of interest and rename ownership
point_select <- point_sample_info %>%
  dplyr::select(c(Prob_of_use, Prop_max_SDI, OWN_GROUP)) %>%
  rename(ownership = OWN_GROUP)

# drop records with no proportion of max SDI information 
point_select <- point_select %>%
  drop_na(Prop_max_SDI)

# rename "NA" in ownership to "private"
point_select[is.na(point_select)] <- "Private"

# rename "ownership" to 
point_select <- point_select %>%
  rename(agency = ownership) 

# convert agency to factor 
point_select$agency <- as.factor(point_select$agency)

# collapse into more meaningful categories 
point_select$ownership <- fct_collapse(point_select$agency,
                          Forest_Service = "USDA Forest Service",
                          National_Park = "National Park Service",
                          Other_federal = c("Bureau of Indian Affairs",
                                            "Bureau of Land Management",
                                            "Bureau of Reclamation",
                                            "Department of Defense", 
                                            "US Fish and Wildlife Service"),
                          State = c("CA Dept. of Fish and Wildlife",
                                    "CA Dept. of Forestry and Fire Protection",
                                    "CA Dept. of Parks and Recreation",
                                    "Local Government",
                                    "Other State Lands"),
                          Private = "Private")

## SAVE RDS ##
saveRDS(point_select, "data/results/point_sample")



# ---
# 3. SDI-Use-Ownership plots ----
# ---

# open as needed 
point_select <- readRDS("data/results/point_sample")

# ---
# # check on scaled threshold (Prob of use threshold = 0.399396911 for SCALED DATA)
# df2 <- data.frame(Prob_of_use = 0.03260122, Prop_max_SDI = 4.756024e-02, 
#                   agency = "NONE", ownership = "NONE")
# test <- rbind(point_select, df2)
# 
# # rescale the probability data 
# test$Prob_of_use_scaled <- rescale(test$Prob_of_use)
# ---

# rescale the probability data 
point_select$Prob_of_use_scaled <- rescale(point_select$Prob_of_use)


ggplot(point_select, aes(Prob_of_use_scaled, Prop_max_SDI, col = ownership)) +
  geom_point()

# look at different ownerships 
NP <- point_select %>%
  filter(ownership == "National_Park")

ggplot(NP, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "brown")


USFS <- point_select %>%
  filter(ownership == "Forest_Service")

ggplot(USFS, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "darkgreen")


Private <- point_select %>%
  filter(ownership == "Private")

ggplot(Private, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "blue")


# three selected categories 
filter <- point_select %>%
  filter(ownership %in% c("Forest_Service", "National_Park", "Private"))

ggplot(filter, aes(Prob_of_use_scaled, Prop_max_SDI, color = ownership)) +
  geom_point() +
  geom_vline(xintercept = 0.40, linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = 0.35, linetype = "dashed", color = "black", size = 1)


##########################
# attempt at averaging 
test <- filter %>%
  group_by(ownership) %>%
  summarize(mean_SDI = mean(Prop_max_SDI),
            mean_prob = mean(Prob_of_use_scaled),
            sd_SDI = sd(Prop_max_SDI),
            sd_prob = sd(Prob_of_use_scaled))

test$ownership <- factor(test$ownership,
                         levels = c("National_Park", "Forest_Service", "Private"))

ggplot(test, aes(mean_prob, mean_SDI, color = ownership)) +
  geom_point() +
  geom_errorbar(aes(ymin= mean_SDI-sd_SDI, ymax= mean_SDI+sd_SDI), width = 0) +
  geom_errorbar(aes(xmin = mean_prob-sd_prob, xmax= mean_prob+sd_prob),  width = 0) +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.40, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 01) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  theme_classic()




# USFS point density plot 
ggplot(USFS, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "darkolivegreen4") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "darkgreen") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.40, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  theme_classic()

# National Park point density plot 
ggplot(NP, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "coral3") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "orangered4") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.4, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  theme_classic()

# Private point density plot 
ggplot(Private, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "blue4") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "darkslateblue") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.4, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  theme_classic()



# ---
# 4. FIGURE CREATION ----
# ---
library(patchwork)
# open as needed 
point_select <- readRDS("data/results/point_sample")


# rescale the probability data 
point_select$Prob_of_use_scaled <- rescale(point_select$Prob_of_use)

# three selected categories 
filter <- point_select %>%
  filter(ownership %in% c("Forest_Service", "National_Park", "Private"))

# FIGURE 4. RESILIENCE-USE BY OWNERSHIP ----

# average the point data 
test <- filter %>%
  group_by(ownership) %>%
  summarize(mean_SDI = mean(Prop_max_SDI),
            mean_prob = mean(Prob_of_use_scaled),
            sd_SDI = sd(Prop_max_SDI),
            sd_prob = sd(Prob_of_use_scaled))

test$ownership <- factor(test$ownership,
                         levels = c("National_Park", "Forest_Service", "Private"))
# look at different ownerships 
NP <- point_select %>%
  filter(ownership == "National_Park")


USFS <- point_select %>%
  filter(ownership == "Forest_Service")



Private <- point_select %>%
  filter(ownership == "Private")


A <- ggplot(test, aes(mean_prob, mean_SDI, color = ownership)) +
  geom_point() +
  geom_errorbar(aes(ymin= mean_SDI-sd_SDI, ymax= mean_SDI+sd_SDI, color = factor(ownership)), width = 0) +
  geom_errorbar(aes(xmin = mean_prob-sd_prob, xmax= mean_prob+sd_prob, color = factor(ownership)),  width = 0) +
  scale_color_manual("ownership", values = c("darkolivegreen4", "coral3", "blue4"), 
                     labels = c("National Park", "National Forest", "Private")) +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.40, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 01) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  guides(color = guide_legend(title = "Ownership")) +
  theme_classic()

## add images ##
# load images for quadrants
stop <- readPNG("figures/images/stop.png")
caution <- readPNG("figures/images/caution.png")
axe <- readPNG("figures/images/axe.png")
tree <- readPNG("figures/images/tree.png")

# add to graph 
A <- A +
  annotate("text", x= 0.875, y= 0.63, label = "imminent mortality threshold", color = "red") +
  annotation_custom(rasterGrob(axe, interpolate = TRUE), xmin = 0, xmax = 0.125,
                    ymin = 0.8, ymax = 1.00) +
  annotation_custom(rasterGrob(stop, interpolate = TRUE), xmin = 0, xmax = 0.125,
                    ymin = 0, ymax = 0.125) +
  annotation_custom(rasterGrob(caution, interpolate = TRUE), xmin = 0.875, xmax = 1.00,
                    ymin = 0.8, ymax = 1.00) +
  annotation_custom(rasterGrob(tree, interpolate = TRUE), xmin = 0.875, xmax = 1.00,
                    ymin = 0, ymax = 0.125)



# USFS point density plot 
B <- ggplot(USFS, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "darkolivegreen4") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "darkgreen") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.40, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  ggtitle("National Forest") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none") 

# National Park point density plot 
C <- ggplot(NP, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "coral3") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "orangered4") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.4, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  ggtitle("National Park") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position="none") 

# Private point density plot 
D <- ggplot(Private, aes(Prob_of_use_scaled, Prop_max_SDI)) +
  geom_point(color = "blue4") +
  stat_density_2d(aes(fill = ..level..), geom= "polygon") +
  scale_fill_gradient(low = "cornsilk", high= "darkslateblue") +
  geom_hline(yintercept = 0.60, linetype = "dotted", color = "red", size = 0.75) +
  geom_vline(xintercept = 0.4, linetype = "solid", color = "gray45", size = .5) +
  geom_hline(yintercept = 0.35, linetype = "solid", color = "gray45", size = .5) +
  ylim(ymin= 0, ymax= 1) +
  xlim(xmin= 0, xmax= 1) +
  ylab("Proportion of Max SDI") +
  xlab("Probability of Use") +
  ggtitle("Private Land") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none") 

bottom <- B | C | D

full <- A / bottom

tiff("figures/SDI_Use_Ownership.tiff", units="in", width=10, height=10, res=300)
full
dev.off()
