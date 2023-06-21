# Script name: 2_xSF_Data_Prep.R
# Author: M. E. Wright
# Date created: November 7, 2022
#

# Notes: 
#      This script does the following:
#       1) Looks at the spatial extent of all points
#       2) Subsets finer grain points to more coarse grain to expand data 
#       3) Draws available points for each xSF analysis 
#       4) Extracts covariate values at points to create xSF-ready data
# ---

# ---
# 1. Loading packages and data ----
# ---
library(tidyverse)
library(dplyr)
library(sf)
library(raster)
library(amt)
library(lubridate)


# load in cleaned GPS data
CSO <- readRDS("data/gps/gps_cleaned")

# get rid of unneeded columns 
CSO <- CSO %>%
  select(-c(date_gmt, time_gmt, attach_type, fix, altitude, duration, dop, 
            speed_knots, search_time, voltage, speed_km_hr, azimuth_dtn,
            pdop, hdop, vdop, sats_in_view, snr_avg_db, delta_s, e_res, 
            is_king_fire_gavin, is_king_fire_southern, is_sierra_nevada_tagging,
            is_sierra_demography, is_lassen_demography, is_yosemite, 
            is_sequoia_kings, telemetry_model))

# separate into different dataframes based on sampling frequency
CSO.night <- CSO %>%
  filter(which.xsf == "rsf")
CSO.hour <- CSO %>%
  filter(which.xsf == "ssf.hr")
CSO.min <- CSO %>%
  filter(which.xsf == "ssf.min")

# Jessalyn already assigned a couple birds from 2017 and 2018 to 'rsf'. To 
# facilitate the way I am going to subsample datasets, we need to remove these 
# birds from our starting dataset (ONLY APPLIES TO THE NIGHTLY HSF DATA)
CSO.night <- CSO.night %>%
  filter(Year < 2017)

# ---
# 2. Look at the spatial extent of points in each category ----
# ---
# load in the Sierra Nevada boundary shapefile
SNC <- st_read("data/shapefiles/SNC_Boundary.shp")

# project to NAD83 UTM 10
SNC.proj <- st_transform(SNC, crs = 32610)

# create spatial object
SNC.sp <- as_Spatial(SNC.proj)
plot(SNC.sp)


# create sf object for owls
owl.sf <- st_as_sf(CSO, coords = c("long_utm", "lat_utm"), crs = 32610)

# filter into unique sampling units 
owl.night <- owl.sf %>%
  filter(which.xsf == "rsf" & Year < 2017)
owl.hour <- owl.sf %>%
  filter(which.xsf == "ssf.hr")
owl.min <- owl.sf %>%
  filter(which.xsf == "ssf.min")

# create spatial objects 
owl.nightSP <- as_Spatial(owl.night)
owl.hourSP <- as_Spatial(owl.hour)
owl.minSP <- as_Spatial(owl.min)

# plot each 
plot(SNC.sp)
plot(owl.nightSP, add = TRUE, pch = 1, col = "purple", lwd = 2)
plot(owl.hourSP, add = TRUE, pch = 0, col = "cyan3", lwd = 2)
plot(owl.minSP, add = TRUE, pch = 2, col = "azure4", lwd = 0.5)


# ---
# 3. Extract covariate data for minute-to-minute telemetry points ----
# ---
# The minute-to-minute data should be the most straight-forward since it does 
# not need to be resampled

# select only the important variables 
CSO.minSUB <- CSO.min[ , c("Territory", "pac_id", "usfws_id", "sex", "land_owner",
                           "forest", "data.origin", "datetime", "long_utm", 
                           "lat_utm", "Year")]

# check for duplicated rows (remove as needed)
dupcheck <- duplicated(CSO.minSUB)
sum(dupcheck, na.rm = TRUE)

# nest dataframes by territory and year 
CSO.minNEST <- nest(CSO.minSUB, data = -c(Territory, Year, usfws_id, pac_id, sex,
                                          forest, land_owner))

# use 'lapply' to make tracks 
CSO.minNEST <- CSO.minNEST %>%
  mutate(track = lapply(data, function(d) {
    res <- make_track(d, long_utm, lat_utm, datetime, crs = 32610)
  }))

# resample track to common time 
# set to every 5 minutes with a 2 minute tolerance (range: every 3-7 minutes)
CSO.minNEST <- CSO.minNEST %>%
  mutate(track_standard = lapply(track, function(track) {
    res <- track_resample(track, rate = minutes(5), tolerance = minutes(2))
  }))

# create steps by burst 
CSO.minNEST <- CSO.minNEST %>%
  mutate(steps = lapply(track_standard, function(track_standard) {
    res <- steps_by_burst(track_standard)
  }))

# create random steps 
# 20 random steps for each observed 
# set random seed to obtain the same random steps each time 
set.seed(20221107)
CSO.minNEST <- CSO.minNEST %>%
  mutate(random_steps = lapply(steps, function(steps) {
    res <- random_steps(steps, n_control = 20)
  }))


# check unique years for this data set (* 2019 and 2020 *)
unique(CSO.min$Year)


# for extracting covariates, separate out into 2019 and 2020 dataframes 
CSO.minNEST2019 <- CSO.minNEST %>%
  filter(Year == 2019)
CSO.minNEST2020 <- CSO.minNEST %>%
  filter(Year == 2020)

# load in the 2019 and 2020 CFO data 
# 2019 
CBH19 <- raster("data/CFO/2019/CFO-California-CanopyBaseHeight-Summer2019-00010m.tif")
CBD19 <- raster("data/CFO/2019/CFO-California-CanopyBulkDensity-Summer2019-00010m.tif")
CLC19 <- raster("data/CFO/2019/CFO-California-CanopyLayerCount-Summer2019-00010m.tif")
LFD19 <- raster("data/CFO/2019/CFO-California-LadderFuelDensity-Summer2019-00010m.tif")
SF19 <- raster("data/CFO/2019/CFO-California-SurfaceFuels-Summer2019-00010m.tif")

# 2020
CBH20 <- raster("data/CFO/2020/CFO-California-CanopyBaseHeight-Summer2020-00010m.tif")
CBD20 <- raster("data/CFO/2020/CFO-California-CanopyBulkDensity-Summer2020-00010m.tif")
CLC20 <- raster("data/CFO/2020/CFO-California-CanopyLayerCount-Summer2020-00010m.tif")
LFD20 <- raster("data/CFO/2020/CFO-California-LadderFuelDensity-Summer2020-00010m.tif")
SF20 <- raster("data/CFO/2020/CFO-California-SurfaceFuels-Summer2020-00010m.tif")


# extract covariate values for the random points (used points are included)
# repeat these steps for ALL layers
CSO.minNEST2019 <- CSO.minNEST2019 %>%
  mutate(CBH = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CBH19, where = "end")
  }))
CSO.minNEST2019 <- CSO.minNEST2019 %>%
  mutate(CBD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CBD19, where = "end")
  }))
CSO.minNEST2019 <- CSO.minNEST2019 %>%
  mutate(CLC = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CLC19, where = "end")
  }))
CSO.minNEST2019 <- CSO.minNEST2019 %>%
  mutate(LFD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LFD19, where = "end")
  }))
CSO.minNEST2019 <- CSO.minNEST2019 %>%
  mutate(SF = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, SF19, where = "end")
  }))



CSO.minNEST2020 <- CSO.minNEST2020 %>%
  mutate(CBH = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CBH20, where = "end")
  }))
CSO.minNEST2020 <- CSO.minNEST2020 %>%
  mutate(CBD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CBD20, where = "end")
  }))
CSO.minNEST2020 <- CSO.minNEST2020 %>%
  mutate(CLC = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, CLC20, where = "end")
  }))
CSO.minNEST2020 <- CSO.minNEST2020 %>%
  mutate(LFD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LFD20, where = "end")
  }))
CSO.minNEST2020 <- CSO.minNEST2020 %>%
  mutate(SF = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, SF20, where = "end")
  }))


# unnest/nest to build dataframes with only the important covariate data 
# IS THERE A BETTER WAY TO DO THIS? -> YES!
# DO I KNOW WHAT IT IS?! -> NO!

# pull out the identifying and covariate data and unnest
cov19 <- CSO.minNEST2019 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                forest, land_owner, CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

# change covariate names for easier manipulation
cov19 <- cov19 %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2019.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2019.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2019.00010m,
         SF = CFO.California.SurfaceFuels.Summer2019.00010m)

# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov19$CBD)
CLC <- as.vector(cov19$CLC)
LFD <- as.vector(cov19$LFD)
SF <- as.vector(cov19$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.min2019 <- CSO.minNEST2019 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2019.00010m)

# add in the rest of the covariate vectors 
CSO.min2019$CanopyBulkDensity <- CBD
CSO.min2019$CanopyLayerCount <- CLC
CSO.min2019$LadderFuelDensity <- LFD
CSO.min2019$SurfaceFuels <- SF



# SAME for 2020!
# pull out the identifying and covariate data and unnest
cov20 <- CSO.minNEST2020 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

# change covariate names for easier manipulation
cov20 <- cov20 %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2020.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2020.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2020.00010m,
         SF = CFO.California.SurfaceFuels.Summer2020.00010m)

# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov20$CBD)
CLC <- as.vector(cov20$CLC)
LFD <- as.vector(cov20$LFD)
SF <- as.vector(cov20$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.min2020 <- CSO.minNEST2020 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2020.00010m)

# add in the rest of the covariate vectors 
CSO.min2020$CanopyBulkDensity <- CBD
CSO.min2020$CanopyLayerCount <- CLC
CSO.min2020$LadderFuelDensity <- LFD
CSO.min2020$SurfaceFuels <- SF


# bind the 2019 and 2020 data back together 
CSO.min <- rbind(CSO.min2019, CSO.min2020)

# add columns for log of step length and cosine of turn angle 
CSO.min <- CSO.min %>%
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

## SAVE FINAL MINUTE SSF DATA ##
saveRDS(CSO.min, file = "data/xSF_data/minute_ssf")



# ---
# 4. Extract covariate data for hourly telemetry points ----
# ---
# For the hourly dataset, we will follow the same data prep steps that we did for
# the minute-to-minute SSfs. I attempted to resample the minute-to-minute data 
# to add in to the hourly, but the minute-to-minute data does not cover enough
# hours per night to be added to this sample.

# select only the important columns
CSO.hourSUB <- CSO.hour[ , c("Territory", "pac_id", "usfws_id", "sex", "land_owner",
                          "forest", "data.origin", "datetime", "long_utm", 
                          "lat_utm", "Year")]

# check for duplicated rows (remove as needed)
dupcheck <- duplicated(CSO.hourSUB)
sum(dupcheck, na.rm = TRUE)


# nest dataframes by territory and year 
CSO.hrNEST <- nest(CSO.hourSUB, data = -c(Territory, Year, usfws_id, pac_id, sex,
                                          forest, land_owner))

# use 'lapply' to make tracks 
CSO.hrNEST <- CSO.hrNEST %>%
  mutate(track = lapply(data, function(d) {
    res <- make_track(d, long_utm, lat_utm, datetime, crs = 32610)
  }))

# resample track to common time 
# set to every 1 hour with a 15 minute tolerance (range: 45 minutes - 1:15)
CSO.hrNEST <- CSO.hrNEST %>%
  mutate(track_standard = lapply(track, function(track) {
    res <- track_resample(track, rate = minutes(60), tolerance = minutes(15))
  }))

# create steps by burst 
CSO.hrNEST <- CSO.hrNEST %>%
  mutate(steps = lapply(track_standard, function(track_standard) {
    res <- steps_by_burst(track_standard)
  }))

# create random steps 
# 20 random steps for each observed 
# set random seed to obtain the same random steps each time 
set.seed(20221107)
CSO.hrNEST <- CSO.hrNEST %>%
  mutate(random_steps = lapply(steps, function(steps) {
    res <- random_steps(steps, n_control = 20)
  }))


# look at the unique years in this dataset for extracting covariate values
unique(CSO.hrNEST$Year)


# for extracting covariates, separate out into 2017 and 2018 dataframes 
CSO.hrNEST2017 <- CSO.hrNEST %>%
  filter(Year == 2017)
CSO.hrNEST2018 <- CSO.hrNEST %>%
  filter(Year == 2018)

# load in the 2017 and 2018 CFO data 
# 2017
BH.2017 <- raster("data/CFO/2017/CFO-California-CanopyBaseHeight-Summer2017-00010m.tif")
BD.2017 <- raster("data//CFO/2017/CFO-California-CanopyBulkDensity-Summer2017-00010m.tif")
LC.2017 <- raster("data//CFO/2017/CFO-California-CanopyLayerCount-Summer2017-00010m.tif")
LFD.2017 <- raster("data/CFO/2017/CFO-California-LadderFuelDensity-Summer2017-00010m.tif")
SF.2017 <- raster("data/CFO/2017/CFO-California-SurfaceFuels-Summer2017-00010m.tif")

# 2018
BH.2018 <- raster("data/CFO/2018/CFO-California-CanopyBaseHeight-Summer2018-00010m.tif")
BD.2018 <- raster("data/CFO/2018/CFO-California-CanopyBulkDensity-Summer2018-00010m.tif")
LC.2018 <- raster("data/CFO/2018/CFO-California-CanopyLayerCount-Summer2018-00010m.tif")
LFD.2018 <- raster("data/CFO/2018/CFO-California-LadderFuelDensity-Summer2018-00010m.tif")
SF.2018 <- raster("data/CFO/2018/CFO-California-SurfaceFuels-Summer2018-00010m.tif")


# extract covariate values for the random points (used points are included)
CSO.hrNEST2017 <- CSO.hrNEST2017 %>%
  mutate(CBH = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, BH.2017, where = "end")
  }))
CSO.hrNEST2017 <- CSO.hrNEST2017 %>%
  mutate(CBD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, BD.2017, where = "end")
  }))
CSO.hrNEST2017 <- CSO.hrNEST2017 %>%
  mutate(CLC = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LC.2017, where = "end")
  }))
CSO.hrNEST2017 <- CSO.hrNEST2017 %>%
  mutate(LFD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LFD.2017, where = "end")
  }))
CSO.hrNEST2017 <- CSO.hrNEST2017 %>%
  mutate(SF = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, SF.2017, where = "end")
  }))




CSO.hrNEST2018 <- CSO.hrNEST2018 %>%
  mutate(CBH = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, BH.2018, where = "end")
  }))
CSO.hrNEST2018 <- CSO.hrNEST2018 %>%
  mutate(CBD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, BD.2018, where = "end")
  }))
CSO.hrNEST2018 <- CSO.hrNEST2018 %>%
  mutate(CLC = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LC.2018, where = "end")
  }))
CSO.hrNEST2018 <- CSO.hrNEST2018 %>%
  mutate(LFD = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, LFD.2018, where = "end")
  }))
CSO.hrNEST2018 <- CSO.hrNEST2018 %>%
  mutate(SF = lapply(random_steps, function(random_steps) {
    res <- extract_covariates(random_steps, SF.2018, where = "end")
  }))


# pull out the identifying and covariate data and unnest
cov17 <- CSO.hrNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

# change covariate names for easier manipulation
cov17 <- cov17 %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2017.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2017.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2017.00010m,
         SF = CFO.California.SurfaceFuels.Summer2017.00010m)

# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov17$CBD)
CLC <- as.vector(cov17$CLC)
LFD <- as.vector(cov17$LFD)
SF <- as.vector(cov17$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.hr2017 <- CSO.hrNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2017.00010m)

# add in the rest of the covariate vectors 
CSO.hr2017$CanopyBulkDensity <- CBD
CSO.hr2017$CanopyLayerCount <- CLC
CSO.hr2017$LadderFuelDensity <- LFD
CSO.hr2017$SurfaceFuels <- SF


# pull out the identifying and covariate data and unnest
cov18 <- CSO.hrNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

# change covariate names for easier manipulation
cov18 <- cov18 %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2018.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2018.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2018.00010m,
         SF = CFO.California.SurfaceFuels.Summer2018.00010m)

# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov18$CBD)
CLC <- as.vector(cov18$CLC)
LFD <- as.vector(cov18$LFD)
SF <- as.vector(cov18$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.hr2018 <- CSO.hrNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2018.00010m)

# add in the rest of the covariate vectors 
CSO.hr2018$CanopyBulkDensity <- CBD
CSO.hr2018$CanopyLayerCount <- CLC
CSO.hr2018$LadderFuelDensity <- LFD
CSO.hr2018$SurfaceFuels <- SF


# bind the 2019 and 2020 data back together 
CSO.hr <- rbind(CSO.hr2017, CSO.hr2018)

# add columns for log of step length and cosine of turn angle 
CSO.hr <- CSO.hr %>%
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

## SAVE FINAL MINUTE SSF DATA ##
saveRDS(CSO.hr, file = "data/xSF_data/hour_ssf")



# ---
# 5. Extract covariate data for nightly telemetry points ----
# ---
# For the nightly dataset, we can resample the minute-to-minute and hourly in 
# order to increase the number of points for the RSF analysis.


# we can use the CSO.minSUB and CSO.hourSUB dataframes to start working through 
# this.

# use lubridate to create a unique column for each day 
CSO.minSUB <- CSO.minSUB %>%
  mutate(day = day(datetime))
CSO.hourSUB <- CSO.hourSUB %>%
  mutate(day = day(datetime))

# now we can group by territory/year and select one random location per night
set.seed(20221130)
CSO.minSUB2 <- CSO.minSUB %>%
 group_by(Territory, Year, day) %>%
 slice_sample(n = 1)
CSO.hourSUB2 <- CSO.hourSUB %>%
  group_by(Territory, Year, day) %>%
  slice_sample(n = 1)

# select only the important columns from each 
CSO.minSUB2 <- CSO.minSUB2[ , c("Territory", "pac_id", "usfws_id", "sex", "land_owner",
                                "forest", "data.origin", "datetime", "long_utm", 
                                "lat_utm", "Year")]
CSO.hourSUB2 <- CSO.hourSUB2[ , c("Territory", "pac_id", "usfws_id", "sex", "land_owner",
                                  "forest", "data.origin", "datetime", "long_utm", 
                                  "lat_utm", "Year")]

# add in a column to designate which original dataset these came from 
CSO.minSUB2 <- CSO.minSUB2 %>%
  mutate(data_source = "minute")
CSO.hourSUB2 <- CSO.hourSUB2 %>%
  mutate(data_source = "hour")

# bind together each subset dataframe
CSO.subset <- rbind(CSO.minSUB2, CSO.hourSUB2)

# check for duplicated rows (remove as needed)
dupcheck <- duplicated(CSO.subset)
sum(dupcheck, na.rm = TRUE)

# now select only the important columns from the nightly data 
CSO.nightSUB <- CSO.night[ , c("Territory", "pac_id", "usfws_id", "sex", "land_owner",
                               "forest", "data.origin", "datetime", "long_utm", 
                               "lat_utm", "Year")]

# add in the data_source column 
CSO.nightSUB <- CSO.nightSUB %>%
  mutate(data_source = "night")

# bind together all data 
CSO.night <- rbind(CSO.subset, CSO.nightSUB)

# check for duplicated rows (remove as needed)
dupcheck <- duplicated(CSO.night)
sum(dupcheck, na.rm = TRUE)


# nest dataframes by territory and year 
CSO.nightNEST <- nest(CSO.night, data = -c(Territory, Year, usfws_id, pac_id, sex,
                                           forest, land_owner, data_source))

# use 'lapply' to make tracks 
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(track = lapply(data, function(d) {
    res <- make_track(d, long_utm, lat_utm, datetime, crs = 32610)
  }))

# minute-to-minute data may not work for nightly RSFs because there are too few
# points (minute-to-minute data collected over a very short number of days)
# filter out minute data and try again
CSO.nightNEST <- CSO.nightNEST %>%
  filter(Year != 2019 & Year != 2020)


# Create a list of counted number of points for each bird 
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(n_pts = lapply(data, function(data) {
    res <- nrow(data)
  }))

# unnest the n_pts data 
CSO.nightNEST <- CSO.nightNEST %>%
  unnest(n_pts) 

# select only birds for which we have <20 pts to move to the HSF
CSO.nightNEST <- CSO.nightNEST %>%
  filter(n_pts >= 20)

# look at the mean number of used points (mean = 49 pts)
mean(unlist(CSO.nightNEST$n_pts))


# fit 95% akde home range
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(akde = lapply(track, function(track) {
    res <- hr_akde(track, levels = 0.95)
  })) 


# look at the number of birds left from subsampling and number of original 
# nighlty telemetry points birds
CSO.nightNEST %>%
  group_by(data_source, Year) %>%
  summarize(n = n())

# look at the mean number of points for each group 
CSO.nightNEST %>%
  group_by(data_source) %>%
  summarize(mean = mean(n_pts))

# to generate the correct number of available points, we need to create a column
# that has the correct number for each bird 
# Muff et al. (2020) uses a 33:66 ratio, so we will copy this
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(n_available = lapply(n_pts, function(n_pts) {
    res <- n_pts * 2 
  }))

# unnest the n_available data 
CSO.nightNEST <- CSO.nightNEST %>%
  unnest(n_available) 


# generate random points using the n_available column as a template for number 
# of available points
set.seed(20221107)
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(random_pts = lapply(akde, function(akde) {
    res <- random_points(akde, n = n_available, type = "random")
  }))

# check the number of random points drawn is equal to the number of available
# points that we calculated 
CSO.nightNEST <- CSO.nightNEST %>%
  mutate(available_pts_CHECK = lapply(random_pts, function(random_pts) {
    res <- nrow(random_pts)
  }))

################################################################################
# FIXING # OF RANDOM POINTS ---- 
# for some reason, the n_available works for MOST birds, but there are some that
# have way more points than needed 
# let's try looking at the akde home range and random points for one of these 
# birds 
# choose line #61 from the nested dataframe
test <- CSO.nightNEST[68, ]

# unnest data and random points 
test.used <- test %>%
  unnest(cols = data) %>%
  dplyr::select(c(Territory, usfws_id, Year, long_utm, lat_utm)) %>%
  mutate(case_ = "TRUE") %>%
  rename(x_ = long_utm, y_ = lat_utm)
test.random <- test %>%
  unnest(cols = random_pts) %>%
  dplyr::select(c(Territory, usfws_id, Year, x_, y_, case_))

# convert case_ in test.random to character
test.random$case_ <- as.character(test.random$case_)

# create sf object 
testrandom.sf <- st_as_sf(test.random, coords = c("x_", "y_"), crs = 32610)
testused.sf <- st_as_sf(test.used, coords = c("x_", "y_"), crs = 32610)
  
# as spatial 
testrandom.sp <- as_Spatial(testrandom.sf)
testused.sp <- as_Spatial(testused.sf)

# unlist akde 
test.akde <- unlist(test$akde)

# pull the rasterlayer
test.akde <- test.akde[[1]]

# plot 
plot(test.akde)
plot(testrandom.sp, add = TRUE)
plot(testused.sp, add = TRUE, pch = 0, col = "cyan3", lwd = 2)


# I have no idea why the n_available works for some birds and not for others, so
# I will try to fix it manually for now (10/80 birds have this issue)
## TURNS OUT THIS IS A LOT EASIER AT THE END OF THIS SCRIPT!!!! ## 

################################################################################

# RETURN TO CODING HSF DATA PREP ---- 

# look at the unique years in this dataset for extracting covariate values
unique(CSO.nightNEST$Year)

# for extracting covariates, separate out into 2015, 2016, 2017, and 2018
# dataframes
# ** there is no CFO data for 2015, but we will use 2016 data for pulling 
#    covariate values **
CSO.nightNEST2016 <- CSO.nightNEST %>%
  filter(Year == 2016 | Year == 2015)
CSO.nightNEST2017 <- CSO.nightNEST %>%
  filter(Year == 2017)
CSO.nightNEST2018 <- CSO.nightNEST %>%
  filter(Year == 2018)


# should already have Rasterlayers for 2017 and 2018, so load in only 2016 data 
# 2016 
BH.2016 <- raster("data/CFO/2016/CFO-California-CanopyBaseHeight-Summer2016-00010m.tif")
BD.2016 <- raster("data/CFO/2016/CFO-California-CanopyBulkDensity-Summer2016-00010m.tif")
LC.2016 <- raster("data/CFO/2016/CFO-California-CanopyLayerCount-Summer2016-00010m.tif")
LFD.2016 <- raster("data/CFO/2016/CFO-California-LadderFuelDensity-Summer2016-00010m.tif")
SF.2016 <- raster("data/CFO/2016/CFO-California-SurfaceFuels-Summer2016-00010m.tif")


# extract covariate values for the random points (used points are NOT included)
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CBH = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BH.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CBD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BD.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CLC = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LC.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(LFD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LFD.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(SF = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, SF.2016)
  }))



CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CBH = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BH.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CBD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BD.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CLC = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LC.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(LFD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LFD.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(SF = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, SF.2017)
  }))



CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CBH = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BH.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CBD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, BD.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CLC = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LC.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(LFD = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, LFD.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(SF = lapply(random_pts, function(random_pts) {
    res <- extract_covariates(random_pts, SF.2018)
  }))



# do the same steps for the tracks (covariate values at used points)
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CBH.track = lapply(track, function(track) {
    res <- extract_covariates(track, BH.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CBD.track = lapply(track, function(track) {
    res <- extract_covariates(track, BD.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(CLC.track = lapply(track, function(track) {
    res <- extract_covariates(track, LC.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(LFD.track = lapply(track, function(track) {
    res <- extract_covariates(track, LFD.2016)
  }))
CSO.nightNEST2016 <- CSO.nightNEST2016 %>%
  mutate(SF.track = lapply(track, function(track) {
    res <- extract_covariates(track, SF.2016)
  }))



CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CBH.track = lapply(track, function(track) {
    res <- extract_covariates(track, BH.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CBD.track = lapply(track, function(track) {
    res <- extract_covariates(track, BD.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(CLC.track = lapply(track, function(track) {
    res <- extract_covariates(track, LC.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(LFD.track = lapply(track, function(track) {
    res <- extract_covariates(track, LFD.2017)
  }))
CSO.nightNEST2017 <- CSO.nightNEST2017 %>%
  mutate(SF.track = lapply(track, function(track) {
    res <- extract_covariates(track, SF.2017)
  }))



CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CBH.track = lapply(track, function(track) {
    res <- extract_covariates(track, BH.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CBD.track = lapply(track, function(track) {
    res <- extract_covariates(track, BD.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(CLC.track = lapply(track, function(track) {
    res <- extract_covariates(track, LC.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(LFD.track = lapply(track, function(track) {
    res <- extract_covariates(track, LFD.2018)
  }))
CSO.nightNEST2018 <- CSO.nightNEST2018 %>%
  mutate(SF.track = lapply(track, function(track) {
    res <- extract_covariates(track, SF.2018)
  }))



## AVAILABLE POINTS ##
# pull out the identifying and covariate data and unnest (AVAILABLE POINTS)
cov16AVAILABLE <- CSO.nightNEST2016 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

cov17AVAILABLE <- CSO.nightNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available, 
                  available_pts_CHECK,
                  CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")

cov18AVAILABLE <- CSO.nightNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH, CBD, CLC, LFD, SF)) %>%
  unnest(cols = c(CBD, CLC, LFD, SF), names_repair = "universal")


# change covariate names for easier manipulation
cov16AVAILABLE <- cov16AVAILABLE %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2016.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2016.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2016.00010m,
         SF = CFO.California.SurfaceFuels.Summer2016.00010m)

cov17AVAILABLE <- cov17AVAILABLE %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2017.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2017.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2017.00010m,
         SF = CFO.California.SurfaceFuels.Summer2017.00010m)

cov18AVAILABLE <- cov18AVAILABLE %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2018.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2018.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2018.00010m,
         SF = CFO.California.SurfaceFuels.Summer2018.00010m)


################################################################################
## DO THE NEXT STEPS IN SEQUENCE FOR EACH YEAR-SPECIFIC DATAFRAME ## 
## 2016 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov16AVAILABLE$CBD)
CLC <- as.vector(cov16AVAILABLE$CLC)
LFD <- as.vector(cov16AVAILABLE$LFD)
SF <- as.vector(cov16AVAILABLE$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2016AVAILABLE <- CSO.nightNEST2016 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK, 
                  CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2016.00010m)

# add in the rest of the covariate vectors 
CSO.night2016AVAILABLE$CanopyBulkDensity <- CBD
CSO.night2016AVAILABLE$CanopyLayerCount <- CLC
CSO.night2016AVAILABLE$LadderFuelDensity <- LFD
CSO.night2016AVAILABLE$SurfaceFuels <- SF


## 2017 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov17AVAILABLE$CBD)
CLC <- as.vector(cov17AVAILABLE$CLC)
LFD <- as.vector(cov17AVAILABLE$LFD)
SF <- as.vector(cov17AVAILABLE$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2017AVAILABLE <- CSO.nightNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2017.00010m)

# add in the rest of the covariate vectors 
CSO.night2017AVAILABLE$CanopyBulkDensity <- CBD
CSO.night2017AVAILABLE$CanopyLayerCount <- CLC
CSO.night2017AVAILABLE$LadderFuelDensity <- LFD
CSO.night2017AVAILABLE$SurfaceFuels <- SF


## 2018 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov18AVAILABLE$CBD)
CLC <- as.vector(cov18AVAILABLE$CLC)
LFD <- as.vector(cov18AVAILABLE$LFD)
SF <- as.vector(cov18AVAILABLE$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2018AVAILABLE <- CSO.nightNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH)) %>%
  unnest(cols = CBH) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2018.00010m)

# add in the rest of the covariate vectors 
CSO.night2018AVAILABLE$CanopyBulkDensity <- CBD
CSO.night2018AVAILABLE$CanopyLayerCount <- CLC
CSO.night2018AVAILABLE$LadderFuelDensity <- LFD
CSO.night2018AVAILABLE$SurfaceFuels <- SF

################################################################################

# bind together all the AVAILABLE data 
CSO.nightAVAILABLE <- rbind(CSO.night2016AVAILABLE, CSO.night2017AVAILABLE, 
                            CSO.night2018AVAILABLE)


################################################################################
## FIX THE RANDOM POINTS TO HAVE 33:66 RATIO ----
# need to do a manual fix since I can't figure out why the random_points code 
# is not working
CSO.NESTMANUALFIX <- CSO.nightAVAILABLE %>%
 filter(n_available != available_pts_CHECK)

# check that there are 10 birds in this group 
CSO.NESTMANUALFIX %>%
  group_by(Territory, Year, usfws_id) %>%
  summarize(n = n())

# still not sure what is going on with this, so I am going to just take a 
# random sample of the generated random points
# I am not sure how to do this in a more elegant way, so I am going to go 
# through it for each bird, messy and inefficient 
bird1 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2015 & usfws_id == 138784884)
bird2 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2015 & usfws_id == 180749189)
bird3 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 138784835)
bird4 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 195723204)
bird5 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2015 & usfws_id == 117722074)
bird6 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2015 & usfws_id == 138758120)
bird7 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 117722074)
bird8 <-  CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 138758120)
bird9 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 195723202)
bird10 <- CSO.NESTMANUALFIX %>%
  filter(Year == 2016 & usfws_id == 195723203)

# check for duplicated rows (remove as needed) -- random duplicate check 
dupcheck <- duplicated(bird9)
sum(dupcheck, na.rm = TRUE)

# unnest random points for each bird and then slice the sample based on the 
# n_available column number
set.seed(20221201)

bird1 <- bird1 %>%
  slice_sample(n = 244)
bird2 <- bird2 %>%
  slice_sample(n = 206)
bird3 <- bird3 %>%
  slice_sample(n = 176)
bird4 <- bird4 %>%
  slice_sample(n = 206)
bird5 <- bird5 %>%
  slice_sample(n = 182)
bird6 <- bird6 %>%
  slice_sample(n = 236)
bird7 <- bird7 %>%
  slice_sample(n = 162)
bird8 <- bird8 %>%
  slice_sample(n = 166)
bird9 <- bird9 %>%
  slice_sample(n = 222)
bird10 <- bird10 %>%
  slice_sample(n = 254)

# bind these together 
problembirds <- rbind(bird1, bird2, bird3, bird4, bird5, bird6, bird7, bird8, 
                      bird9, bird10)
# reduce the CSO.nightAVAILABLE dataframe down to only the birds where 
# n_available and available_pts_CHECK are in agreement
CSO.nightAVAILABLE2 <- CSO.nightAVAILABLE %>%
  filter(n_available == available_pts_CHECK)

# add in the cleaned problem birds 
CSO.nightAVAILABLEclean <- rbind(CSO.nightAVAILABLE2, problembirds)

# drop the available_pts_CHECK column (no longer accurate)
CSO.nightAVAILABLEclean <- CSO.nightAVAILABLEclean %>%
  dplyr::select(-c(available_pts_CHECK))

# nest dataframe so that we can do another available points check 
CSO.nightAVAILABLENEST <- nest(CSO.nightAVAILABLEclean, 
                               random_pts = c(case_, x_, y_, CanopyBaseHeight,
                                              CanopyBulkDensity, CanopyLayerCount,
                                              LadderFuelDensity, SurfaceFuels))

# check the number of random points drawn is equal to the number of available
# points that we calculated 
CSO.nightAVAILABLENEST <- CSO.nightAVAILABLENEST %>%
  mutate(available_pts_CHECK = lapply(random_pts, function(random_pts) {
    res <- nrow(random_pts)
  }))

## NOW THEY MATCH!!! ##
################################################################################

# unnest that random points data 
CSO.nightAVAILABLE <- unnest(CSO.nightAVAILABLENEST, cols = random_pts)


## USED POINTS ##
# pull out the identifying and covariate data and unnest (USED POINTS)
cov16USED <- CSO.nightNEST2016 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH.track, CBD.track, CLC.track, LFD.track, SF.track)) %>%
  unnest(cols = c(CBD.track, CLC.track, LFD.track, SF.track),
         names_repair = "universal")

cov17USED <- CSO.nightNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH.track, CBD.track, CLC.track, LFD.track, SF.track)) %>%
  unnest(cols = c(CBH.track, CBD.track, CLC.track, LFD.track, SF.track), 
         names_repair = "universal")

cov18USED <- CSO.nightNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,
                  CBH.track, CBD.track, CLC.track, LFD.track, SF.track)) %>%
  unnest(cols = c(CBH.track, CBD.track, CLC.track, LFD.track, SF.track), 
         names_repair = "universal")


# change covariate names for easier manipulation
cov16USED<- cov16USED %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2016.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2016.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2016.00010m,
         SF = CFO.California.SurfaceFuels.Summer2016.00010m)

cov17USED <- cov17USED %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2017.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2017.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2017.00010m,
         SF = CFO.California.SurfaceFuels.Summer2017.00010m)

cov18USED <- cov18USED %>%
  rename(CBD = CFO.California.CanopyBulkDensity.Summer2018.00010m,
         CLC = CFO.California.CanopyLayerCount.Summer2018.00010m,
         LFD = CFO.California.LadderFuelDensity.Summer2018.00010m,
         SF = CFO.California.SurfaceFuels.Summer2018.00010m)


################################################################################
## DO THE NEXT STEPS IN SEQUENCE FOR EACH YEAR-SPECIFIC DATAFRAME ## 
## 2016 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov16USED$CBD)
CLC <- as.vector(cov16USED$CLC)
LFD <- as.vector(cov16USED$LFD)
SF <- as.vector(cov16USED$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2016USED <- CSO.nightNEST2016 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK, CBH.track)) %>%
  unnest(cols = CBH.track) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2016.00010m)

# add in the rest of the covariate vectors 
CSO.night2016USED$CanopyBulkDensity <- CBD
CSO.night2016USED$CanopyLayerCount <- CLC
CSO.night2016USED$LadderFuelDensity <- LFD
CSO.night2016USED$SurfaceFuels <- SF


## 2017 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov17USED$CBD)
CLC <- as.vector(cov17USED$CLC)
LFD <- as.vector(cov17USED$LFD)
SF <- as.vector(cov17USED$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2017USED <- CSO.nightNEST2017 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK, CBH.track)) %>%
  unnest(cols = CBH.track) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2017.00010m)

# add in the rest of the covariate vectors 
CSO.night2017USED$CanopyBulkDensity <- CBD
CSO.night2017USED$CanopyLayerCount <- CLC
CSO.night2017USED$LadderFuelDensity <- LFD
CSO.night2017USED$SurfaceFuels <- SF


## 2018 ##
# pull out only the covariate columns and save as vectors (skip CBH for now)
CBD <- as.vector(cov18USED$CBD)
CLC <- as.vector(cov18USED$CLC)
LFD <- as.vector(cov18USED$LFD)
SF <- as.vector(cov18USED$SF)

# create a dataframe with the identifying data and CBH (unnest CBH so we have 
# movement data related to the SSF)
CSO.night2018USED <- CSO.nightNEST2018 %>%
  dplyr::select(c(Territory, Year, usfws_id, pac_id, sex,
                  forest, land_owner, data_source, n_pts, n_available,
                  available_pts_CHECK,CBH.track)) %>%
  unnest(cols = CBH.track) %>%
  rename(CanopyBaseHeight = CFO.California.CanopyBaseHeight.Summer2018.00010m)

# add in the rest of the covariate vectors 
CSO.night2018USED$CanopyBulkDensity <- CBD
CSO.night2018USED$CanopyLayerCount <- CLC
CSO.night2018USED$LadderFuelDensity <- LFD
CSO.night2018USED$SurfaceFuels <- SF

################################################################################

# bind together all the USED data 
CSO.nightUSED <- rbind(CSO.night2016USED, CSO.night2017USED, CSO.night2018USED)


# add a 'case_' column to used and delete the 't_' column
CSO.nightUSED <- CSO.nightUSED %>%
  mutate(case_ = "TRUE") %>%
  dplyr::select(-t_)

# make the 'case_' column logical 
CSO.nightUSED$case_ <- as.logical(CSO.nightUSED$case_)

# bind together used and available points
CSO.night <- rbind(CSO.nightAVAILABLE, CSO.nightUSED)

## SAVE FINAL NIGHT RSF DATA ##
saveRDS(CSO.night, file = "data/xSF_data/night_rsf")



