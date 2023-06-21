# Script name: 6_Predictions.R
# Author: M. E. Wright
# Date created: February 14, 2023
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
library(glmmTMB)
library(magrittr)

# Load in original 'night' dataframe
night <- readRDS("data/xSF_data/night_rsf")

# convert logical column 'case_' to an integer column of 1's and 0's 'STATUS'
night$STATUS <- as.integer(ifelse(night$case_ == "TRUE", 1, 0))

# center and scale all variables
night$scaled_CanopyBaseHeight <- scale(night$CanopyBaseHeight, center = TRUE, scale = TRUE)
night$scaled_CanopyBulkDensity <- scale(night$CanopyBulkDensity, center = TRUE, scale = TRUE)
night$scaled_CanopyLayerCount <- scale(night$CanopyLayerCount, center = TRUE, scale = TRUE)
night$scaled_LadderFuelDensity <- scale(night$LadderFuelDensity, center = TRUE, scale = TRUE)
night$scaled_SurfaceFuels <- scale(night$SurfaceFuels, center = TRUE, scale = TRUE)

# Load in the rsf models (fixed variance mixed effects models)
CBH.rsf <- readRDS("data/model_objects/night_hsf/fixed_variance/CanopyBaseHeight")
CBD.rsf <- readRDS("data/model_objects/night_hsf/fixed_variance/CanopyBulkDensity")
CLC.rsf <- readRDS("data/model_objects/night_hsf/fixed_variance/CanopyLayerCount")
LFD.rsf <- readRDS("data/model_objects/night_hsf/fixed_variance/LadderFuelDensity")
SF.rsf  <- readRDS("data/model_objects/night_hsf/fixed_variance/SurfaceFuels")


# ---
# 2. Predict preference profiles for each variable ----
# ---

# Canopy Base Height (SCALED: min= -1.626, max= 8.690 for RSF data) ----
# make prediction dataframe
newdat.CBH <- expand.grid(CanopyBaseHeight= seq(-2, 9, by= 0.5),
                          ID="dummy", weight=c(0,5000))

# design matrix (fixed effects)
mm.CBH <- model.matrix(delete.response(terms(CBH.rsf)),newdat.CBH)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CBH$STATUS <- plogis(drop(mm.CBH %*% fixef(CBH.rsf)[["cond"]]))

# Skip the var-cov matrix for now (complicated with the high fixed variance in 
# the model)
#predvar <- diag(mm.CBH %*% vcov(CBH.rsf)[["cond"]] %*% t(mm.CBH))
#newdat.CBH$SE <- sqrt(predvar)
#newdat.CBH$SE2 <- sqrt(predvar+sigma(CBH.rsf)^2)

ggplot(newdat.CBH, aes(x = CanopyBaseHeight, y = STATUS)) +
  geom_line() +
  ylim(0, 1)



# Canopy Bulk Density (SCALED: min= -1.486, max= 6.203 for RSF data) ----
# make prediction dataframe
newdat.CBD <- expand.grid(CanopyBulkDensity= seq(from= -1, to= 6, by= 0.5), 
                          ID="dummy", weight=c(0,5000))
# design matrix (fixed effects)
mm.CBD <- model.matrix(delete.response(terms(CBD.rsf)),newdat.CBD)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CBD$STATUS <- plogis(drop(mm.CBD %*% fixef(CBD.rsf)[["cond"]]))


ggplot(newdat.CBD, aes(x = CanopyBulkDensity, y = STATUS)) +
  geom_line()+
  ylim(0, 1)



# Canopy Layer Count (SCALED: min= -2.252, max= 2.844 for RSF data) ----
# make prediction dataframe
newdat.CLC <- expand.grid(CanopyLayerCount= seq(from= -2, to= 3, by= 0.5), 
                          ID="dummy", weight=c(0,5000))
# design matrix (fixed effects)
mm.CLC <- model.matrix(delete.response(terms(CLC.rsf)),newdat.CLC)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CLC$STATUS <- plogis(drop(mm.CLC %*% fixef(CLC.rsf)[["cond"]]))


ggplot(newdat.CLC, aes(x = CanopyLayerCount, y = STATUS)) +
  geom_line() +
  ylim(0, 1)



# Ladder Fuel Density (SCALED: min= -2.099, max= 5.943 for RSF data) ----
# make prediction dataframe
newdat.LFD <- expand.grid(LadderFuelDensity= seq(from= -2, to= 5, by= 0.5), 
                          ID="dummy", weight=c(0,5000))
# design matrix (fixed effects)
mm.LFD <- model.matrix(delete.response(terms(LFD.rsf)),newdat.LFD)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.LFD$STATUS <- plogis(drop(mm.LFD %*% fixef(LFD.rsf)[["cond"]]))


ggplot(newdat.LFD, aes(x = LadderFuelDensity, y = STATUS)) +
  geom_line() +
  ylim(0, 1)



# Surface Fuels (SCALED: min= -3.140, max= 1.359 for RSF data) ----
# make prediction dataframe
newdat.SF <- expand.grid(SurfaceFuel= seq(from= -3, to= 1.5, by= 0.5), 
                         ID="dummy", weight=c(0,5000))
# design matrix (fixed effects)
mm.SF <- model.matrix(delete.response(terms(SF.rsf)),newdat.SF)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.SF$STATUS <- plogis(drop(mm.SF %*% fixef(SF.rsf)[["cond"]]))


ggplot(newdat.SF, aes(x = SurfaceFuel, y = STATUS)) +
  geom_line() +
  ylim(0, 1)



# For each predictive dataframe, backtransform the centered/scaled variables

# create sd and mean objects
sd.CBH <- sd(night$CanopyBaseHeight)
mean.CBH <- mean(night$CanopyBaseHeight)

sd.CBD <- sd(night$CanopyBulkDensity)
mean.CBD <- mean(night$CanopyBulkDensity)

sd.CLC <- sd(night$CanopyLayerCount)
mean.CLC <- mean(night$CanopyLayerCount)

sd.LFD <- sd(night$LadderFuelDensity)
mean.LFD <- mean(night$LadderFuelDensity)

sd.SF <- sd(night$SurfaceFuels)
mean.SF <- mean(night$SurfaceFuels)

newdat.CBH$trans_var <- (newdat.CBH$CanopyBaseHeight * sd.CBH) + mean.CBH


# ---
# 3. Create centered and scaled rasters ----
# ---
# 
# Load in CFO rasters from 2020
CBH <- rast("data/CFO/2020/CFO-California-CanopyBaseHeight-Summer2020-00010m.tif")
CBD <- rast("data/CFO/2020/CFO-California-CanopyBulkDensity-Summer2020-00010m.tif")
CLC <- rast("data/CFO/2020/CFO-California-CanopyLayerCount-Summer2020-00010m.tif")
LFD <- rast("data/CFO/2020/CFO-California-LadderFuelDensity-Summer2020-00010m.tif")
SF <- rast("data/CFO/2020/CFO-California-SurfaceFuels-Summer2020-00010m.tif")


# Prep Sierra Nevada boundary ----
Sierra <- vect("data/shapefiles/SNC_boundary.shp")

# Match the projection
Sierra <- terra::project(Sierra, CBH)



# Canopy base height scaled ----
# mask with Sierra boundary
CBH.mask <- mask(CBH, Sierra)

# Center and scale all raster values
CBH.center <- CBH.mask - mean.CBH
CBH.scale  <- CBH.center / sd.CBH

# save scaled raster
writeRaster(CBH.scale, "data/scaled_rasters/CBH_scaled.tif", overwrite= TRUE)




# Canopy bulk density scaled ----
CBD.mask <- mask(CBD, Sierra)

# Center and scale all raster values
CBD.center <- CBD.mask - mean.CBD
CBD.scale  <- CBD.center / sd.CBD

# save scaled raster
writeRaster(CBD.scale, "data/scaled_rasters/CBD_scaled.tif", overwrite= TRUE)




# Canopy layer count scaled ----
CLC.mask <- mask(CLC, Sierra)

# Center and scale all raster values
CLC.center <- CLC.mask - mean.CLC
CLC.scale  <- CLC.center / sd.CLC

# save scaled raster
writeRaster(CLC.scale, "data/scaled_rasters/CLC_scaled.tif", overwrite = TRUE)




# Ladder fuel density scaled ----
LFD.mask <- mask(LFD, Sierra)

# Center and scale all raster values
LFD.center <- LFD.mask - mean.LFD
LFD.scale  <- LFD.center / sd.LFD

# save scaled raster
writeRaster(LFD.scale, "data/scaled_rasters/LFD_scaled.tif", overwrite = TRUE)




# Surface fuels scaled ----
SF.mask <- mask(SF, Sierra)

# Center and scale all raster values
SF.center <- SF.mask - mean.SF
SF.scale  <- SF.center / sd.SF

# save scaled raster
writeRaster(SF.scale, "data/scaled_rasters/SF_scaled.tif", overwrite = TRUE)

################################################################################

# ---
# 4. Create individual probability rasters ---- 
# ---

# load in scaled rasters (if needed)*********
CBH.scaled <- rast("data/scaled_rasters/CBH_scaled.tif")
CBD.scaled <- rast("data/scaled_rasters/CBD_scaled.tif")
CLC.scaled <- rast("data/scaled_rasters/CLC_scaled.tif")
LFD.scaled <- rast("data/scaled_rasters/LFD_scaled.tif")
SF.scaled  <- rast("data/scaled_rasters/SF_scaled.tif")


# Canopy base height probability raster ----
# pull unique values from scaled raster 
CBH.vals <- unique(CBH.scaled)

# rename 
CBH.vals <- CBH.vals %>%
  rename(CanopyBaseHeight = `CFO-California-CanopyBaseHeight-Summer2020-00010m`)

# use these values in a new predictions dataframe
newdat.CBHraster <- expand.grid(CanopyBaseHeight= CBH.vals$CanopyBaseHeight,
                                ID = "dummy", weight= 5000)

# design matrix (fixed effects)
mm.CBHraster <- model.matrix(delete.response(terms(CBH.rsf)),newdat.CBHraster)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CBHraster$STATUS <- plogis(drop(mm.CBHraster %*% fixef(CBH.rsf)[["cond"]]))

# create a reclassification matrix from the predicted 'STATUS' values 
rclmat.CBH <- matrix(c(CBH.vals$CanopyBaseHeight, newdat.CBHraster$STATUS),
                     ncol= 2)

# reclassify raster
CBH.probability <- classify(CBH.scaled, rclmat.CBH)

# save probability raster #
writeRaster(CBH.probability, "data/probability_rasters/CBH_probability.tif",
            overwrite = TRUE)




# Canopy bulk density probability raster ----
# CBD layer has a LOT of unique values b/c of decimal points... round to only 
# 2 decimal spaces over the entire raster
CBD.round <- round(CBD.scaled, digits = 2)

# pull unique values from raster (omit NA's)
CBD.vals <- unique(CBD.round)

# rename 
CBD.vals <- CBD.vals %>%
  rename(CanopyBulkDensity = `CFO-California-CanopyBulkDensity-Summer2020-00010m`)

# # round CBD values to 3 decimal places and then take 'unique()' again
# CBD.vals$CanopyBulkDensity <- round(CBD.vals$CanopyBulkDensity, digits = 3)
# 
# CBD.vals <- unique(CBD.vals$CanopyBulkDensity)
# 
# CBD.vals <- as.data.frame(CBD.vals)
# 
# CBD.vals <- CBD.vals %>%
#   rename(CanopyBulkDensity = CBD.vals)


# use these values in a new predictions dataframe
newdat.CBDraster <- expand.grid(CanopyBulkDensity= CBD.vals$CanopyBulkDensity,
                                ID = "dummy", weight= 5000)

# design matrix (fixed effects)
mm.CBDraster <- model.matrix(delete.response(terms(CBD.rsf)),newdat.CBDraster)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CBDraster$STATUS <- plogis(drop(mm.CBDraster %*% fixef(CBD.rsf)[["cond"]]))

# create a reclassification matrix from the predicted 'STATUS' values 
rclmat.CBD <- matrix(c(CBD.vals$CanopyBulkDensity, newdat.CBDraster$STATUS),
                     ncol= 2)

# reclassify raster
CBD.probability <- classify(CBD.round, rclmat.CBD)

# save probability raster #
writeRaster(CBD.probability, "data/probability_rasters/CBD_probability.tif",
            overwrite = TRUE)



# Canopy layer count probability raster ----
# pull unique values from raster (omit NA's)
CLC.vals <- unique(CLC.scaled)

# rename 
CLC.vals <- CLC.vals %>%
  rename(CanopyLayerCount = `CFO-California-CanopyLayerCount-Summer2020-00010m`)

# use these values in a new predictions dataframe
newdat.CLCraster <- expand.grid(CanopyLayerCount= CLC.vals$CanopyLayerCount,
                                ID = "dummy", weight= 5000)

# design matrix (fixed effects)
mm.CLCraster <- model.matrix(delete.response(terms(CLC.rsf)),newdat.CLCraster)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.CLCraster$STATUS <- plogis(drop(mm.CLCraster %*% fixef(CLC.rsf)[["cond"]]))

# create a reclassification matrix from the predicted 'STATUS' values 
rclmat.CLC <- matrix(c(CLC.vals$CanopyLayerCount, newdat.CLCraster$STATUS),
                     ncol= 2)

# reclassify raster
CLC.probability <- classify(CLC.scaled, rclmat.CLC)

# save probability raster #
writeRaster(CLC.probability, "data/probability_rasters/CLC_probability.tif",
            overwrite = TRUE)



# Ladder fuel density probability raster ----
# pull unique values from raster (omit NA's)
LFD.vals <- unique(LFD.scaled)

# rename 
LFD.vals <- LFD.vals %>%
  rename(LadderFuelDensity = `CFO-California-LadderFuelDensity-Summer2020-00010m`)

# use these values in a new predictions dataframe
newdat.LFDraster <- expand.grid(LadderFuelDensity= LFD.vals$LadderFuelDensity,
                                ID = "dummy", weight= 5000)

# design matrix (fixed effects)
mm.LFDraster <- model.matrix(delete.response(terms(LFD.rsf)),newdat.LFDraster)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.LFDraster$STATUS <- plogis(drop(mm.LFDraster %*% fixef(LFD.rsf)[["cond"]]))

# create a reclassification matrix from the predicted 'STATUS' values 
rclmat.LFD <- matrix(c(LFD.vals$LadderFuelDensity, newdat.LFDraster$STATUS),
                     ncol= 2)

# reclassify raster
LFD.probability <- classify(LFD.scaled, rclmat.LFD)

# save probability raster #
writeRaster(LFD.probability, "data/probability_rasters/LFD_probability.tif",
            overwrite = TRUE)



# Surface fuels probability raster ----
# pull unique values from raster (omit NA's)
SF.vals <- unique(SF.scaled)

# rename 
SF.vals <- SF.vals %>%
  rename(SurfaceFuels = `CFO-California-SurfaceFuels-Summer2020-00010m`)

# use these values in a new predictions dataframe
newdat.SFraster <- expand.grid(SurfaceFuels= SF.vals$SurfaceFuels,
                               ID = "dummy", weight= 5000)

# design matrix (fixed effects)
mm.SFraster <- model.matrix(delete.response(terms(SF.rsf)),newdat.SFraster)

# linear predictor (for GLMMs, backtransform w/ the inverse link function - for 
#         binomial - plogis())
newdat.SFraster$STATUS <- plogis(drop(mm.SFraster %*% fixef(SF.rsf)[["cond"]]))

# create a reclassification matrix from the predicted 'STATUS' values 
rclmat.SF <- matrix(c(SF.vals$SurfaceFuels, newdat.SFraster$STATUS),
                    ncol= 2)

# reclassify raster
SF.probability <- classify(SF.scaled, rclmat.SF)

# save probability raster #
writeRaster(SF.probability, "data/probability_rasters/SF_probability.tif",
            overwrite = TRUE)



# ---
# 5. Create composite habitat probability raster ----
# ---
# Load in individual probability rasters (if needed)
CBH <- rast("data/probability_rasters/CBH_probability.tif")
CBD <- rast("data/probability_rasters/CBD_probability.tif")
CLC <- rast("data/probability_rasters/CLC_probability.tif")
LFD <- rast("data/probability_rasters/LFD_probability.tif")
SF <- rast("data/probability_rasters/SF_probability.tif")

# multiply all probability rasters together
composite <- CBH * CBD * CLC * LFD * SF

# save #
writeRaster(composite, "data/probability_rasters/composite_probability.tif", 
            overwrite = TRUE)


# ---
# 6. Determine threshold for probability of use (good vs. bad habitat) ----
# ---

################################################################################
## This code identifies a threshold in CSO habitat models that can best separate
# "Good" from "Bad" habitat

## Code by Kyle Rodman (Ecological Restoration Institute) 
# 4/12/2023

################################################################################
### Bring in necessary packages
package.list <- c("here", "tidyverse", "sf", "terra", "pROC")

## Installing them if they aren't already on the computer
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

################################################################################
### Read in files

## Owl used and available points (from RSF)
night_rsf <- readRDS("data/xSF_data/night_rsf")


# ## Owl - pseudo-absence points
# absPts <- st_read(here("Data", "Locations", "MSO_availLocs.shp"))
# 
# ## Owl - occurrence points
# presPts <- read.csv(here("Data", "Locations", "masterMSO_20200807_geeRevised_v3.csv"))

## Owl probability of use 
prob <- rast("data/probability_rasters/composite_probability.tif")

################################################################################
### Extract habitat values

# create sf object for all points 
night_rsf.sf <- st_as_sf(night_rsf, coords = c("x_", "y_"), crs = crs(prob))

# create a 0/1 column for used and available 
night_rsf.sf <- night_rsf.sf %>%
  mutate(case = ifelse(case_ == "FALSE", 0, 1))

## Convert to SpatVector and extract sdm values
night <- vect(night_rsf.sf)
owlSDMValues <- terra::extract(prob, night)

## And convert back to DF for pROC calculations
night <- st_as_sf(night) %>%
  st_drop_geometry() %>%
  mutate(sdmVal = owlSDMValues$`CFO-California-CanopyBaseHeight-Summer2020-00010m`,
         presence = factor(case))
boxplot(night$sdmVal~night$case)

################################################################################
### pROC calculations to get data threshold

## Create ROC for observed (presence/absence) and predicted (habitat SDM)
owlROC <- roc(data = night, response = case, predictor = sdmVal)
plot(owlROC)

## Sum of sensitivity and specificity - the true skill statistic
sensSpec <- owlROC$sensitivities + owlROC$specificities

## Getting the threshold that corresponds to the max TSS
owlROC$thresholds[which.max(sensSpec)] # 0.0326 FOR UNSCALED DATA, 0.399396911 FOR SCALED DATA

