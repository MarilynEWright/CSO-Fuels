# Script name: 3_xSF_Model_Fit.R
# Author: M. E. Wright
# Date created: November 8, 2022
#

# Notes: 
#      This script does the following fits models for habitat selection and step
#       selection functions for fuels metrics
# ---

# ---
# 1. Loading packages and data ----
# ---
library(lubridate)
library(raster)
library(move)
library(amt) 
library(broom)
library(glmmTMB)
library(tidyverse)
library(dplyr)
library(broom)

# turn off scientific notation 
options(scipen=999)

# load in data 
night <- readRDS("data/xSF_data/night_rsf")

# ---
# 2. Center and scale variables, create weight column ----
# ---

# center and scale all variables
night$CanopyBaseHeight <- scale(night$CanopyBaseHeight, center = TRUE, scale = TRUE)
night$CanopyBulkDensity <- scale(night$CanopyBulkDensity, center = TRUE, scale = TRUE)
night$CanopyLayerCount <- scale(night$CanopyLayerCount, center = TRUE, scale = TRUE)
night$LadderFuelDensity <- scale(night$LadderFuelDensity, center = TRUE, scale = TRUE)
night$SurfaceFuels <- scale(night$SurfaceFuels, center = TRUE, scale = TRUE)

# create a weight column 
night$weight <- ifelse(night$case_ == "TRUE", 1, 5000)


# create unique ID column that combines territory and year 
night$id <- paste(night$Territory, night$Year, night$usfws_id, sep = "_")

# Due to high autocorrelation, these models will be fit as univariate models
# with a random slope for individuals 


# ---
# 3. Fit FIXED VARIANCE HSF models for nightly data for each variable ----
# ---

# Look at the proportion of used and available by animals 
with(night, prop.table(table(id, case_), 1))

# univariate models with random slopes and intercepts with fixed intercept 
# variance
# (Muff et al. (2018) suggest decision to fix or estimate the intercept for HSFs
# is not critical like it is for SSFs, but it seemed like a safe option since 
# we have some individuals with ~30 used points only)
# variance (WEIGHT: used points = 1, available points = 5000)
# fixed effects: each variable
# random effects: intercepts and slopes (variable)



# FIT HSF MODELS (FIXED INTERCEPT VARIANCE) ---- 
# 1) set up the models, but do not fit 
# 2) fix the standard deviation of the first random term (1 | ID) using σ=103 ,
#     which corresponds to a variance of 10^6
# 3) tell 'glmmTMB' not to change the entry of the vector of variances and give 
#     all other variances another indicator to make sure they can be be freely 
#     estimated
# 4) fit the models
# 5) save the results 

str(night)

# try to match the structure of Muff et al. 2018 

# convert 'id' to integer column 'ID'
night$ID <- as.integer(factor(night$id, 
                              levels=unique(night$id)))

# convert logical column 'case_' to an integer column of 1's and 0's 'STATUS'
night$STATUS <- as.integer(ifelse(night$case_ == "TRUE", 1, 0))



# HSF: canopy base height ----
# set up model
temp.cbh <- glmmTMB(STATUS ~ CanopyBaseHeight + (1 | ID) + 
                     (0 + CanopyBaseHeight | ID), data = night, 
                 family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.cbh$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.cbh$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.cbh <- glmmTMB::fitTMB(temp.cbh)
summary(hsf.cbh)

## SAVE MODEL OBJECT ##
saveRDS(hsf.cbh, file = "data/model_objects/night_hsf/fixed_variance/CanopyBaseHeight")

################################################################################

# HSF: canopy bulk density ----
# set up model
temp.cbd <- glmmTMB(STATUS ~ CanopyBulkDensity + (1 | ID) +
                     (0 + CanopyBulkDensity | ID), data = night, 
                 family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.cbd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.cbd$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.cbd <- glmmTMB::fitTMB(temp.cbd)
summary(hsf.cbd)

## SAVE MODEL OBJECT ##
saveRDS(hsf.cbd, file = "data/model_objects/night_hsf/fixed_variance/CanopyBulkDensity")

################################################################################


# HSF: canopy layer count ----
# set up model
temp.clc <- glmmTMB(STATUS ~ CanopyLayerCount + (1 | ID) +
                     (0 + CanopyLayerCount | ID), data = night,
                 family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.clc$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.clc$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.clc <- glmmTMB::fitTMB(temp.clc)
summary(hsf.clc)

## SAVE MODEL OBJECT ##
saveRDS(hsf.clc, file = "data/model_objects/night_hsf/fixed_variance/CanopyLayerCount")


################################################################################

# HSF: ladder fuel density ----
# set up model
temp.lfd <- glmmTMB(STATUS ~ LadderFuelDensity + (1 | ID) +
                     (0 + LadderFuelDensity | ID), data = night,
                 family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.lfd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.lfd$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.lfd <- glmmTMB::fitTMB(temp.lfd)
summary(hsf.lfd)

## SAVE MODEL OBJECT ##
saveRDS(hsf.lfd, file = "data/model_objects/night_hsf/fixed_variance/LadderFuelDensity")

################################################################################

# HSF: surface fuels ----
# set up model
temp.sf <- glmmTMB(STATUS ~ SurfaceFuels + (1 | ID) +
                    (0 + SurfaceFuels | ID), data = night, 
                family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.sf$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.sf$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.sf <- glmmTMB::fitTMB(temp.sf)
summary(hsf.sf)

## SAVE MODEL OBJECT ##
saveRDS(hsf.sf, file = "data/model_objects/night_hsf/fixed_variance/SurfaceFuels")

################################################################################



# FIT HSF MODELS (FIXED INTERCEPT VARIANCE), ORIGINAL DATA ONLY  ---- 
# create a filtered data set with only the original nightly telemetry data
night.OG <- night %>%
  filter(data_source == "night")

# HSF: canopy base height ----
# set up model
temp.cbhRED <- glmmTMB(STATUS ~ CanopyBaseHeight + (1 | ID) + 
                      (0 + CanopyBaseHeight | ID), data = night.OG, 
                    family = binomial(), doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.cbhRED$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.cbhRED$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.cbhRED <- glmmTMB::fitTMB(temp.cbhRED)
summary(hsf.cbhRED)

## SAVE MODEL OBJECT ##
saveRDS(hsf.cbhRED, file = "data/model_objects/night_hsf/fixed_variance/original_data_only/CanopyBaseHeight")

################################################################################

# HSF: canopy bulk density ----
# set up model
temp.cbdRED <- glmmTMB(STATUS ~ CanopyBulkDensity + (1 | ID) +
                      (0 + CanopyBulkDensity | ID), data = night.OG, 
                    family = binomial, doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.cbdRED$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.cbdRED$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.cbdRED <- glmmTMB::fitTMB(temp.cbdRED)
summary(hsf.cbdRED)

## SAVE MODEL OBJECT ##
saveRDS(hsf.cbdRED, file = "data/model_objects/night_hsf/fixed_variance/original_data_only/CanopyBulkDensity")

################################################################################


# HSF: canopy layer count ----
# set up model
temp.clcRED <- glmmTMB(STATUS ~ CanopyLayerCount + (1 | ID) +
                      (0 + CanopyLayerCount | ID), data = night.OG,
                    family = binomial, doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.clcRED$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.clcRED$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.clcRED <- glmmTMB::fitTMB(temp.clcRED)
summary(hsf.clcRED)

## SAVE MODEL OBJECT ##
saveRDS(hsf.clcRED, file = "data/model_objects/night_hsf/fixed_variance/original_data_only/CanopyLayerCount")


################################################################################

# HSF: ladder fuel density ----
# set up model
temp.lfdRED <- glmmTMB(STATUS ~ LadderFuelDensity + (1 | ID) +
                      (0 + LadderFuelDensity | ID), data = night.OG,
                    family = binomial, doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.lfdRED$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.lfdRED$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.lfdRED <- glmmTMB::fitTMB(temp.lfdRED)
summary(hsf.lfdRED)

## SAVE MODEL OBJECT ##
saveRDS(hsf.lfdRED, file = "data/model_objects/night_hsf/fixed_variance/original_data_only/LadderFuelDensity")

################################################################################

# HSF: surface fuels ----
# set up model
temp.sfRED <- glmmTMB(STATUS ~ SurfaceFuels + (1 | ID) +
                     (0 + SurfaceFuels | ID), data = night.OG, 
                   family = binomial, doFit = F, weights = weight)

# fix standard deviation of 1st random term 
temp.sfRED$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first entry of the vector of variances, and 
# give all other variances another indicator to make sure they can be freely 
# estimated
temp.sfRED$mapArg = list(theta=factor(c(NA, 1)))

# fit the model 
hsf.sfRED <- glmmTMB::fitTMB(temp.sfRED)
summary(hsf.sfRED)

## SAVE MODEL OBJECT ##
saveRDS(hsf.sfRED, file = "data/model_objects/night_hsf/fixed_variance/original_data_only/SurfaceFuels")

################################################################################

# ---
# 4. Fit ESTIMATED VARIANCE HSF models for nightly data for each variable ----
# ---

# FIT HSF MODELS (INTERCEPT VARIANCE ESTIMATED) ---- 
# HSF: canopy base height ----
hsf.cbh <- glmmTMB(STATUS ~ CanopyBaseHeight + (1 | ID) + 
                      (0 + CanopyBaseHeight | ID), data = night, 
                    family = binomial(), weights = weight)


# HSF: canopy bulk density ----
hsf.cbd <- glmmTMB(STATUS ~ CanopyBulkDensity + (1 | ID) +
                      (0 + CanopyBulkDensity | ID), data = night, 
                    family = binomial, weights = weight)


# HSF: canopy layer count ----
hsf.clc <- glmmTMB(STATUS ~ CanopyLayerCount + (1 | ID) +
                      (0 + CanopyLayerCount | ID), data = night,
                    family = binomial, weights = weight)


# HSF: ladder fuel density ----
hsf.lfd <- glmmTMB(STATUS ~ LadderFuelDensity + (1 | ID) +
                      (0 + LadderFuelDensity | ID), data = night,
                    family = binomial, weights = weight)


# HSF: surface fuels ----
hsf.sf <- glmmTMB(STATUS ~ SurfaceFuels + (1 | ID) +
                     (0 + SurfaceFuels | ID), data = night, 
                   family = binomial, weights = weight)

## SAVE MODEL OBJECTS ## 
saveRDS(hsf.cbh, file = "data/model_objects/night_hsf/estimated_variance/CanopyBaseHeight")
saveRDS(hsf.cbd, file = "data/model_objects/night_hsf/estimated_variance/CanopyBulkDensity")
saveRDS(hsf.clc, file = "data/model_objects/night_hsf/estimated_variance/CanopyLayerCount")
saveRDS(hsf.lfd, file = "data/model_objects/night_hsf/estimated_variance/LadderFuelDensity")
saveRDS(hsf.sf, file = "data/model_objects/night_hsf/estimated_variance/SurfaceFuels")

################################################################################

# FIT HSF MODELS (INTERCEPT VARIANCE ESTIMATED), ORIGINAL DATA ONLY  ---- 
# use the filtered data set with only the original nightly telemetry data
# 'night.OG'

# HSF: canopy base height ----
hsf.cbhRED <- glmmTMB(STATUS ~ CanopyBaseHeight + (1 | ID) + 
                     (0 + CanopyBaseHeight | ID), data = night.OG, 
                   family = binomial(), weights = weight)


# HSF: canopy bulk density ----
hsf.cbdRED <- glmmTMB(STATUS ~ CanopyBulkDensity + (1 | ID) +
                     (0 + CanopyBulkDensity | ID), data = night.OG, 
                   family = binomial, weights = weight)


# HSF: canopy layer count ----
hsf.clcRED <- glmmTMB(STATUS ~ CanopyLayerCount + (1 | ID) +
                     (0 + CanopyLayerCount | ID), data = night.OG,
                   family = binomial, weights = weight)


# HSF: ladder fuel density ----
hsf.lfdRED <- glmmTMB(STATUS ~ LadderFuelDensity + (1 | ID) +
                     (0 + LadderFuelDensity | ID), data = night.OG,
                   family = binomial, weights = weight)


# HSF: surface fuels ----
hsf.sfRED <- glmmTMB(STATUS ~ SurfaceFuels + (1 | ID) +
                    (0 + SurfaceFuels | ID), data = night.OG, 
                  family = binomial, weights = weight)

## SAVE MODEL OBJECTS ## 
saveRDS(hsf.cbhRED, file = "data/model_objects/night_hsf/estimated_variance/original_data_only/CanopyBaseHeight")
saveRDS(hsf.cbdRED, file = "data/model_objects/night_hsf/estimated_variance/original_data_only/CanopyBulkDensity")
saveRDS(hsf.clcRED, file = "data/model_objects/night_hsf/estimated_variance/original_data_only/CanopyLayerCount")
saveRDS(hsf.lfdRED, file = "data/model_objects/night_hsf/estimated_variance/original_data_only/LadderFuelDensity")
saveRDS(hsf.sfRED, file = "data/model_objects/night_hsf/estimated_variance/original_data_only/SurfaceFuels")


