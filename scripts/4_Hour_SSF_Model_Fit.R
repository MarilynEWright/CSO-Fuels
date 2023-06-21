# Script name: 4_Hour_SSF_Model_Fit.R
# Author: M. E. Wright
# Date created: December 2, 2022
#

# Notes: 
#      This script does the following fits models for hourly step
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
hour  <- readRDS("data/xSF_data/hour_ssf")


# ---
# 2. Center and scale variables, create weight column ----
# ---

# center and scale all variables
hour$CanopyBaseHeight <- scale(hour$CanopyBaseHeight, center = TRUE, scale = TRUE)
hour$CanopyBulkDensity <- scale(hour$CanopyBulkDensity, center = TRUE, scale = TRUE)
hour$CanopyLayerCount <- scale(hour$CanopyLayerCount, center = TRUE, scale = TRUE)
hour$LadderFuelDensity <- scale(hour$LadderFuelDensity, center = TRUE, scale = TRUE)
hour$SurfaceFuels <- scale(hour$SurfaceFuels, center = TRUE, scale = TRUE)

# create unique ID column that combines territory and year 
hour$id  <- paste(hour$Territory, hour$Year, hour$usfws_id, sep = "_")

# create a weight column 
hour$weight <- ifelse(hour$case_ == "TRUE", 1, 5000)

str(hour)

# try to match the structure of Muff et al. 2018 

# convert 'id' to integer column 'ID'
hour$ID <- as.integer(factor(hour$id, 
                              levels=unique(hour$id)))

# convert logical column 'case_' to an integer column of 1's and 0's 'STATUS'
hour$STATUS <- as.integer(ifelse(hour$case_ == "TRUE", 1, 0))


# ---
# 4. Fit SSF models for hourly data for each variable ----
# ---
# For the step selection functions, we will fit again fit univariate models with
# random intercepts, random slopes, and including movement parameter 
# (step length). The models will also include step_id as the stratum and fixed 
# variance

# FIT HOURLY SSF MODELS (FIXED INTERCEPT VARIANCE) ---- 
# 1) set up the models, but do not fit 
# 2) fix the standard deviation of the first random term (1 | step_id_) using σ=103 ,
#     which corresponds to a variance of 10^6
# 3) tell 'glmmTMB' not to change the entry of the vector of variances and give 
#     all other variances another indicator to make sure they can be be freely 
#     estimated
# 4) fit the models
# 5) save the results 

# HOUR SSF: CANOPY BASE HEIGHT ----
# set up the model 
TEMP.cbh <- glmmTMB(STATUS ~  -1 + CanopyBaseHeight + (1 | step_id_) +
                            (0 + CanopyBaseHeight | ID), 
                          family = poisson, data = hour, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.cbh$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.cbh$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
HRssf.cbh <- glmmTMB::fitTMB(TEMP.cbh)
summary(HRssf.cbh)

# 95% CI for fixed and random effects 
confint(HRssf.cbh)


## SAVE MODEL OBJECT ##
saveRDS(HRssf.cbh, file = "data/model_objects/hourly_ssf/CanopyBaseHeight")

################################################################################


# HOUR SSF: CANOPY BULK DENSITY ----
TEMP.cbd <- glmmTMB(STATUS ~  -1 + CanopyBulkDensity + (1 | step_id_) +
                      (0 + CanopyBulkDensity | ID), 
                    family = poisson, data = hour, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.cbd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.cbd$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
HRssf.cbd <- glmmTMB::fitTMB(TEMP.cbd)
summary(HRssf.cbd)

# 95% CI for fixed and random effects 
confint(HRssf.cbd)


## SAVE MODEL OBJECT ##
saveRDS(HRssf.cbd, file = "data/model_objects/hourly_ssf/CanopyBulkDensity")


################################################################################


# HOUR SSF: CANOPY LAYER COUNT ----
TEMP.clc <- glmmTMB(STATUS ~  -1 + CanopyLayerCount + (1 | step_id_) +
                      (0 + CanopyLayerCount | ID), 
                    family = poisson, data = hour, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.clc$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.clc$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
HRssf.clc <- glmmTMB::fitTMB(TEMP.clc)
summary(HRssf.clc)

# 95% CI for fixed and random effects 
confint(HRssf.clc)


## SAVE MODEL OBJECT ##
saveRDS(HRssf.clc, file = "data/model_objects/hourly_ssf/CanopyLayerCount")


################################################################################


# HOUR SSF: LADDER FUEL DENSITY ----
TEMP.lfd <- glmmTMB(STATUS ~  -1 + LadderFuelDensity + (1 | step_id_) +
                      (0 + LadderFuelDensity | ID), 
                    family = poisson, data = hour, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.lfd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.lfd$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
HRssf.lfd <- glmmTMB::fitTMB(TEMP.lfd)
summary(HRssf.lfd)

# 95% CI for fixed and random effects 
confint(HRssf.lfd)


## SAVE MODEL OBJECT ##
saveRDS(HRssf.lfd, file = "data/model_objects/hourly_ssf/LadderFuelDensity")


################################################################################


# HOUR SSF: CANOPY BULK DENSITY ----
TEMP.sf <- glmmTMB(STATUS ~  -1 + SurfaceFuels + (1 | step_id_) +
                      (0 + SurfaceFuels | ID), 
                    family = poisson, data = hour, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.sf$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.sf$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
HRssf.sf <- glmmTMB::fitTMB(TEMP.sf)
summary(HRssf.sf)

# 95% CI for fixed and random effects 
confint(HRssf.sf)


## SAVE MODEL OBJECT ##
saveRDS(HRssf.sf, file = "data/model_objects/hourly_ssf/SurfaceFuels")

