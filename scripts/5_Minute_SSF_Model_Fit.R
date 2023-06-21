# Script name: 5_Minute_SSF_Model_Fit.R
# Author: M. E. Wright
# Date created: December 2, 2022
#

# Notes: 
#      This script does the following fits models for minute-to-minute step
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
min  <- readRDS("data/xSF_data/minute_ssf")


# ---
# 2. Center and scale variables, create weight column ----
# ---

# center and scale all variables
min$CanopyBaseHeight <- scale(min$CanopyBaseHeight, center = TRUE, scale = TRUE)
min$CanopyBulkDensity <- scale(min$CanopyBulkDensity, center = TRUE, scale = TRUE)
min$CanopyLayerCount <- scale(min$CanopyLayerCount, center = TRUE, scale = TRUE)
min$LadderFuelDensity <- scale(min$LadderFuelDensity, center = TRUE, scale = TRUE)
min$SurfaceFuels <- scale(min$SurfaceFuels, center = TRUE, scale = TRUE)

# create unique ID column that combines territory and year 
min$id  <- paste(min$Territory, min$Year, min$usfws_id, sep = "_")

# create a weight column 
min$weight <- ifelse(min$case_ == "TRUE", 1, 5000)

str(min)

# try to match the structure of Muff et al. 2018 

# convert 'id' to integer column 'ID'
min$ID <- as.integer(factor(min$id, 
                             levels=unique(min$id)))

# convert logical column 'case_' to an integer column of 1's and 0's 'STATUS'
min$STATUS <- as.integer(ifelse(min$case_ == "TRUE", 1, 0))


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
                    family = poisson, data = min, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.cbh$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.cbh$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
MINssf.cbh <- glmmTMB::fitTMB(TEMP.cbh)
summary(MINssf.cbh)

# 95% CI for fixed and random effects 
confint(MINssf.cbh)


## SAVE MODEL OBJECT ##
saveRDS(MINssf.cbh, file = "data/model_objects/minute_ssf/CanopyBaseHeight")

################################################################################


# HOUR SSF: CANOPY BULK DENSITY ----
TEMP.cbd <- glmmTMB(STATUS ~  -1 + CanopyBulkDensity + (1 | step_id_) +
                      (0 + CanopyBulkDensity | ID), 
                    family = poisson, data = min, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.cbd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.cbd$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
MINssf.cbd <- glmmTMB::fitTMB(TEMP.cbd)
summary(MINssf.cbd)

# 95% CI for fixed and random effects 
confint(MINssf.cbd)


## SAVE MODEL OBJECT ##
saveRDS(MINssf.cbd, file = "data/model_objects/minute_ssf/CanopyBulkDensity")


################################################################################


# HOUR SSF: CANOPY LAYER COUNT ----
TEMP.clc <- glmmTMB(STATUS ~  -1 + CanopyLayerCount + (1 | step_id_) +
                      (0 + CanopyLayerCount | ID), 
                    family = poisson, data = min, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.clc$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.clc$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
MINssf.clc <- glmmTMB::fitTMB(TEMP.clc)
summary(MINssf.clc)

# 95% CI for fixed and random effects 
confint(MINssf.clc)


## SAVE MODEL OBJECT ##
saveRDS(MINssf.clc, file = "data/model_objects/minute_ssf/CanopyLayerCount")


################################################################################


# HOUR SSF: LADDER FUEL DENSITY ----
TEMP.lfd <- glmmTMB(STATUS ~  -1 + LadderFuelDensity + (1 | step_id_) +
                      (0 + LadderFuelDensity | ID), 
                    family = poisson, data = min, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.lfd$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.lfd$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
MINssf.lfd <- glmmTMB::fitTMB(TEMP.lfd)
summary(MINssf.lfd)

# 95% CI for fixed and random effects 
confint(MINssf.lfd)


## SAVE MODEL OBJECT ##
saveRDS(MINssf.lfd, file = "data/model_objects/minute_ssf/LadderFuelDensity")


################################################################################


# HOUR SSF: CANOPY BULK DENSITY ----
TEMP.sf <- glmmTMB(STATUS ~  -1 + SurfaceFuels + (1 | step_id_) +
                     (0 + SurfaceFuels | ID), 
                   family = poisson, data = min, doFit = FALSE)


# set the value of standard deviation for the first random effect (1 | step_id_)
TEMP.sf$parameters$theta[1] = log(1e3)

# tell glmmTMB not to change the first standard deviation, 
# all other values are freely estimated (and are different from each other)
TEMP.sf$mapArg = list(theta = factor(c(NA, 1)))

# fit the model 
MINssf.sf <- glmmTMB::fitTMB(TEMP.sf)
summary(MINssf.sf)

# 95% CI for fixed and random effects 
confint(MINssf.sf)


## SAVE MODEL OBJECT ##
saveRDS(MINssf.sf, file = "data/model_objects/minute_ssf/SurfaceFuels")

