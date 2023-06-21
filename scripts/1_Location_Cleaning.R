# Script name: 1_Location_Cleaning.R
# Author: Jessalyn Ayars
# Date created: October 20, 2021
#

# Notes: This code script is the same as the cleaning code script for Jessalyn's 
# GPS-fitness project with some added annotation and minor changes as needed. 
# Script format the raw GPS telemetry data for proper use with HSF and SSF. 
# ---     


# ---
# 1. Loading packages and data ----
# ---
# Load packages
library(tidyverse)
library(here)
library(readxl)
library(lubridate)
library(sf)
library(tibbletime)
library(amt)

# Load in GPS data
# load csv with SPI data
test_data_SPI <- read_csv(here('data', 'gps', "sierra_gps_SPI.csv"))

# drop ID column 
# rename columns
# select only columns needed
# add a column ('data.origin') labeling these points as 'SPI' data
SPI <- test_data_SPI %>%
  dplyr::select(-'id...1') %>%
  rename("Territory" = 'short_territory_name', 'tag_number' = "tag_number...21", 'telemetry_model' = "telemetry_model...23", 'gps_record_id' = "id...24", 'num_sats' = 'satellite', ) %>%
  dplyr::select(Territory, pac_id, usfws_id, sex, date_gmt, time_gmt, tag_number, attach_type, telemetry_model, gps_record_id, latitude, longitude, fix, altitude, duration, dop, num_sats, speed_knots, search_time, voltage, speed_km_hr, azimuth_dtn, pdop, hdop, vdop, sats_in_view, snr_avg_db, delta_s, e_res, eldorado_demography_type, is_king_fire_gavin, is_king_fire_southern, is_sierra_nevada_tagging, is_sierra_demography, is_lassen_demography, is_yosemite, is_sequoia_kings, county, land_owner, forest, district) %>%
  mutate(data.origin = "SPI")

# load csv with other CSO telemetry data
# add a data origing column
test_data_other <- read_csv(here('data', 'gps', "sierra_gps_noSPI.csv")) %>%
  mutate(data.origin = "not.SPI")

# bind together SPI and non-SPI data
test_data <- rbind(SPI, test_data_other) 

# load in telemetry unit model data
tmfc <- read_xlsx(here('data', 'gps', "telemetry_model_field_corrected.xlsx"))

# add corrected telemetry model data to gps records
data1 <- left_join(test_data, tmfc, by = "gps_record_id") %>%
  mutate(telemetry_model = ifelse(is.na(telemetry_model.y), telemetry_model.x, telemetry_model.y)) %>%
  dplyr::select(-telemetry_model.x, -telemetry_model.y)




# ---
# 2. Cleaning location data ----
# --- 
# Properly formats data, getting rid of fixes without recorded GPS locations, 
# adding in a datetime variable in proper format for use with 'amt' package, 
# releveling variables as needed, and dropping fixes of poor quality (i.e. low
# voltage, not enough satellites, ect.). Duplicate records are also dropped from 
# the dataset. 

# drops all records where there is no recorded gps location
lon_0 <- as.integer(which(data1$longitude == 0))
data1 <- data1[-lon_0,]


data <- data1 %>%
  # parsing dates/times into lubridate format
  mutate(time_gmt = as.character(time_gmt)) %>%
  mutate(datetime_gmt = mdy_hms(str_c(date_gmt, " ", time_gmt), tz = "GMT")) %>%
  # switching to pacific time ### someday you need to switch this out of daylight savings ###
  mutate(datetime = with_tz(datetime_gmt, "US/Pacific")) %>%
  # switching to factors
  mutate(across(where(is.character), as.factor)) %>%
  # make usfws id a factor
  mutate(usfws_id = as.factor(usfws_id)) %>%
  # time of day w/o date (in progress)
  mutate(hour_daily = hour(datetime)) %>%
  
  # month to see if daylight savings will be an issue (nope)
  mutate(month_obs = as.factor(month(datetime))) %>%
  
  # getting the number of fixes by bird and type of telemetry model (assuming one device per bird per season/year)
  group_by(usfws_id, telemetry_model) %>%
  mutate(n_fixes = n()) %>%
  
  # releveling telemetry model with redundant names
  mutate(telemetry_model = recode_factor(telemetry_model, 
                                         "SWIFT PP 120" = "Swift PP 120")) %>%
  
  # relevel territory MIDMD to PLCWA
  mutate(Territory = recode_factor(Territory, "MIDMD" = "PLCWA")) %>%
  
  # filter out bad locations based on DOP < 5 /voltage > 3/6 /numsats >= 3. removed ~3,000 obs
  dplyr::filter(is.na(num_sats)|num_sats >= 3) %>%
  dplyr::filter(is.na(dop)|dop < 5) %>%
  dplyr::filter(is.na(hdop)|hdop < 5) %>%
  dplyr::filter(is.na(vdop)|vdop < 5) %>%
  dplyr::filter(is.na(voltage)|voltage > 3.6) %>%
  
  # make an sf object for plotting purposes
  drop_na(latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(4326), remove = FALSE) %>%
  st_transform(32610) %>% 
  # getting utm columns so that track can use it
  mutate(long_utm = unlist(map(geometry,1)),
         lat_utm = unlist(map(geometry,2))) %>%
  dplyr::select(-geometry) 

# determine where there are duplicated records and drop these
duplicated <- data %>%
  filter(usfws_id == "195722595" | usfws_id == "138784870" | usfws_id == "195738005") %>%
  distinct(datetime, .keep_all = TRUE) 

# keep unduplicated records 
unduplicated <- data %>%
  filter(usfws_id != "195722595" & usfws_id != "138784870" & usfws_id != "195738005") %>%
  rbind(duplicated)



# --- 
# 3. Filtering out daylight hours ----
# ---
# Consider only the hours of peak activity for CSO 

# create a vector of unique study days
study_days <- unique(date(unduplicated$datetime))

# create a 'to_remove' dataframe that includes all daylight hour fixes between
# 0700 - 1900 hours# doing it just by filtering from 7-19:00
gps_data_simplified_post <- unduplicated %>%
  arrange(datetime) %>%
  as_tbl_time(index = datetime)
to_remove <- slice(gps_data_simplified_post, 0)

for (study_day in study_days) {
  study_day = as.Date(study_day, origin = "1970-01-01")
  morning <- ymd_hms(str_c(as.character(study_day), "07:00:00"), tz = "US/Pacific")
  evening <- ymd_hms(str_c(as.character(study_day), "19:00:00"), tz = "US/Pacific")
  to_remove.temp <- filter_time(gps_data_simplified_post, morning ~ evening)
  to_remove <- rbind(to_remove, to_remove.temp)
}

gps_data_simplified_post <- anti_join(gps_data_simplified_post, to_remove,
                                      by = "gps_record_id")




# ---
# 4. Making an effort variable based on lag time ----
# ---

# smaller lag --> more effort
# 1/(average lag time per owl per year?)

gps_data_effort <- gps_data_simplified_post %>%
  mutate(Year = year(datetime)) %>%
  arrange(datetime) %>%
  group_by(Year, usfws_id) %>%
  mutate(lag1 = as.numeric(difftime(datetime, lag(datetime), units = "mins"))) %>% # some are NAs because they are the beginning
  mutate(effort = 1/(mean(lag1, na.rm = TRUE))) 



# ---
# 5. Assigning owls to each analysis (HSF or SSF) ----
# ---
# Uses sampling rates to determine which owls will be in which analysis group 

test.track <- gps_data_effort %>%
  ungroup() %>%
  make_track(long_utm, lat_utm, datetime, Territory, usfws_id, Year)

# summarizing sampling rates to determine how i should divvy this up
test.rates <- test.track %>%
  summarize_sampling_rate_many(c("usfws_id", "Year"), time_unit = "min")

test.rates.tojoin <- test.rates %>%
  dplyr::select(usfws_id, median, Year)

gps_data.xsf <- gps_data_effort %>%
  left_join(test.rates.tojoin, by = c("usfws_id", "Year")) %>%
  rename(median.sampling.rate = median) %>%
  mutate(which.xsf = NA)

rsf <- which(gps_data.xsf$median.sampling.rate > 500)
ssf.hr <- which(gps_data.xsf$median.sampling.rate < 500 & gps_data.xsf$median.sampling.rate > 30)
ssf.min <- which(gps_data.xsf$median.sampling.rate < 30)

gps_data.xsf[rsf,]$which.xsf <- "rsf"
gps_data.xsf[ssf.hr,]$which.xsf <- "ssf.hr"
gps_data.xsf[ssf.min,]$which.xsf <- "ssf.min"



# --- 
# 6. Removing potential issue locations from 2015 ----
# ---
# CFO data only goes as far back as 2016
sl <- st_read(here("data", "shapefiles", "20190411_Anu_digitized_salvage_logging_2016NAIP_erase_GNNopen_unburned.shp")) %>%
  st_transform(32610) %>%
  st_union()
territories <- st_read(here("Data", "shapefiles", "territory_polylgons.shp")) %>%
  st_transform(32610) %>%
  st_union()

gps.2015.bad2 <- gps_data.xsf %>%
  dplyr::filter(Year == 2015) %>%
  st_as_sf(coords = c("long_utm", "lat_utm"), crs = 32610, remove = FALSE) %>%
  ungroup()

gps.2015.bad <- gps.2015.bad2 %>%
  mutate(in.sl = as.numeric(st_intersects(gps.2015.bad2, sl))) %>%
  drop_na(in.sl)

gps.data <- anti_join(gps_data.xsf, gps.2015.bad)



## SAVE ## 
saveRDS(gps.data, file = "data/gps/gps_cleaned")
