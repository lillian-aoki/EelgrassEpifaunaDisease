# 01-DataPrep
## Last updated: 2024-03-19 by Lillian Aoki

## This script aggregates data for analysis of grazing impacts on eelgrass wasting disease
## Input datasets are derived elsewhere (i.e. project data)
## Output datasets are derived from this script and used in subsequent scripts
## Data dictionary is available for the output datasets

library(tidyverse)

# outputs:
# epiphyte_SEM_data_all_large.csv [used in TaxaInteractionSEMwithGrazing, GrazingDirectionSEM, partial_effect_plots, grazing_fig]
# epifauna_for_region_specific_models_no_epiphyte.csv [used in TaxaComparisonAllData]
# epifauna_site_for_plotting.csv [used in latitude_figures, maps]

# data inputs ####
# site metadata (for latitude)
meta <- read_csv("data/input/combined_site_metadata.csv")

# temperature data (meadow-scale, monthly mean and anomaly)
mur <- read_csv("data/input/monthly_temp_anomaly_mur_9y_sites.csv")
ghr <- read_csv("data/input/monthly_temp_anomaly_g1sst_9y_sites.csv")
monthly <- read_csv("data/input/monthly_mean_temps_mur_g1sst_9y_sites.csv")

# seagrass data (transect level)
sg <- read_csv("data/input/combined_transect_survey_metrics.csv")

# epifauna data (transect-level, large with BB included)
epi <- read_csv("data/input/EGWD_transect_data_v20230307_big_epi_with_BB_transect.csv")

# grazing scars for individual leaves
gz <- read_csv("data/input/grazing_scars_compiled.csv")

# leaf-level data (disease, leaf area, epiphyte load)
dis <- read_csv("data/input/meter_level_shoot_metrics_with_disease.csv")

# prep main dataset for modeling ####
# keep only June SST data
mur <- select(mur, c(Year, Meadow, Region, Site, TempAnomWarm_June))
ghr <- select(ghr, c(Year, Meadow, Region, Site, TempAnomWarm_June))
monthly$MonthW <- format(monthly$Month, "%B")
monthly <- subset(monthly, MonthW=="June")
monthly$Year <- as.numeric(format(monthly$Month, "%Y"))
temps <- rbind(mur, ghr)
temps <- left_join(temps, monthly)

# combine with a right-join to limit to sites with a temperature record
sg_temps <- right_join(sg, temps, by = c("Year", "Region", "SiteCode"="Site"))

# add epifauna - again limit to prevent any repeating rows. note, the epifauna are limited to temperature sites only
sg_temps$transect_unique_code <- paste(sg_temps$Region, ".", sg_temps$SiteCode, ".", 
                                       sg_temps$TidalHeight, sg_temps$Transect, ".", sg_temps$Year, sep="")
sg_temps_epi <- inner_join(sg_temps, epi)

# adjust grazing scars - only 5 blades in 2019 except at SD, WA analysis was post-hoc 
gz <- gz[-c(which(gz$Region=="SD" & gz$Year=="2019" & gz$Blade %in% c(2,3,5,6,7,9,10,11,13,14,15))),]
gz <- gz[-c(which(gz$Region=="WA" & gz$Year=="2019" & (gz$Blade<16 | gz$Blade>21))),]
gz <- subset(gz, !is.na(gz$GrazingScars))
# this keeps all the grazing for 2021 
gz_summ <- gz %>%
  group_by(Region, SiteCode, Year, Transect, TidalHeight) %>%
  summarise(GrazingScarsMeanTransect=mean(GrazingScars, na.rm=TRUE), CountScars=length(!is.na(GrazingScars)))
sg_temps_epi_gz <- left_join(sg_temps_epi, gz_summ)
# extract transect locations from the meta data
meta19 <- subset(meta, Year=="2019")
loc19 <- select(meta19, c("Region","SiteCode", "TidalHeight", "Transect", "TransectBeginDecimalLatitude",
                          "TransectBeginDecimalLongitude", "TransectEndDecimalLatitude", "TransectEndDecimalLongitude"))
# fill in missing transect location data (based on 2019 locations unless updated)
full_meta <- meta %>%
  rows_update(loc19, by = c("Region", "SiteCode", "TidalHeight", "Transect"))

meta_site <- full_meta %>%
  group_by(Year, Region, SiteCode) %>%
  summarise(Latitude=mean(as.numeric(TransectBeginDecimalLatitude)), Longitude=mean(TransectBeginDecimalLongitude))

# add meta data
sg_temps_epi_gz_meta <- left_join(sg_temps_epi_gz, full_meta)
# final step is to add a few calculated seagrass metrics that might be useful for the SEM - Canopy height and Structure
sg_temps_epi_gz_meta$CanopyHeight <- sg_temps_epi_gz_meta$LongestBladeLengthMean + sg_temps_epi_gz_meta$SheathLengthMean
# convert EpiphytePerAreaMean units from g per cm2 to mg per cm2
sg_temps_epi_gz_meta$EpiphytePerAreaMeanMg <- sg_temps_epi_gz_meta$EpiphytePerAreaMean*1000
# calculate log of shoot density
sg_temps_epi_gz_meta$DensityLog <- log10(sg_temps_epi_gz_meta$DensityShootsMean)
# compiled dataset is now ready for modeling etc. 
# write_csv(sg_temps_epi_gz_meta, "data/full_seagrass_epifauna_for_SEM.csv")

tran_vars <- select(sg_temps_epi_gz_meta, "PercentOtherMean", "PercentBareMean","epifauna_log_per_g_transect_all",
                    "epifauna_log_per_g_transect_large", "MonthlyMeanTemp_June" = "MonthlyMeanTemp",
                    "Meadow", "Year", "TidalHeight", "Region", "SiteCode", "Transect", "TransectBeginDecimalLatitude",
                    "TempAnomWarm_June",
                    "PrevalenceMean", "SeverityMean", "LesionAreaMean", "Latitude"="TransectBeginDecimalLatitude",
                    "CanopyHeight", "DensityShootsMean", "SheathLengthMean", "DensityLog", "EpiphytePerAreaMeanMg",
                    "Ampithoid" = "ampithoid_log_per_g_transect_all", "Lacuna" = "lacuna_log_per_g_transect_all",
                    "Idoteid"="idoteid_log_per_g_transect_all", "Richness"="richness_site_all",
                    "Ampithoid_large" = "ampithoid_log_per_g_transect_large", "Lacuna_large" = "lacuna_log_per_g_transect_large",
                    "Idoteid_large"="idoteid_log_per_g_transect_large", "Richness_large"="richness_site_large")

site <- tran_vars %>%
  group_by(Year, Region, SiteCode, Meadow) %>%
  summarise(PercentOther=mean(PercentOtherMean), PercentBare=mean(PercentBareMean), 
            Latitude=mean(Latitude),
            Density=mean(DensityShootsMean), SheathLength=mean(SheathLengthMean), CanopyHeight=mean(CanopyHeight),
            Epifauna_all=mean(epifauna_log_per_g_transect_all), 
            Epifauna_large = mean(epifauna_log_per_g_transect_large),
            DensityLog=log10(Density), Latitude=mean(Latitude),
            EpiphytePerAreaMeadowMg=mean(EpiphytePerAreaMeanMg),
            TempAnomWarm_June=mean(TempAnomWarm_June),
            MonthlyMeanTemp_June=mean(MonthlyMeanTemp_June),
            Prevalence=mean(PrevalenceMean), LesionArea=mean(LesionAreaMean),
            Ampithoid_all=mean(Ampithoid, na.rm=TRUE), Lacuna_all=mean(Lacuna, na.rm=TRUE),
            Idoteid_all=mean(Idoteid, na.rm=TRUE), Richness_all=max(Richness),
            Ampithoid_large=mean(Ampithoid_large, na.rm=TRUE), Lacuna_large=mean(Lacuna_large, na.rm=TRUE),
            Idoteid_large=mean(Idoteid_large, na.rm=TRUE), Richness_large=max(Richness_large))

site_sem <- select(site, c(CanopyHeight, Epifauna_all, Epifauna_large, EpiphytePerAreaMeadowMg,
                           TempAnomWarm_June, MonthlyMeanTemp_June, DensityLog, 
                           Ampithoid_all, Ampithoid_large, Lacuna_all, Lacuna_large,
                           Idoteid_all, Idoteid_large,Latitude,
                           Richness_all, Richness_large))
# divide canopy height by 1000 to convert from mm to m
site_sem$CanopyHeight <- site_sem$CanopyHeight/1000

# add blade-level epiphytes
dis$EpiphytePerAreagcm2 <- dis$EpiphyteDryMass/dis$BladeArea # blade area is already in cm2 for this
dis$EpiphytePerAreamgcm2 <- dis$EpiphyteDryMass/dis$BladeArea*1000 # convert to mg to get better numbers
dis <- dis[-which(dis$EpiphyteDryMass<0),] # removes erroneous mass values
dis_site <- left_join(dis, site_sem)
# add binary variables for year and tidal height
dis_site$TidalHeightBinary <- ifelse(dis_site$TidalHeight=="L", 0, 1)
# restirct to 2019 and 2021
dis_site <- subset(dis_site, Year!=2020)
dis_site$YearBinary <- ifelse(dis_site$Year==2019, 0, 1)
# this is the data set for models that will build SEM
# exclude sites without epifauna (for various reasons, see supplemental table)
dis_site <- dis_site[-which(is.na(dis_site$Epifauna_all)),] # reduces dataset 
dis_site$Meadow <- paste(dis_site$Region, dis_site$SiteCode, sep="_")
write_csv(dis_site, "data/output/epiphyte_SEM_data_all_large.csv")

## alternative data set without temperature ####
## include all sites and years but not temperature for region-specific analyses
## predictors at site level except blade area and epiphyte at blade level
## responses at blade level (disease)

epi_vars <- select(epi, c("year", "region", "site_unique_code", "transect_unique_code", 
                          "Epifauna_large" = "epifauna_log_per_g_transect_large", 
                          "Ampithoid_large" = "ampithoid_log_per_g_transect_large", 
                          "Lacuna_large" = "lacuna_log_per_g_transect_large",
                          "Idoteid_large"="idoteid_log_per_g_transect_large", 
                          "Richness_large"="richness_site_large"))

epi_site <- epi_vars %>%
  group_by(year, region, site_unique_code) %>%
  summarise(Epifauna_large = mean(Epifauna_large),
            Ampithoid_large=mean(Ampithoid_large, na.rm=TRUE), Lacuna_large=mean(Lacuna_large, na.rm=TRUE),
            Idoteid_large=mean(Idoteid_large, na.rm=TRUE), Richness_large=max(Richness_large))
epi_site$Meadow <- str_extract(epi_site$site_unique_code, "[:graph:]{4}")
epi_site$Meadow <- str_replace(epi_site$Meadow, "[.]", "_")

sg$CanopyHeight <- sg$SheathLengthMean + sg$LongestBladeLengthMean

sg_site <- sg %>%
  group_by(Year, Region, SiteCode) %>%
  summarise(Density=mean(DensityShootsMean, na.rm=T), CanopyHeight=mean(CanopyHeight, na.rm=T)/1000,
            DensityLog=log10(Density)) %>%
  mutate(Meadow=paste(Region, SiteCode, sep="_"))

epi_sg_site <- left_join(epi_site, sg_site, by=c("year" = "Year", "region" = "Region", "Meadow"))
epi_sg_site <- rename(epi_sg_site, Region = region, Year = year)

epi_sg_site <- full_join(epi_sg_site, meta_site)

blade_vars <- dis %>%
  select(c("Year", "Region", "SiteCode", "TidalHeight", "Transect", "Blade", "Prevalence", "LesionArea", "BladeArea")) %>%
  mutate(BladeAreaLog=log10(BladeArea), Meadow=paste(Region, SiteCode, sep = "_"))

all_epi <- inner_join(blade_vars, epi_sg_site, by=c("Year", "Region", "Meadow", "SiteCode"))
all_epi <- na.omit(all_epi)
# 2370
write_csv(all_epi, "data/output/epifauna_for_region_specific_models_no_epiphyte.csv")

## meadow scale variables for plots ####

sg_site2 <- sg %>%
  group_by(Year, Region, SiteCode) %>%
  summarise(Density=mean(DensityShootsMean, na.rm=T), CanopyHeight=mean(CanopyHeight, na.rm=T)/1000,
            DensityLog=log10(Density), Prevalence=mean(PrevalenceMean), LesionArea=mean(LesionAreaMean),
            LesionAreaLog=log10(LesionArea)) %>%
  mutate(Meadow=paste(Region, SiteCode, sep="_"), MeadowYear=paste(Meadow, Year, sep="_"))

all_site <- full_join(sg_site2, epi_site, by=c("Year" = "year", "Region" = "region", "Meadow"))
all_site <- full_join(meta_site, all_site)

# monthly$MonthW <- format(monthly$Month, "%B")
june <- subset(monthly, MonthW=="June")
june$Year <- as.integer(format(june$Month, "%Y"))
june <- select(june, c(Year, Region, SiteCode=Site, Meadow, MonthlyMeanTemp_June=MonthlyMeanTemp))

all_site <- full_join(june, all_site)
write_csv(all_site, "data/output/epifauna_site_for_plotting.csv")
