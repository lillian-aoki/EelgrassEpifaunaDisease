## 02-Supplemental Figures 

# Last updated: 2024-03-21 by Lillian Aoki

## This script produces supplemental figures for the SEM manuscript: 
## (1) the distribution of epifauna across latitudes,
## (2) comparison of seagrass structure and disease prevalence, 
## (3) a map of study locations with temperature 
## (4) relative abundance of eelgrass grazers within epifauna community
## (5) mosaic plot of grazing scar and disease prevalence
## (6) pathogen loads in lesioned tissue via qPCR analysis

# libraries ####
library(tidyverse)
library(patchwork)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggmosaic)
library(broom)
# data input ####
region_order <- c("AK", "BC", "WA", "OR", "BB", "SD")
site_dat <- read_csv("data/output/epifauna_site_for_plotting.csv")
site_dat$fYear <- ordered(as.factor(site_dat$Year), levels = c("2019", "2020", "2021"))
site_dat$Region <- ordered(site_dat$Region, levels=region_order)

epi <- read_csv("data/input/EGWD_transect_data_v20230307_big_epi_with_BB_transect.csv")
meta <- read_csv("data/input/combined_site_metadata.csv")
sem_dat <- read_csv("data/output/epiphyte_SEM_data_all_large.csv")
conservative <- read_csv("data/input/2021_lesion_qPCR_conservative_results.csv")

# (1) epifauna distributions ####
amp_plot <- ggplot(site_dat, aes(x=Latitude, y=Ampithoid_large, color=fYear))+
  geom_point(size=1.5, alpha=0.75)+
  # xlab("")+
  xlab("Latitude (ºN)")+
  ylab("Log ampithoid amphipods \n(abundance per g macrophytes)")+
  scale_y_continuous(limits=c(-2, 1))+
  scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60))+
  # scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60), trans = "reverse")+
  scale_color_manual(values=c('#CCBB44', '#66CCEE', '#AA3377'))+
  # scale_color_viridis_d()+
  # scale_x_reverse(breaks = c(60, 55, 50, 45, 40, 35, 30))+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", legend.box.spacing = unit(0, "pt"),legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"))
amp_plot

lac_plot <- ggplot(site_dat, aes(x=Latitude, y=Lacuna_large, color=fYear))+geom_point(size=1.5, alpha=0.75)+
  xlab("")+
  ylab(expression(atop(paste("Log ", italic("Lacuna"), " snails"), "(abundance per g macrophytes)")))+
  # scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60), trans = "reverse")+
  scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60))+
  scale_y_continuous(limits=c(-2, 1))+
  scale_color_manual(values=c('#CCBB44', '#66CCEE', '#AA3377'))+
  # scale_color_viridis_d()+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", legend.box.spacing = unit(0, "pt"),legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"))

lac_plot

ido_plot <- ggplot(site_dat, aes(x=Latitude, y=Idoteid_large, color=fYear))+geom_point(size=1.5, alpha=0.75)+
  annotate(geom="text", x = c(33, 38.5, 43.5, 48, 52, 55),
           y=0.4, color = "grey30", size = 2.5,
           label = c("San\nDiego", "Bodega\nBay", "Oregon", "Washington","British\nColumbia", "Alaska"),
           angle = 45)+
  ylab("Log idoteid isopods \n(abundance per g marcophytes)")+
  xlab("")+
  # xlab("Latitude (ºN)")+
  scale_y_continuous(limits=c(-2, 1))+
  # scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60), trans = "reverse")+
  scale_x_continuous(breaks = c(30,35,40, 45, 50, 55, 60))+
  # scale_color_viridis_d()+
  scale_color_manual(values=c('#CCBB44', '#66CCEE', '#AA3377'))+
  # scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377'))+
  # scale_color_manual(values=c('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99'))+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(), 
        legend.position = "bottom", legend.box.spacing = unit(0, "pt"), legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt"))

ido_plot

lac_plot / ido_plot / amp_plot / guide_area() + plot_layout(guides="collect", heights = c(1,1,1,0.1))
# ggsave("figures/taxa_latitude_colors.jpg", width = 4, height =8.25 )
ggsave("figures/taxa_latitude_colors.tiff", width = 4, height =8.25 )

# (2) seagrass structure ####

ggplot(site_dat, aes(x=Density, y=CanopyHeight, color=Region, shape=fYear, size=Prevalence))+geom_point(alpha=0.8)+
  scale_x_continuous(trans = "log10")+
  xlab("Shoot Density (shoots per m2)")+
  ylab("Canopy Height (m)")+
  scale_size(breaks=c(0.25, 0.50, 0.75),labels = c("25%", "50%", "75%"), name="Disease \n prevalence")+
  scale_shape(name="Year")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.direction = "vertical")

# ggsave("figures/seagrass_structure_prevalence.jpg", width = 4, height = 6)
ggsave("figures/seagrass_structure_prevalence.tiff", width = 4, height = 6)

# (3) site map ####
sites <- na.omit(site_dat)
sites <- subset(sites, Year=="2019")
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant")

world <- ne_countries(scale = "medium", returnclass = "sf")

mapFull <- ggplot(data = world) + 
  geom_sf(fill="grey90")+
  # geom_sf(data=sites_sf, size=2, aes(color=MonthlyMeanTemp_June), alpha=0.5, nudge_y = 0.1) + 
  geom_sf(data=sites_sf, size=2, aes(color=MonthlyMeanTemp_June), alpha=0.7) + 
  # geom
  coord_sf(xlim = c(max(sites$Longitude), min(sites$Longitude-4)), 
           ylim = c(max(sites$Latitude), min(sites$Latitude)),
           default_crs = sf::st_crs(4326))+
  scale_color_continuous(type="viridis", name="Mean Daily\nSea Surface Temperature\n(June 2019)",
                         breaks=c(10, 12, 14, 16, 18, 20))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.position = "bottom", plot.margin = unit(c(1,2,1,2), "pt"), legend.title.align = 0.5,
        legend.title = element_text(size=11))

mapFull + 
  # geom_sf_label(data=regionsf, aes(label=label))
  annotate(geom = "label", color = "grey30", size = 3.5,
           x = c(-136, -133, -128, -127, -127, -120.5), 
           y = c(55, 51, 48, 44, 38, 32.5),
           fontface = "italic",
           label = c("Alaska", "British Columbia", "Washington", "Oregon", "Bodega Bay", "San Diego"))

ggsave("figures/map_study_sites.jpg", width = 4, height = 6)
ggsave("figures/map_study_sites.tiff", width = 4, height = 6)

# (4) relative abundance of eelgrass grazers ####
#transect level epifauna - need to collapse to site
meta19 <- subset(meta, Year=="2019")
loc19 <- select(meta19, c("Region","SiteCode", "TidalHeight", "Transect", "TransectBeginDecimalLatitude",
                          "TransectBeginDecimalLongitude", "TransectEndDecimalLatitude", "TransectEndDecimalLongitude"))
full_meta <- meta %>%
  rows_update(loc19, by = c("Region", "SiteCode", "TidalHeight", "Transect"))
meta_site <- full_meta %>%
  group_by(Year, Region, SiteCode) %>%
  summarise(Latitude=mean(as.numeric(TransectBeginDecimalLatitude)), Longitude=mean(TransectBeginDecimalLongitude))

epi$site <- str_extract(epi$site_unique_code, "\\.[A-F]\\.")
epi$site <- gsub("\\.", "", epi$site)
epi <- left_join(epi, meta_site, by=c("year"="Year", "region" = "Region", "site" = "SiteCode"))
epi$region <- ordered(epi$region, levels=region_order)
epi_summ <- epi %>%
  group_by(year, region, site, site_unique_code, Latitude, meadow=paste(region, site, sep="_")) %>%
  mutate(meadow=paste(region, site, sep="_"))%>%
  summarise(total_abun_large=mean(epifauna_per_g_transect_large),
            lac_abun_large=mean(lacuna_per_g_transect_large),
            amp_abun_large=mean(ampithoid_per_g_transect_large),
            ido_abun_large=mean(idoteid_per_g_transect_large),
            other_abun_large=total_abun_large-sum(lac_abun_large, amp_abun_large, ido_abun_large))
epi_summ_long <- pivot_longer(epi_summ, cols = c("other_abun_large", "lac_abun_large", "amp_abun_large", "ido_abun_large"),
                              names_to = "taxon") %>%
  arrange(Latitude)
meadow_order <- c("AK_A", "AK_B", "AK_C", "AK_D", "AK_E", "AK_F",
                  "BC_A", "BC_B", "BC_C", "BC_D", "BC_E",
                  "WA_A", "WA_B", "WA_C", "WA_D", "WA_E", 
                  "OR_A", "OR_B", "OR_C", "OR_D", "OR_E", 
                  "BB_A", "BB_B", "BB_C", "BB_D", "BB_E", "BB_F",
                  "SD_A", "SD_B", "SD_C", "SD_D", "SD_E")
epi_summ_long$meadow <- ordered(epi_summ_long$meadow, levels=meadow_order)
epi_summ$other_per <- epi_summ$other_abun_large/epi_summ$total_abun_large
epi_summ$gz_per <- 1-epi_summ$other_per
epi_summ$lac_per <- epi_summ$lac_abun_large/epi_summ$total_abun_large
epi_summ$ido_per <- epi_summ$ido_abun_large/epi_summ$total_abun_large
epi_summ$amp_per <- epi_summ$amp_abun_large/epi_summ$total_abun_large
summary(epi_summ$gz_per)
epi_summ %>% filter(gz_per>0.5) %>% filter(total_abun_large>1) %>% select(c(site_unique_code, lac_per, ido_per, amp_per))
epi_summ %>% filter(lac_per>0.15) %>% filter(total_abun_large>1) %>% select(c(site_unique_code, lac_per, ido_per, amp_per))

ggplot(epi_summ_long, aes(x=meadow, y=value, fill=taxon))+
  geom_rect(aes(xmin="BB_A", xmax="BB_F", ymin=0, ymax=20), linetype="dashed", fill="grey90")+
  geom_col()+
  facet_grid(rows = "year")+
  ylab("Count of animals per g macrophytes")+
  scale_fill_viridis_d(labels=c("ampithiod \namphipods", "idoteid \nisopods", expression(atop(italic("Lacuna"), "snails")), "other \ntaxa"))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw(base_size=14)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, vjust = 1, size=9, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(hjust = 0))

ggsave("figures/rel_abun.tiff", width=8, height = 6)

# (5) grazing-prevalence mosaic ####
sem_dat$fGrazing <- as.factor(sem_dat$GrazingScars)
sem_dat$fGrazing <- recode_factor(sem_dat$fGrazing, "0" = "Absent", "1" = "Present")
sem_dat$fYear <- as.factor(sem_dat$Year)
sem_dat$fPrevalence <- as.factor(sem_dat$Prevalence)
sem_dat$fPrevalence <- recode_factor(sem_dat$fPrevalence, "0" = "Healthy", "1" = "Diseased")
sem_dat <- drop_na(sem_dat, c("Prevalence", "GrazingScars"))

mosaic <- ggplot(sem_dat)+
  geom_mosaic(aes(x=product(fPrevalence, fGrazing), fill= fPrevalence))+
  geom_mosaic_text(aes(x = product(fPrevalence, fGrazing), label = after_stat(.wt)), as.label=TRUE, size = 4)+
  theme_bw(base_size=14)+
  scale_fill_manual(labels=c("Healthy", "Diseased"), values = c('#228833', '#CCBB44'))+
  xlab("Grazing scars")+
  ylab("Wasting disease infection")+
  labs(tag = "(b)")+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        # plot.margin = unit(0, "pt"),
        legend.position = "bottom",
        legend.text = element_text(size=10))
mosaic / guide_area() + plot_layout(guides= "collect", heights = c(1, 0.1))

# ggsave("figures/4b_grazing_scar_mosaic_colors.jpg", width = 4, height = 3.5)
ggsave("figures/4b_grazing_scar_mosaic_colors.tiff", width = 4, height = 3.5)

# proportion test for grazing and disease over full dataet
length(which(sem_dat$GrazingScars==0)) # n of 'trials' without grazing
length(which(sem_dat$GrazingScars==1)) # n of 'trials' with grazing
length(which(sem_dat$GrazingScars==0 & sem_dat$Prevalence==1)) # n of 'successes' i.e. disease present without grazing
length(which(sem_dat$GrazingScars==1 & sem_dat$Prevalence==1)) # n of 'successes' i.e. disease present with grazing
prop.test(x=c(327, 247), n=c(931, 384))

# (6) pathogen loads ####
con_summ <- conservative %>%
  group_by(Region, Site, Tissue) %>%
  summarise(mean_cells_per_mg=mean(cells_per_mg), count_positive=length(which(cells_per_mg>0)), count_samples=length(cells_per_mg),
            sd_cells_per_mg=sd(cells_per_mg), se_cells_per_mg=sd_cells_per_mg/sqrt(count_samples)) %>%
  mutate(prop_positive=count_positive/count_samples, se_prop_positive = sqrt(prop_positive*(1-prop_positive)/count_samples))
con_summ$Region <- ordered(con_summ$Region, levels=region_order)

lz_prev <- ggplot(con_summ[con_summ$Tissue=="L",], aes(x=Region, color=Site, y=prop_positive, ymin=prop_positive-se_prop_positive, ymax=prop_positive+se_prop_positive))+
  geom_errorbar(position=position_dodge(width = 0.5))+
  geom_point(position=position_dodge(width = 0.5), size=4)+
  scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377'))+
  scale_x_discrete(labels=c("San\nDiego","Bodega\nBay","Oregon","Washington","British\nColumbia","Alaska"), drop=FALSE, limits=rev)+
  ylab(expression(paste(italic("L. zosterae")," prevalence")))+
  labs(tag = "(a)   ")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text.x = element_text(size=10, angle = 45, vjust=0.75))
lz_prev

lz_int <- ggplot(con_summ[con_summ$Tissue=="L",], aes(x=Region, color= Site, y=mean_cells_per_mg))+
  geom_hline(yintercept = mean(conservative$cells_per_mg[conservative$Tissue=="L"]), linetype="dashed", color="darkgrey")+ # global mean of dataset
  geom_errorbar(data= con_summ[con_summ$Tissue=="L",], aes(ymax=mean_cells_per_mg+se_cells_per_mg, ymin=mean_cells_per_mg-se_cells_per_mg), position=position_dodge(width = 0.5))+
  geom_point(data= con_summ[con_summ$Tissue=="L",], position=position_dodge(width = 0.5), size=4)+
  scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377'))+
  scale_x_discrete(labels=c("San\nDiego","Bodega\nBay","Oregon","Washington","British\nColumbia","Alaska"), drop=FALSE, limits=rev)+
  ylab(expression(paste(italic("L. zosterae "), "cells per mg leaf tissue")))+
  labs(tag = "(b)   ")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text.x = element_text(size=10, angle = 45, vjust=0.75))
lz_int

lz_prev / lz_int + plot_layout(guides="collect")
ggsave(filename = "figures/lz_qpcr_results_site.tiff", width = 6, height = 6)
