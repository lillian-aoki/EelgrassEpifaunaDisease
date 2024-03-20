## Figures for epifauna paper

# Last updated: 2024-03-19 by Lillian Aoki

## This script produces figures for (1) the distribution of epifauna across latitudes,
## (2) comparison of seagrass structure and disease prevalence, 
## and (3) a map of study locations with temperature. 

library(tidyverse)
library(patchwork)
library(tidyverse)
library(sf)
library(rnaturalearth)

# data input ####
region_order <- c("AK", "BC", "WA", "OR", "BB", "SD")
site_dat <- read_csv("data/output/epifauna_site_for_plotting.csv")
site_dat$fYear <- ordered(as.factor(site_dat$Year), levels = c("2019", "2020", "2021"))
site_dat$Region <- ordered(site_dat$Region, levels=region_order)

# epifauna distributions ####
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

# seagrass structure ####

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

## site map ####
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
  scale_color_continuous(type="viridis", name="Mean Daily Temperature\n(June 2019)",
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