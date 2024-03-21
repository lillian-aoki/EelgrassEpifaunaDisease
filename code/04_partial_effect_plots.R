# 04_partial_effect_plots

# This script plots partial effects for component models of the meadow-scale SEM
# outputs are partial effect plots for the main manuscript

# Partial predictor plots

# librarys ####
library(lme4)
library(dplyr)
library(readr)
library(DHARMa)
library(ggplot2)
library(performance)
library(effects)
library(ggeffects)
library(patchwork)

# data ###
region_order <- c("AK", "BC", "WA", "OR", "BB", "SD")
dis <- read_csv("data/output/epiphyte_SEM_data_all_large.csv")
dis$Region <- ordered(dis$Region, levels=region_order)
# updating the SEM to compare the effects of using large vs all animals
dis$BladeAreaLog <- log10(dis$BladeArea)
dis$EpiphyteLog <- log10(dis$EpiphytePerAreamgcm2+0.01)

dis1_large <- select(dis, c(Epifauna = Epifauna_large, TempAnomWarm_June, MonthlyMeanTemp_June, CanopyHeight, 
                            DensityLog, YearBinary, Year, Meadow, Region, Transect, Blade, BladeAreaLog, TidalHeightBinary, GrazingScars,
                            Prevalence, LesionArea, EpiphyteLog, Lacuna = Lacuna_large, 
                            Ampithoid = Ampithoid_large, Idoteid = Idoteid_large, Richness = Richness_large))
dis1_large$GrazingScars <- as.factor(dis1_large$GrazingScars)
dis_large <- na.omit(dis1_large)
dis_large$Meadow_Year <- paste(dis_large$Meadow, dis_large$Year, sep = "_")
site_large <- distinct(dis_large, Meadow_Year, .keep_all = T)
site_large <- select(site_large, -c(Prevalence, LesionArea, EpiphyteLog, TidalHeightBinary, BladeAreaLog))
les_large <- subset(dis_large, LesionArea>0)
les_large$LesionAreaLog <- log10(les_large$LesionArea)
les_large <- na.omit(les_large)

site_les <- distinct(les_large, Meadow_Year, .keep_all = T)

# prev #####
prev_epi_1 <- glmer(Prevalence ~ BladeAreaLog + EpiphyteLog + GrazingScars + Epifauna + CanopyHeight + DensityLog + 
                      TempAnomWarm_June + MonthlyMeanTemp_June + 
                      TidalHeightBinary + YearBinary +
                      (1|Region) +(1|Meadow), 
                    data=dis_large,
                    family="binomial")
epi1 <- predictorEffect("Epifauna", prev_epi_1, partial.residuals=T)
epi2 <- as.data.frame(epi1)

epi_prev <- ggplot(epi2)+geom_ribbon(aes(x=Epifauna, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Epifauna, y=fit))+
  scale_y_continuous(limits=c(-0, 1), expand=c(0,0))+
  xlab("")+
  ylab("Probability of disease")+
  scale_x_continuous(limits=c(-1,1))+
  labs(tag = "(a)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))
epi_prev

epi_prev_dist <- ggplot()+geom_histogram(data=dis_large, aes(x=Epifauna, fill=Region),bins = 40)+
  xlab("Log epifauna \n(abundance per g macrophytes)")+
  ylab("Count\nleaves")+
  scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  scale_x_continuous(limits=c(-1,1))+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))

a <- epi_prev / epi_prev_dist + plot_layout(heights = c(1,0.3), guides = "collect")
a <- a + theme(legend.position = "")
a
# ggsave(filename = "figures/2a_gz_sem_partial_epifauna_prev.jpg", width=4, height=3)
ggsave(filename = "figures/2a_gz_sem_partial_epifauna_prev.tiff", width=4, height=3)
# okay use ggplot approach - it's a lot nicer than trying to get base r plot to work
prev_lac_1 <- glmer(Prevalence ~ BladeAreaLog + EpiphyteLog + GrazingScars + Lacuna + CanopyHeight + DensityLog + 
                      TempAnomWarm_June + MonthlyMeanTemp_June + 
                      TidalHeightBinary + YearBinary +
                      (1|Region) +(1|Meadow), 
                    data=dis_large,
                    family="binomial")
lac1 <- predictorEffect("Lacuna", prev_lac_1, partial.residuals=T)
lac2 <- as.data.frame(lac1)

lac_prev <- ggplot(lac2)+geom_ribbon(aes(x=Lacuna, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Lacuna, y=fit))+
  scale_y_continuous(limits=c(-0, 1), expand=c(0,0))+
  xlab("")+
  ylab("Probability of disease")+
  labs(tag = "(b)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))
lac_prev

lac_prev_dist <- ggplot()+geom_histogram(data=dis_large, aes(x=Lacuna, fill=Region),bins = 40)+
  xlab(expression(atop(paste("Log ", italic("Lacuna"), " snails"), "(abundance per g macrophytes)")))+
  ylab("Count\nleaves")+
  scale_y_continuous(breaks=c(0,75,150), limits=c(0,150))+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))

b <- lac_prev / lac_prev_dist + plot_layout(heights = c(1,0.3), guides = "collect")
b <- b + theme(legend.position = "")
b
# ggsave("figures/2b_gz_sem_partial_lac.jpg", width=4, height = 3)
ggsave("figures/2b_gz_sem_partial_lac.tiff", width=4, height = 3)

prev_ido_1 <- glmer(Prevalence ~ BladeAreaLog + EpiphyteLog + GrazingScars + Idoteid + CanopyHeight + DensityLog + 
                      TempAnomWarm_June + MonthlyMeanTemp_June + 
                      TidalHeightBinary + YearBinary +
                      (1|Region) +(1|Meadow), 
                    data=dis_large,
                    family="binomial")
ido1 <- predictorEffect("Idoteid", prev_ido_1, partial.residuals=T)
ido2 <- as.data.frame(ido1)

ido_prev <- ggplot(ido2)+geom_ribbon(aes(x=Idoteid, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Idoteid, y=fit))+
  scale_y_continuous(limits=c(-0, 1), expand=c(0,0))+
  xlab("")+
  ylab("Probabilty of disease")+
  scale_x_continuous(breaks=c(-2.0, -1.5, -1, -0.5))+
  labs(tag = "(c)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))
ido_prev

ido_prev_dist <- ggplot()+geom_histogram(data=dis_large, aes(x=Idoteid, fill=Region),bins = 40)+
  xlab("Log idoteid isopods \n(abundance per g macrophytes)")+
  ylab("Count\nleaves")+
  scale_y_continuous(breaks=c(0,100,200))+
  scale_x_continuous(breaks=c(-2.0, -1.5, -1, -0.5))+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(2,2,2,2), "pt"))

c <- ido_prev / ido_prev_dist + plot_layout(heights = c(1,0.3))
c <- c + theme(legend.position = "")
c
# ggsave("figures/2c_gz_sem_partial_ido.jpg", width=4, height=3)
ggsave("figures/2c_gz_sem_partial_ido.tiff", width=4, height=3)

prev_amp_1 <- glmer(Prevalence ~ BladeAreaLog + EpiphyteLog + GrazingScars + Ampithoid + CanopyHeight + DensityLog + 
                      TempAnomWarm_June + MonthlyMeanTemp_June + 
                      TidalHeightBinary + YearBinary +
                      (1|Region) +(1|Meadow), 
                    data=dis_large,
                    family="binomial")
amp1 <- predictorEffect("Ampithoid", prev_amp_1, partial.residuals=T)
amp2 <- as.data.frame(amp1)
amp_prev <- ggplot(amp2)+geom_ribbon(aes(x=Ampithoid, y=fit, ymax=upper, ymin=lower),fill="grey95")+
  geom_line(aes(x=Ampithoid, y=fit), linetype="dashed", color="grey80")+
  # geom_point(data=dis_large, aes(x=Ampithoid, y=0.01, color=Region), shape="|", size=3)+
  geom_text(aes(x=-1, y=0.85), label="N.S.")+
  scale_y_continuous(limits=c(-0, 1), expand=c(0,0))+
  # xlab("Log ampithoid \n(abundance per g macrophytes)")+
  xlab("")+
  ylab("Probability of disease")+
  labs(tag = "(d)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(2,2,2,2), "pt"))
amp_prev

amp_prev_dist <- ggplot()+geom_histogram(data=dis_large, aes(x=Ampithoid, fill=Region),bins = 40)+
  xlab("Log ampithoid amphipods \n(abundance per g macrophytes)")+
  ylab("Count\nleaves")+
  scale_y_continuous(breaks=c(0,100,200))+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        plot.margin = unit(c(2,2,2,2), "pt"))

d <- amp_prev / amp_prev_dist + plot_layout(heights = c(1,0.3))
d <- d + theme(legend.position = "bottom")
d
# ggsave("figures/2d_gz_sem_partial_amp.jpg", width = 4, height=4)
ggsave("figures/2d_gz_sem_partial_amp.tiff", width = 4, height=4)

# lesion area ####
les_epi_1 <- lmer(LesionAreaLog ~ BladeAreaLog + EpiphyteLog + GrazingScars + Epifauna + CanopyHeight + DensityLog + 
                    TempAnomWarm_June + MonthlyMeanTemp_June + 
                    TidalHeightBinary + YearBinary +
                    (1|Region) +(1|Meadow), 
                  data=les_large)

les_lac_1 <- lmer(LesionAreaLog ~ BladeAreaLog + EpiphyteLog + GrazingScars + Lacuna + CanopyHeight + DensityLog + 
                    TempAnomWarm_June + MonthlyMeanTemp_June + 
                    TidalHeightBinary + YearBinary +
                    (1|Region) +(1|Meadow), 
                  data=les_large)

les_amp_1 <- lmer(LesionAreaLog ~ BladeAreaLog + EpiphyteLog + GrazingScars + Ampithoid + CanopyHeight + DensityLog + 
                    TempAnomWarm_June + MonthlyMeanTemp_June + 
                    TidalHeightBinary + YearBinary +
                    (1|Region) +(1|Meadow), 
                  data=les_large)

les_ido_1 <- lmer(LesionAreaLog ~ BladeAreaLog + EpiphyteLog + GrazingScars + Idoteid + CanopyHeight + DensityLog + 
                    TempAnomWarm_June + MonthlyMeanTemp_June + 
                    TidalHeightBinary + YearBinary +
                    (1|Region) +(1|Meadow), 
                  data=les_large)

epi3 <- predictorEffect("Epifauna", les_epi_1, partial.residuals=T)
epi4 <- as.data.frame(epi3)

lac3 <- predictorEffect("Lacuna", les_lac_1, partial.residuals=T)
lac4 <- as.data.frame(lac3)

amp3 <- predictorEffect("Ampithoid", les_amp_1, partial.residuals=T)
amp4 <- as.data.frame(amp3)

ido3 <- predictorEffect("Idoteid", les_ido_1, partial.residuals=T)
ido4 <- as.data.frame(ido3)

epi_les <- ggplot(epi4)+geom_ribbon(aes(x=Epifauna, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Epifauna, y=fit))+
  geom_jitter(data=les_large, aes(x=Epifauna, y=LesionAreaLog, color = Region), width=0.05, height = 0.05, alpha = 0.5)+
  # scale_y_continuous(trans = "logit")+
  xlab("Log epifauna \n(abundance per g macrophytes)")+
  ylab(expression(paste("Log lesion area (",mm^2, ")")))+
  labs(tag = "(e)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank())
epi_les
e <- epi_les + theme(legend.position = "")
e
# ggsave("figures/2e_gz_sem_les_epi.jpg", width = 4, height=3)
ggsave("figures/2e_gz_sem_les_epi.tiff", width = 4, height=3)

ples_lac <- ggplot(lac4)+geom_ribbon(aes(x=Lacuna, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Lacuna, y=fit), color="black")+
  geom_jitter(data=les_large, aes(x=Lacuna, y=LesionAreaLog, color=Region), width=0.05, height = 0.05, alpha = 0.5)+
  xlab(expression(atop(paste("Log ", italic("Lacuna"), " snails"), "(abundance per g macrophytes)")))+
  # xlab("Log lacuna \n(abundance per g macrophytes)")+
  ylab(expression(paste("Log lesion area (",mm^2, ")")))+
  labs(tag = "(f)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        legend.position="bottom")
ples_lac
f <- ples_lac + theme(legend.position = "")
f
# ggsave("figures/2f_gz_sem_les_lac.jpg", width = 4, height=3)
ggsave("figures/2f_gz_sem_les_lac.tiff", width = 4, height=3)

ples_ido <- ggplot(ido4)+geom_ribbon(aes(x=Idoteid, y=fit, ymax=upper, ymin=lower),fill="lightgrey")+
  geom_line(aes(x=Idoteid, y=fit), color="black")+
  geom_jitter(data=les_large, aes(x=Idoteid, y=LesionAreaLog, color=Region), width=0.05, height = 0.05, alpha = 0.5)+
  xlab("Log idoteid isopods \n(abundance per g macrophytes)")+
  ylab(expression(paste("Log lesion area (",mm^2, ")")))+
  # ylab("Log lesion area (mm2)")+
  labs(tag = "(g)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        legend.position="bottom")
ples_ido
g <- ples_ido + theme(legend.position = "blank")
g
# ggsave("figures/2g_gz_sem_les_ido.jpg", width = 4, height=3)
ggsave("figures/2g_gz_sem_les_ido.tiff", width = 4, height=3)

ples_amp <- ggplot(amp4)+geom_ribbon(aes(x=Ampithoid, y=fit, ymax=upper, ymin=lower),fill="grey95")+
  geom_line(aes(x=Ampithoid, y=fit), color="grey80", linetype="dashed")+
  geom_jitter(data=les_large, aes(x=Ampithoid, y=LesionAreaLog, color=Region), width=0.05, height = 0.05, alpha = 0.5)+
  geom_text(aes(x=-1, y=0.85), label="N.S.")+
  xlab("Log ampithoid amphipods \n(abundance per g macrophytes)")+
  ylab(expression(paste("Log lesion area (",mm^2, ")")))+
  labs(tag = "(h)   ")+
  theme_bw(base_size=12)+
  theme(panel.grid = element_blank(),
        legend.position="bottom")
ples_amp
h <- ples_amp # + theme(legend.position = "")
h
# ggsave("figures/2h_gz_sem_les_amp.jpg", width = 4, height=4)
ggsave("figures/2h_gz_sem_les_amp.tiff", width = 4, height=4)
