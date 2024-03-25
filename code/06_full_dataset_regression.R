# 06_full_dataset_regression

# This script runs a regression analysis using all data collected but without temperature predictors
# Regression analysis confirms the relationships observed in the SEM

library(lme4)
library(dplyr)
library(readr)
library(performance)
library(effects)
library(ggplot2)
library(sjPlot)
library(ggeffects)
library(knitr)
library(DHARMa)
# library(ggplot2)
# data input ####
all_epi <- read_csv("data/output/epifauna_for_region_specific_models_no_epiphyte.csv")
all_epi$fYear <- as.factor(all_epi$Year)
all_epi$Meadow_year <- paste(all_epi$Meadow, all_epi$Year, sep= "_")
length(unique(all_epi$Meadow_year))
unique(all_epi$Meadow_year)
les <- subset(all_epi, LesionArea>0)
les$LesionAreaLog <- log10(les$LesionArea)

# prevalence model comparison 


prev1 <- glmer(Prevalence ~ BladeAreaLog + 
                 Ampithoid_large + Lacuna_large + Idoteid_large +
                 DensityLog + CanopyHeight +
                 fYear +(1|Region) +(1|Meadow), 
               family = "binomial", 
               data=all_epi)

prev2 <- glmer(Prevalence ~ BladeAreaLog + 
                 Ampithoid_large + Lacuna_large + Idoteid_large +
                 
                 fYear +(1|Region) +(1|Meadow), 
               family = "binomial", 
               data=all_epi)

prev3 <- glmer(Prevalence ~ 
                 Ampithoid_large + Lacuna_large + Idoteid_large +
                 DensityLog + CanopyHeight +
                 
                 fYear +(1|Region) +(1|Meadow), 
               family = "binomial", 
               data=all_epi)

prev4 <- glmer(Prevalence ~ BladeAreaLog + 
                 DensityLog + CanopyHeight +
                 fYear +(1|Region) +(1|Meadow), 
               family = "binomial", 
               data=all_epi)

AIC(prev1, prev2, prev3, prev4)


df.AIC <- AIC(prev1, prev2, prev3, prev4)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- format(df.AIC$likelihood/sum(df.AIC$likelihood), scientific = FALSE)
df.AIC$predictors <- c("Full model (individual leaf area, all taxa abundances, seagrass structure)", 
                       "Drop seagrass structure from full model",
                       "Drop individual leaf area from full model", 
                       "Drop taxa abundances from full model") 
df.AIC
kable(df.AIC, digits=4, caption = "Model comparison for prevalence")

prev1_E <- simulateResiduals(prev1)
plot(prev1_E)

plot_model(prev1, type="std", show.p=T, show.values=T, title = "",value.offset = 0.5)
plot(predictorEffects(prev1))
plot(predictorEffects(prev1, partial.residuals=T))
performance(prev1)

prev_names <- c("Year\n 2020", "Lacuna\n snails","Leaf\n area","Year\n 2021","Shoot\n density", "Canopy\n height", 
                "ampithoid\n amphipods", "idoteid\n isopods")
prev_names <- rev(prev_names)
prev_plot <- plot_model(prev1, 
                        type="std", sort.est = TRUE,
                        p.threshold = c(0.05,0,0),
                        show.p=T, 
                        show.values=T, 
                        title = "", 
                        value.offset = 0.3,
                        value.size = 3,
                        axis.labels = prev_names,
                        axis.lim = c(0.5,5))
prev_fig <- prev_plot + theme_bw()+
  geom_hline(yintercept=1, linetype="dashed", color="darkgrey")+
  scale_y_log10(limits=c(0.4,2.5), breaks=c(0.5,1,2))+
  scale_color_viridis_d()+
  # scale_color_manual(values=c('#AA3377', '#66CCEE'))+
  # scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377'))+
  # ylab("Something")+
  # labs(ylab="Scaled estimates of \ndisease prevalence \nodds ratio")+
  labs(tag = "(a)"  )+
  ylab("Scaled estimates of \ndisease prevalence \nodds ratio")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=10))
prev_fig

ggsave("figures/3a_3y_prev_effect_size.tiff", width = 3.35, height = 3.5)

# lesion area ####
les1 <- lmer(LesionAreaLog ~ BladeAreaLog +
               Ampithoid_large + Lacuna_large + Idoteid_large +
               DensityLog + CanopyHeight + 
               fYear + (1|Meadow) + (1|Region), 
             data=les)

les2 <- lmer(LesionAreaLog ~ BladeAreaLog +
               Ampithoid_large + Lacuna_large + Idoteid_large +
               fYear + (1|Meadow) + (1|Region), 
             data=les)

les3 <- lmer(LesionAreaLog ~ 
               Ampithoid_large + Lacuna_large + Idoteid_large +
               DensityLog + CanopyHeight + 
               fYear + (1|Meadow) + (1|Region), 
             data=les)
les4 <- lmer(LesionAreaLog ~ BladeAreaLog +
               DensityLog + CanopyHeight + 
               fYear + (1|Meadow) + (1|Region), 
             data=les)
AIC(les1, les2, les3, les4)


df.AIC <- AIC(les1, les2, les3, les4)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC$predictors <- c("Full model (individual leaf area, all taxa abundances, seagrass structure)", 
                       "Drop seagrass structure from full model",
                       "Drop individual leaf area from full model", 
                       "Drop taxa abundances from full model")
kable(df.AIC, digits = 4, caption = "Model comparison for lesion area")

les2_E <- simulateResiduals(les2)
plot(les2_E)

plot_model(les2, type="std", show.p=T, show.values=T)
plot(predictorEffects(les2, partial.residuals=T))
plot(predictorEffects(les2))
performance(les2)

les2 <- lmer(LesionAreaLog ~ BladeAreaLog +
               Ampithoid_large + Lacuna_large + Idoteid_large +
               fYear + (1|Meadow) + (1|Region), 
             data=les)
les_names <- c("ampithoid\n amphipods", "Year\n 2021", "idoteid\n isopods", "Leaf\n area", "Lacuna\n snails", "Year\n 2020")
les_plot <- plot_model(les2, 
                       type="std",
                       sort.est = T,
                       p.threshold = c(0.05,0,0),
                       show.p=T, 
                       show.values=T,
                       axis.labels = les_names,
                       title = "", 
                       value.offset = 0.3,
                       value.size = 3)
les_fig <- les_plot + theme_bw()+
  geom_hline(yintercept=0, linetype="dashed", color="darkgrey")+
  scale_y_continuous(limits=c(-0.5,0.5))+
  scale_color_viridis_d(begin = 0.6, end = 0)+
  ylab("Scaled estimates of \nlog lesion area\n")+
  labs(tag = "(b)  ")+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=10))
les_fig
ggsave("figures/3b_3y_les_effect_size.tiff", width = 3.35, height = 3.5)

