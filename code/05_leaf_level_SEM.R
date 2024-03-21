# 05_leaf_level_SEM

# This script analyses the direction of the relationship between grazing scar presence
# and disease metrics in a simplified leaf-level SEM

# libraries ####
library(here)
# note, the two files helpers.R and coefs.R are from the piecewiseSEM package (available on github: http://jslefche.github.io/piecewiseSEM/)
# this analysis was conducted with piecewiseSEM version 2.1.0; however the developer versions of these files
# were required to calculate standardized coefficients for GLMMs. 
# this issue is likely superseded in new versions of the package but for continuity the scripts are included here. 
source(here("code/helpers.R"))
source(here("code/coefs.R"))
library(lme4)
library(dplyr)
library(readr)
library(piecewiseSEM)
# data input ####
dis <- read_csv("data/output/epiphyte_SEM_data_all_large.csv")
# same data as for complex SEM because no grazing scars in 2020 except WA
dis$BladeAreaLog <- log10(dis$BladeArea)
dis_gz <- select(dis, c(Year, Meadow, Region, Blade, BladeAreaLog, TidalHeightBinary, YearBinary, GrazingScars, Prevalence, LesionArea))
dis_gz1 <- na.omit(dis_gz)
dis_gz1$Meadow_Year <- paste(dis_gz1$Meadow, dis_gz1$Year, sep="_")
print(paste("meadow-year combos in prevalence SEM = ", length(unique(dis_gz1$Meadow_Year))))
unique(dis_gz1$Meadow_Year)
print(paste("nrows of disease data (sample size) = ", nrow(dis_gz1)))
les_gz <- subset(dis_gz1, LesionArea>0)
les_gz$LesionAreaLog <- log10(les_gz$LesionArea)
les_gz <- na.omit(les_gz)
print(paste("sample size for lesion models = ", nrow(les_gz)))

# prevalence: grazing to disease ####
sem_prev_gz <- psem(
  lmer(BladeAreaLog ~ TidalHeightBinary + YearBinary + 
         (1|Region) + (1|Meadow),
       data=dis_gz1),
  glmer(GrazingScars ~ BladeAreaLog + 
          YearBinary + 
          (1|Region) + (1|Meadow),
        data=dis_gz1,
        family = "binomial"),
  glmer(Prevalence ~ BladeAreaLog + GrazingScars +
          TidalHeightBinary + YearBinary + 
          (1|Region) + (1|Meadow),
        data=dis_gz1,
        family = "binomial")
)
summary(sem_prev_gz, conserve = T)
coefs(sem_prev_gz)

# prevalence: disease to grazing ####
sem_prev_dis <- psem(
  lmer(BladeAreaLog ~ TidalHeightBinary + YearBinary + 
         (1|Region) + (1|Meadow),
       data=dis_gz1),
  glmer(GrazingScars ~ BladeAreaLog + Prevalence + 
          YearBinary + 
          (1|Region) + (1|Meadow),
        data=dis_gz1,
        family = "binomial"),
  glmer(Prevalence ~ BladeAreaLog +
          TidalHeightBinary + YearBinary + 
          (1|Region) + (1|Meadow),
        data=dis_gz1,
        family = "binomial")
)
summary(sem_prev_dis)
coefs(sem_prev_dis)

# AIC is greater for dis --> gz (40.78) vs gz --> dis (37.06) - 
# suggests that the gz --> disease model better fits the data. 
# Standardized coefficients are similar (0.2227 vs 0.2672); 
# p-values are equivalent (0.0000). Positive relationship - 
# either disease infection increases likelihood of grazing, 
# or presence of grazing increases likelihood of disease. 
# Based on the AIC improvement, there is some evidence for 
# more support for the pathway from grazing to disease. 

# lesion area: grazing to disease ####
sem_les_dis <- psem(
  lmer(BladeAreaLog ~ TidalHeightBinary + YearBinary +
         (1|Region) + (1|Meadow),
       data=les_gz),
  glmer(GrazingScars ~ BladeAreaLog + 
          YearBinary + 
          (1|Region) + (1|Meadow),
        data=les_gz,
        family = "binomial"),
  lmer(LesionAreaLog ~ BladeAreaLog + GrazingScars + 
         TidalHeightBinary + YearBinary + 
         (1|Region) + (1|Meadow),
       data=les_gz)
)
summary(sem_les_dis)
coefs(sem_les_dis)
# lesion area: disease to grazing ####
sem_les_gz <- psem(
  lmer(BladeAreaLog ~ TidalHeightBinary + YearBinary +
         (1|Region) + (1|Meadow),
       data=les_gz),
  glmer(GrazingScars ~ BladeAreaLog + LesionAreaLog + 
          YearBinary + 
          (1|Region) + (1|Meadow),
        data=les_gz,
        family = "binomial"),
  lmer(LesionAreaLog ~ BladeAreaLog + 
         TidalHeightBinary + YearBinary + 
         (1|Region) + (1|Meadow),
       data=les_gz)
)
summary(sem_les_gz)
coefs(sem_les_gz)