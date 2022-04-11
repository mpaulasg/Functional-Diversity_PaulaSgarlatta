################################################################################
##
## Script for FD statistical analysis
##
## 
## Code by Paula Sgarlatta
##
################################################################################

### Using data from sites

rm(list=ls()) # cleaning memory

##Load packages

library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(here)


load(here::here("data", "temporal_alpha_kelp_biomass.RData") )
load(here::here("data", "spatial_alpha_biomass.RData") )

## Preparing data for stats

kelp_stats <- temporal_alpha_kelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  dplyr::select(-Site1)

kelp_stats$Year_cont <- as.numeric(kelp_stats$Year)

spatial_stats <- spatial_alpha %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Habitat=sub(".*_", "", Site1), .before="Site1")%>% 
  dplyr::select(-Site1)


####### GLMM

### TEMPORAL ###

## Species richness

kelp_s <- glmmTMB(sp_richn ~ Year + (1|Site), data=kelp_stats, family = poisson())

summary(kelp_s)
Anova(kelp_s)#Significant - p=0.0002

kelp_res_s <- simulateResiduals(kelp_s)
plot(kelp_res_s) #Good

## FRic 

kelp_fric <- glmmTMB(fric ~ Year + (1|Site), data=kelp_stats, family=beta_family())

summary(kelp_fric)
Anova(kelp_fric)#Significant - p=0.004

kelp_res_fric <- simulateResiduals(kelp_fric)
plot(kelp_res_fric) # Good

## FDis 

kelp_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=kelp_stats, family = beta_family())

summary(kelp_fdis)
Anova(kelp_fdis)#Significative - p=0.006

kelp_res_fdis <- simulateResiduals(kelp_fdis)
plot(kelp_res_fdis) # Good

## FIde1

kelp_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide1)
Anova(kelp_fide1) #Significant - p=0.03

kelp_res_fide1 <- simulateResiduals(kelp_fide1)
plot(kelp_res_fide1) # Good

## FIde2

kelp_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide2)
Anova(kelp_fide2) #Significant - p=0.001

kelp_res_fide2 <- simulateResiduals(kelp_fide2)
plot(kelp_res_fide2) # Good

## FIde3

kelp_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide3)
Anova(kelp_fide3) #Significant - p=0.0003

kelp_res_fide3 <- simulateResiduals(kelp_fide3)
plot(kelp_res_fide3) # Good

### SPACE ###

## Species richness

space_s <- glmmTMB (sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = poisson()) 

summary(space_s)
Anova(space_s)#Significant - p=0.025

space_res_s <- simulateResiduals(space_s)
plot(space_res_s) #Good 

## FRic

spatial_fric <- glmmTMB (fric ~ Habitat + (1|Site) , data=spatial_stats, family = beta_family())

summary(spatial_fric)
Anova(spatial_fric)#Significant - p<0.000005

spatial_res_fric <- simulateResiduals(spatial_fric)
plot(spatial_res_fric) # Good

## Fdis

spatial_fdis <- glmmTMB (fdis ~ Habitat, data=spatial_stats, family = beta_family())

summary(spatial_fdis)
Anova(spatial_fdis)# No diferences - p=0.4495

spatial_res_fdis <- simulateResiduals(spatial_fdis)
plot(spatial_res_fdis) # Good

## FIde1

spatial_fide1 <- glmmTMB (fide_PC1 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide1)
Anova(spatial_fide1)# Not significant - p=0.11

spatial_res_fide1 <- simulateResiduals(spatial_fide1)
plot(spatial_res_fide1) # Good

## FIde2 

spatial_fide2 <- glmmTMB (fide_PC2 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide2)
Anova(spatial_fide2) #Not significant - p=0.62

spatial_res_fide2 <- simulateResiduals(spatial_fide2)
plot(spatial_res_fide2) # Good

## FIde3

spatial_fide3 <- glmmTMB (fide_PC3 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide3)
Anova(spatial_fide3) #Significant - p=0.0005

spatial_res_fide3 <- simulateResiduals(spatial_fide3)
plot(spatial_res_fide3) # Good

#################### end of code ##########################################################