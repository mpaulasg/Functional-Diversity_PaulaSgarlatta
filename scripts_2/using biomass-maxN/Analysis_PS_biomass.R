################################################################################
##
## Script for FD statistical analysis
##
## 
## Paula Sgarlatta
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


load(here::here("outputs", "using biomass-maxN",  "temporal_alpha_kelp_biomass.RData") )
load(here::here("outputs", "using biomass-maxN", "spatial_alpha_biomass.RData") )
load(here::here("outputs", "using biomass-maxN",  "shift3D_kelp_stats.RData") )
load(here::here("outputs", "using biomass-maxN", "shift3D_space_stats.RData") )

## Preparing data for stats

kelp_stats <- temporal_alpha_kelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)

kelp_stats$Year_cont <- as.numeric(kelp_stats$Year)

spatial_stats <- spatial_alpha %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Habitat=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)


####### GLMM

### KELP ###

kelp_s <- glmmTMB(sp_richn ~ Year + (1|Site), data=kelp_stats, family = poisson())

summary(kelp_s)
Anova(kelp_s)#Significant - p=0.0002

system.time(mod_1 <- drop1(kelp_s, test = "Chisq", all.cols=TRUE))
print(mod_1)

kelp_res_s <- simulateResiduals(kelp_s)
plot(kelp_res_s) #Good


### SPACE ###

space_s <- glmmTMB (sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = poisson()) 

summary(space_s)
Anova(space_s)#Significant - p=0.025

space_res_s <- simulateResiduals(space_s)
plot(space_res_s) #Good 

system.time(space_mod_s <- drop1(space_s, test = "Chisq", all.cols=TRUE))
print(space_mod_s)

### KELP - FRic ###

kelp_fric <- glmmTMB(fric ~ Year + (1|Site), data=kelp_stats, family=beta_family())

summary(kelp_fric)
Anova(kelp_fric)#Significant - p=0.004

kelp_res_fric <- simulateResiduals(kelp_fric)
plot(kelp_res_fric) # Good

system.time(kelp_mod_fric <- drop1(kelp_fric, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fric) #p=0.01

### SPATIAL - FRic ###

spatial_fric <- glmmTMB (fric ~ Habitat + (1|Site) , data=spatial_stats, family = beta_family())

summary(spatial_fric)
Anova(spatial_fric)#Significant - p<0.000005

spatial_res_fric <- simulateResiduals(spatial_fric)
plot(spatial_res_fric) # Good

system.time(spatial_mod_fric <- drop1(spatial_fric, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fric)# p<0.00005

### KELP - FDis ###

kelp_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=kelp_stats, family = beta_family())

summary(kelp_fdis)
Anova(kelp_fdis)#Significative - p=0.006

kelp_res_fdis <- simulateResiduals(kelp_fdis)
plot(kelp_res_fdis) # Good

system.time(kelp_mod_fdis <- drop1(kelp_fdis, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fdis) # p=0.02

### SPATIAL - Fdis ###

spatial_fdis <- glmmTMB (fdis ~ Habitat, data=spatial_stats, family = beta_family())

summary(spatial_fdis)
Anova(spatial_fdis)# No diferences - p=0.4495

spatial_res_fdis <- simulateResiduals(spatial_fdis)
plot(spatial_res_fdis) # Good

system.time(spatial_mod_fdis <- drop1(spatial_fdis, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fdis) # p=0.47


### KELP - Fide1 ###

kelp_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide1)
Anova(kelp_fide1) #Significant - p=0.03

kelp_res_fide1 <- simulateResiduals(kelp_fide1)
plot(kelp_res_fide1) # Good

system.time(kelp_mod_fide1 <- drop1(kelp_fide1, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide1) #p=0.08 - not significant

### SPATIAL - FIde1 ###

spatial_fide1 <- glmmTMB (fide_PC1 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide1)
Anova(spatial_fide1)# Not significant - p=0.11

spatial_res_fide1 <- simulateResiduals(spatial_fide1)
plot(spatial_res_fide1) # Good


system.time(spatial_mod_fide1 <- drop1(spatial_fide1, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide1)

 
### KELP - Fide2 ###

kelp_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide2)
Anova(kelp_fide2) #Significant - p=0.001

kelp_res_fide2 <- simulateResiduals(kelp_fide2)
plot(kelp_res_fide2) # Good

system.time(kelp_mod_fide2 <- drop1(kelp_fide2, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide2) #p=0.012

### SPATIAL - FIde2 ###

spatial_fide2 <- glmmTMB (fide_PC2 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide2)
Anova(spatial_fide2) #Not significant - p=0.62

spatial_res_fide2 <- simulateResiduals(spatial_fide2)
plot(spatial_res_fide2) # Good

system.time(spatial_mod_fide2 <- drop1(spatial_fide2, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide2) # p=0.64

### KELP - Fide3 ###

kelp_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide3)
Anova(kelp_fide3) #Significant - p=0.0003

kelp_res_fide3 <- simulateResiduals(kelp_fide3)
plot(kelp_res_fide3) # Good

system.time(kelp_mod_fide3 <- drop1(kelp_fide3, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide3) #p=0.004

### SPATIAL - FIde3 ###

spatial_fide3 <- glmmTMB (fide_PC3 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide3)
Anova(spatial_fide3) #Significant - p=0.0005

spatial_res_fide3 <- simulateResiduals(spatial_fide3)
plot(spatial_res_fide3) # Good

system.time(spatial_mod_fide3 <- drop1(spatial_fide3, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide3)#p=0.01

## Shift FIde 


#Kelp

shift3D_kelp_stats$Year_cont <- as.numeric(shift3D_kelp_stats$Year)

kelp_shift <- lmer (shift3D ~ Year_cont + (1|Site1), data=shift3D_kelp_stats)
summary(kelp_shift)
Anova(kelp_shift)

system.time(shift_mod <- drop1(kelp_shift, test = "Chisq", all.cols=TRUE))
print(shift_mod)

#Space

#space_shift <- lmer (shift3D ~ Habitat + (1|Site)  , data=shift3D_space_stats)

summary(space_shift)

Anova(space_shift)

system.time(shift_mod_space <- drop1(space_shift, test = "Chisq", all.cols=TRUE))
print(shift_mod_space)



#################### end of code ##########################################################




### Using transects

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


load(here::here("outputs", "using biomass-maxN", "temporal_alpha_kelp_transect.RData") )
load(here::here("outputs", "using biomass-maxN", "spatial_alpha_transect.RData") )
# load(here::here("outputs", "shift3D_kelp_stats.RData") )
# load(here::here("outputs", "shift3D_space_stats.RData") )

kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv"))
spatial_metadata <- read.csv(here::here("from_paula", "SpatialUVC_metadata_transect.csv"))

## Preparing data for stats

kelp_stats <- temporal_alpha_kelp_transect %>%
  rownames_to_column(var = "Code") %>% 
  inner_join(kelp_metadata[ , c("Code", "Site", "Year")], by = c("Code"), all.x=TRUE) %>% 
  select(-Code)
kelp_stats$Year_cont <- as.numeric(kelp_stats$Year)
spatial_stats <- spatial_alpha_transect %>%
  rownames_to_column(var = "Code") %>% 
  inner_join(spatial_metadata[ , c("Code", "Site", "Habitat")], by = c("Code"), all.x=TRUE) %>% 
  select(-Code)


####### GLMM

### KELP ###

kelp_s <- glmmTMB(sp_richn ~ Year + (1|Site), data=kelp_stats, family = poisson())

summary(kelp_s)
Anova(kelp_s)#Significant - p=<0.0000005

kelp_res_s <- simulateResiduals(kelp_s)
plot(kelp_res_s) #Good


### SPACE ###

space_s <- glmmTMB (sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = nbinom2()) 

summary(space_s)
Anova(space_s)#Not significant - p=0.12

largeModel <- glmer.nb(sp_richn ~Habitat + (1|Site), data=spatial_stats)

smallModel <- glmer.nb(sp_richn ~1+(1|Site), data=spatial_stats)



library(pbkrtest)
PBmodcomp(largeModel, smallModel, nsim = 50)

space_res_s <- simulateResiduals(space_s)
plot(space_res_s) # Good 

### KELP - FRic ###

kelp_fric <- glmmTMB(fric ~ Year + (1|Site), data=kelp_stats, family=beta_family())

summary(kelp_fric)
Anova(kelp_fric)#Significant - p<0.0000005

kelp_res_fric <- simulateResiduals(kelp_fric)
plot(kelp_res_fric) # Good

### SPATIAL - FRic ###

spatial_fric <- glmmTMB (fric ~ Habitat + (1|Site) , data=spatial_stats, family = beta_family())

summary(spatial_fric)
Anova(spatial_fric)#Not significant - p=0.16

spatial_res_fric <- simulateResiduals(spatial_fric)
plot(spatial_res_fric) # Good

### KELP - FDis ###

kelp_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=kelp_stats, family = beta_family())

summary(kelp_fdis)
Anova(kelp_fdis)#No Significative - p=0.5

kelp_res_fdis <- simulateResiduals(kelp_fdis)
plot(kelp_res_fdis) # Not good

resid <- residuals(simulateResiduals(kelp_fdis), quantileFunction = qnorm, outliers = c(-7,7))
plot(resid~predict(kelp_fdis))

#Post-hoc tests

pairs(emmeans(kelp_fdis, spec=~Year, type="response"))

#No differences, strange


### SPATIAL - Fdis ###

spatial_fdis <- glmmTMB (fdis ~ Habitat, data=spatial_stats, family = beta_family())

summary(spatial_fdis)
Anova(spatial_fdis)# No diferences - p=0.88

spatial_res_fdis <- simulateResiduals(spatial_fdis)
plot(spatial_res_fdis) # Good

### KELP - Fide1 ###

kelp_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide1)
Anova(kelp_fide1) #No significant - p=0.39

kelp_res_fide1 <- simulateResiduals(kelp_fide1)
plot(kelp_res_fide1) # Good

### SPATIAL - FIde1 ###

spatial_fide1 <- glmmTMB (fide_PC1 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide1)
Anova(spatial_fide1)# Not significant - p=0.25

spatial_res_fide1 <- simulateResiduals(spatial_fide1)
plot(spatial_res_fide1) # Good


### KELP - Fide2 ###

kelp_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide2)
Anova(kelp_fide2) #Not significant - p=0.5

kelp_res_fide2 <- simulateResiduals(kelp_fide2)
plot(kelp_res_fide2) # Good

### SPATIAL - FIde2 ###

spatial_fide2 <- glmmTMB (fide_PC2 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide2)
Anova(spatial_fide2) #Not significant - p=0.65

spatial_res_fide2 <- simulateResiduals(spatial_fide2)
plot(spatial_res_fide2) # Good

### KELP - Fide3 ###

kelp_fide3 <- glmmTMB(fide_PC3 ~ Year_cont + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide3)
Anova(kelp_fide3) #Significant - p=0.02

kelp_res_fide3 <- simulateResiduals(kelp_fide3)
plot(kelp_res_fide3) # Good

### SPATIAL - FIde3 ###

spatial_fide3 <- glmmTMB (fide_PC3 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide3)
Anova(spatial_fide3) #Significant - p=0.01

spatial_res_fide3 <- simulateResiduals(spatial_fide3)
plot(spatial_res_fide3) # Good

## Shift FIde 


#Kelp

shift3D_kelp_stats$Year_cont <- as.numeric(shift3D_kelp_stats$Year)

kelp_shift <- lmer (shift3D ~ Year_cont + (1|Site1), data=shift3D_kelp_stats)
summary(kelp_shift)
Anova(kelp_shift)



#Space

space_shift <- glmmTMB (shift3D ~ Habitat , data=shift3D_space_stats)

summary(space_shift)

Anova(space_shift)

#################### end of code ##########################################################
