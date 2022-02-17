################################################################################
##
## Script for FD statistical analysis
##
## 
## Paula Sgarlatta
##
################################################################################

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


load(here::here("outputs", "temporal_alpha_kelp.RData") )
load(here::here("outputs", "spatial_alpha.RData") )
load(here::here("outputs", "shift3D_kelp_stats.RData") )
load(here::here("outputs", "shift3D_space_stats.RData") )

## Preparing data for stats

kelp_stats <- temporal_alpha_kelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)


spatial_stats <- spatial_alpha %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Habitat=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)


####### GLMM

### KELP ###

kelp_s <- glmmTMB(sp_richn ~ Year + (1|Site), data=kelp_stats, family = poisson())

summary(kelp_s)
Anova(kelp_s)#Significant

system.time(mod_1 <- drop1(kelp_s, test = "Chisq", all.cols=TRUE))
print(mod_1)

kelp_res_s <- simulateResiduals(kelp_s)
plot(kelp_res_s) #Good


### SPACE ###

space_s <- glmmTMB (sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = poisson()) 

summary(space_s)
Anova(space_s)#Significant

space_res_s <- simulateResiduals(space_s)
plot(space_res_s) #Good 

system.time(space_mod_s <- drop1(space_s, test = "Chisq", all.cols=TRUE))
print(space_mod_s)

### KELP - FRic ###

kelp_fric <- glmmTMB(fric ~ Year + (1|Site), data=kelp_stats, family=beta_family())

summary(kelp_fric)
Anova(kelp_fric)#Significant

kelp_res_fric <- simulateResiduals(kelp_fric)
plot(kelp_res_fric) # Good

system.time(kelp_mod_fric <- drop1(kelp_fric, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fric)

### SPATIAL - FRic ###

spatial_fric <- glmmTMB (fric ~ Habitat + (1|Site) , data=spatial_stats, family = beta_family())

summary(spatial_fric)
Anova(spatial_fric)#Significant

spatial_res_fric <- simulateResiduals(spatial_fric)
plot(spatial_res_fric) # Good

system.time(spatial_mod_fric <- drop1(spatial_fric, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fric)

### KELP - FDis ###

kelp_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=kelp_stats, family = beta_family())

summary(kelp_fdis)
Anova(kelp_fdis)#Significative

kelp_res_fdis <- simulateResiduals(kelp_fdis)
plot(kelp_res_fdis) # Good

system.time(kelp_mod_fdis <- drop1(kelp_fdis, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fdis)

### SPATIAL - Fdis ###

spatial_fdis <- glmmTMB (fdis ~ Habitat, data=spatial_stats, family = beta_family())

summary(spatial_fdis)
Anova(spatial_fdis)# No diferences - p=0.1097

spatial_res_fdis <- simulateResiduals(spatial_fdis)
plot(spatial_res_fdis) # Good

system.time(spatial_mod_fdis <- drop1(spatial_fdis, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fdis)

### KELP - Fide1 ###

kelp_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide1)
Anova(kelp_fide1) #No significant - p=0.2234

kelp_res_fide1 <- simulateResiduals(kelp_fide1)
plot(kelp_res_fide1) # Good

system.time(kelp_mod_fide1 <- drop1(kelp_fide1, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide1)

### SPATIAL - FIde1 ###

spatial_fide1 <- glmmTMB (fide_PC1 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

spatial_fide1 <- lm(fide_PC1 ~ Habitat, data=spatial_stats)

summary(spatial_fide1)
Anova(spatial_fide1)# Not significant - p=0.25

spatial_res_fide1 <- simulateResiduals(spatial_fide1)
plot(spatial_res_fide1) # Good

system.time(spatial_mod_fide1 <- drop1(spatial_fide1, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide1)

 
### KELP - Fide2 ###

kelp_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

#Not working

kelp_fide2 <- lmer (fide_PC2 ~ Year + (1|Site), data=kelp_stats)

summary(kelp_fide2)
Anova(kelp_fide2) #Significant

kelp_res_fide2 <- simulateResiduals(kelp_fide2)
plot(kelp_res_fide2) # Good

system.time(kelp_mod_fide2 <- drop1(kelp_fide2, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide2)

### SPATIAL - FIde2 ###

spatial_fide2 <- glmmTMB (fide_PC2 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide2)
Anova(spatial_fide2) #Not significant - p=0.831

spatial_res_fide2 <- simulateResiduals(spatial_fide2)
plot(spatial_res_fide2) # Good

system.time(spatial_mod_fide2 <- drop1(spatial_fide2, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide2)

### KELP - Fide3 ###

kelp_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide3)
Anova(kelp_fide3) #Not significant - p=0.1684

kelp_res_fide3 <- simulateResiduals(kelp_fide3)
plot(kelp_res_fide3) # Good

system.time(kelp_mod_fide3 <- drop1(kelp_fide3, test = "Chisq", all.cols=TRUE))
print(kelp_mod_fide3)

### SPATIAL - FIde3 ###

spatial_fide3 <- glmmTMB (fide_PC3 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide3)
Anova(spatial_fide3) #Not significant - p=0.25

spatial_res_fide3 <- simulateResiduals(spatial_fide3)
plot(spatial_res_fide3) # Good

system.time(spatial_mod_fide3 <- drop1(spatial_fide3, test = "Chisq", all.cols=TRUE))
print(spatial_mod_fide3)

## Shift FIde 


#Kelp

shift3D_kelp_stats$Year_cont <- as.numeric(shift3D_kelp_stats$Year)

kelp_shift <- lmer (shift3D ~ Year_cont + (1|Site1), data=shift3D_kelp_stats)
summary(kelp_shift)
Anova(kelp_shift)

system.time(shift_mod <- drop1(kelp_shift, test = "Chisq", all.cols=TRUE))
print(shift_mod)

#Space

space_shift <- glmmTMB (shift3D ~ Habitat , data=shift3D_space_stats)

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


load(here::here("outputs", "temporal_alpha_kelp_transect.RData") )
load(here::here("outputs", "spatial_alpha_transect.RData") )
load(here::here("outputs", "shift3D_kelp_stats.RData") )
load(here::here("outputs", "shift3D_space_stats.RData") )

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
Anova(kelp_s)#Significant

kelp_res_s <- simulateResiduals(kelp_s)
plot(kelp_res_s) #Good


### SPACE ###

space_s <- glmmTMB (sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = nbinom1()) 

space_lm_s <- lmer(sp_richn ~ Habitat + (1|Site), data = spatial_stats)

Anova(space_lm_s)

summary(space_s)
Anova(space_s)#Not significant

space_res_s <- simulateResiduals(space_lm_s)
plot(space_res_s) #Good 

#Post-hoc tests

pairs(emmeans(space_s, spec=~Habitat, type="response"))

### KELP - FRic ###

kelp_fric <- glmmTMB(fric ~ Year_cont + (1|Site), data=kelp_stats, family=beta_family())

summary(kelp_fric)
Anova(kelp_fric)#Significant

kelp_res_fric <- simulateResiduals(kelp_fric)
plot(kelp_res_fric) # Good

### SPATIAL - FRic ###

spatial_fric <- glmmTMB (fric ~ Habitat + (1|Site) , data=spatial_stats, family = beta_family())

summary(spatial_fric)
Anova(spatial_fric)#Significant

spatial_res_fric <- simulateResiduals(spatial_fric)
plot(spatial_res_fric) # Good

### KELP - FDis ###

kelp_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=kelp_stats, family = beta_family())

summary(kelp_fdis)
Anova(kelp_fdis)#Significative

kelp_res_fdis <- simulateResiduals(kelp_fdis)
plot(kelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(kelp_fdis, spec=~Year, type="response"))

#No differences, strange


### SPATIAL - Fdis ###

spatial_fdis <- glmmTMB (fdis ~ Habitat, data=spatial_stats, family = beta_family())

summary(spatial_fdis)
Anova(spatial_fdis)# No diferences - p=0.1097

spatial_res_fdis <- simulateResiduals(spatial_fdis)
plot(spatial_res_fdis) # Good

### KELP - Fide1 ###

kelp_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide1)
Anova(kelp_fide1) #No significant - p=0.2234

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

#Not working

kelp_fide2 <- lmer (fide_PC2 ~ Year + (1|Site), data=kelp_stats)

summary(kelp_fide2)
Anova(kelp_fide2) #Significant

kelp_res_fide2 <- simulateResiduals(kelp_fide2)
plot(kelp_res_fide2) # Good

### SPATIAL - FIde2 ###

spatial_fide2 <- glmmTMB (fide_PC2 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide2)
Anova(spatial_fide2) #Not significant - p=0.831

spatial_res_fide2 <- simulateResiduals(spatial_fide2)
plot(spatial_res_fide2) # Good

### KELP - Fide3 ###

kelp_fide3 <- glmmTMB(fide_PC3 ~ Year_cont + (1|Site), data=kelp_stats, family = gaussian())

summary(kelp_fide3)
Anova(kelp_fide3) #Not significant - p=0.1684

kelp_res_fide3 <- simulateResiduals(kelp_fide3)
plot(kelp_res_fide3) # Good

### SPATIAL - FIde3 ###

spatial_fide3 <- glmmTMB (fide_PC3 ~ Habitat + (1|Site),data=spatial_stats, family = gaussian())

summary(spatial_fide3)
Anova(spatial_fide3) #Not significant - p=0.25

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
