################################################################################
##
## Script for FD statistical analysis
##
## 
## Paula Sgarlatta
##
################################################################################

######## Using site data

rm(list=ls()) # cleaning memory

##Load packages

library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(vegan)
library(performance)
library(here)

load(here::here("outputs", "temporal_alpha_nokelp.RData") )
load(here::here("outputs", "temporal_alpha_kelp.RData") )
load(here::here("outputs", "spatial_alpha.RData") )
load(here::here("outputs", "shift3D_nokelp_stats.RData") )
load(here::here("outputs", "shift3D_kelp_stats.RData") )
load(here::here("outputs", "shift3D_space_eucl.RData") )

## Preparing data for stats

nokelp_stats <- temporal_alpha_nokelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)

kelp_stats <- temporal_alpha_kelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)

spatial_stats <- spatial_alpha %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Habitat=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)


shift_nokelp <- shift3D_nokelp_stats

####### Linear mixed model

### NO KELP ###

nokelp_lm_s <- lmer(sp_richn ~ Year + (1|Site) ,data=nokelp_stats)

summary(nokelp_lm_s)
Anova(nokelp_lm_s)

mod.res <- simulateResiduals(nokelp_lm_s)
plot(mod.res) #Good

#Post-hoc tests

emmeans(nokelp_lm_s, pairwise~Year)

### KELP ###

kelp_lm_s <- lmer(sp_richn ~ Year + (1|Site),data=kelp_stats)

summary(kelp_lm_s)
Anova(kelp_lm_s)

kelp_res_s <- simulateResiduals(kelp_lm_s)
plot(kelp_res_s) #Good

#Post-hoc tests

pairs(emmeans(kelp_lm_s, spec=~Year, type="response"))

### SPACE ###

space_lm_s <- lm (sp_richn ~ Habitat, data=spatial_stats) 

#lmer not working here: Error: number of levels of each 
#grouping factor must be < number of observations (problems: Site)

summary(space_lm_s)
Anova(space_lm_s)

space_res_s <- simulateResiduals(space_lm_s)
plot(space_res_s) #Good? 

#Post-hoc tests

pairs(emmeans(space_lm_s, spec=~Habitat, type="response"))

### NO KELP - FRic ###

nokelp_lm_fric <- lmer(fric ~ Year + (1|Site), data=nokelp_stats)

summary(nokelp_lm_fric)
Anova(nokelp_lm_fric)

nokelp_res_fric <- simulateResiduals(nokelp_lm_fric)
plot(nokelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(nokelp_lm_fric, spec=~Year, type="response"))


### KELP - FRic ###

kelp_lm_fric <- lmer (fric ~ Year + (1|Site) ,data=kelp_stats)

summary(kelp_lm_fric)
Anova(kelp_lm_fric)

kelp_res_fric <- simulateResiduals(kelp_lm_fric)
plot(kelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(kelp_lm_fric, spec=~Year, type="response"))

### SPATIAL - FRic ###

spatial_lm_fric <- lm (fric ~ Habitat , data=spatial_stats)

summary(spatial_lm_fric)
Anova(spatial_lm_fric)

spatial_res_fric <- simulateResiduals(spatial_lm_fric)
plot(spatial_res_fric) # Good

#Post-hoc tests

pairs(emmeans(spatial_lm_fric, spec=~Habitat, type="response"))


### NO KELP - FDis ###

nokelp_lm_fdis <- lmer(fdis ~ Year + (1|Site), data=nokelp_stats)

summary(nokelp_lm_fdis)
Anova(nokelp_lm_fdis)

nokelp_res_fdis <- simulateResiduals(nokelp_lm_fdis)
plot(nokelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(nokelp_lm_fdis, spec=~Year, type="response"))

### KELP - FDis ###

kelp_lm_fdis <-lmer(fdis ~ Year + (1|Site),data=kelp_stats)

summary(kelp_lm_fdis)
Anova(kelp_lm_fdis)

kelp_res_fdis <- simulateResiduals(kelp_lm_fdis)
plot(kelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(kelp_lm_fdis, spec=~Year, type="response"))

#No differences, strange


### SPATIAL - Fdis ###

spatial_lm_fdis <- lm(fdis ~ Habitat ,data=spatial_stats)

summary(spatial_lm_fdis)
Anova(spatial_lm_fdis)

spatial_res_fdis <- simulateResiduals(spatial_lm_fdis)
plot(spatial_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(spatial_lm_fdis, spec=~Habitat, type="response")) # Not significant



### NO KELP - FIde1 ###

nokelp_lm_fide1 <- lmer(fide_PC1 ~ Year + (1|Site), data=nokelp_stats)

summary(nokelp_lm_fide1)
Anova(nokelp_lm_fide1)

nokelp_res_fide <- simulateResiduals(nokelp_lm_fide1)
plot(nokelp_res_fide) # Good

#Post-hoc tests

pairs(emmeans(nokelp_lm_fide1, spec=~Year, type="response"))

### KELP - Fide1 ###

kelp_lm_fide1 <-lmer(fide_PC1 ~ Year + (1|Site),data=kelp_stats)

summary(kelp_lm_fide1)
Anova(kelp_lm_fide1)

kelp_res_fide1 <- simulateResiduals(kelp_lm_fide1)
plot(kelp_res_fide1) # Good

### SPATIAL - FIde1 ###

spatial_lm_fide1 <- lm(fide_PC1 ~ Habitat ,data=spatial_stats)

summary(spatial_lm_fide1)
Anova(spatial_lm_fide1)

spatial_res_fide1 <- simulateResiduals(spatial_lm_fide1)
plot(spatial_res_fide1) # Good

 # Not significant


## Shift FIde 

## No kelp

nokelp_perma <- adonis(shift3D ~ Year, data=shift_nokelp)
print(nokelp_perma)

#Kelp

kelp_perma <- adonis(shift3D ~ Year, data=shift3D_kelp_stats)
print(kelp_perma)

#Space

Habitat <- spatial_alpha %>% 
  rownames_to_column("Site") %>% 
  mutate(Habitat=sub(".*_", "", Site)) %>% 
  select(Habitat, Site)
  
  
space_perma <- adonis(shift3D_space_eucl ~ Habitat, data=Habitat)

print(space_perma)

########################### Using transect data

rm(list=ls()) # cleaning memory

##Load packages

library(tidyverse)
library(dplyr)
library(car)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(vegan)
library(performance)
library(here)

load(here::here("outputs", "temporal_alpha_nokelp_all.RData") )
load(here::here("outputs", "temporal_alpha_kelp_all.RData") )
load(here::here("outputs", "spatial_alpha_all.RData") )
load(here::here("outputs", "shift3D_nokelp_eucl_1.RData") )
load(here::here("outputs", "shift3D_kelp_eucl_1 .RData") )
load(here::here("outputs", "shift3D_space_eucl_1 .RData") )
load(here::here("data", "kelp_metadata_all.RData") )
load(here::here("data", "nokelp_metadata_all.RData") )
load(here::here("data", "spatial_metadata_all.RData") )
load(here::here("outputs", "temporal_alpha_nokelp_all.RData") )


## Preparing data for stats

nokelp_stats <- temporal_alpha_nokelp_all %>%
  rownames_to_column(var = "Site") %>% 
  left_join(nokelp_metadata_all, by=c("Site"="Code")) %>% 
  mutate(Sites=Site.y, .before="sp_richn") %>% 
  select(-Site, -Site.y)

kelp_stats <- temporal_alpha_kelp_all %>%
  rownames_to_column(var = "Site") %>% 
  left_join(kelp_metadata_all, by=c("Site"="Code")) %>% 
  mutate(Sites=Site.y, .before="sp_richn") %>% 
  select(-Site, -Site.y, -YearID)

spatial_stats <- spatial_alpha_all %>%
  rownames_to_column(var = "Site") %>% 
  left_join(spatial_metadata_all, by=c("Site"="Code")) %>% 
  mutate(Sites=Site.y, .before="sp_richn") %>% 
  select(-Site, -Site.y, -HabitatID)


# Generalized Linear Mixed Modesl (GLMM) for species richnes, FRic, FDis and FIde

### NO KELP ###

nokelp_glmm_s <- glmmTMB(sp_richn ~ Year + (1|Sites:Replicate) ,data=nokelp_stats, family=poisson())


summary(nokelp_glmm_s)
Anova(nokelp_glmm_s)

mod.res <- simulateResiduals(nokelp_glmm_s)
plot(mod.res) #Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_s, spec=~Year, type="response"))

#Not working

#Check overdispertion

check_overdispersion(nokelp_glmm_s) #No overdispersion (dispertion ratio = 0.786)

#Not working

### KELP ###

kelp_glmm_s <- glmmTMB(sp_richn ~ Year + (1|Sites) ,data=kelp_stats, family=poisson())

#Not working with replicates 

summary(kelp_glmm_s)
Anova(kelp_glmm_s)

kelp_res_s <- simulateResiduals(kelp_glmm_s)
plot(kelp_res_s) #Good

check_overdispersion(kelp_glmm_s) #No overdispersion (dispertion ratio = 0.847)

#Post-hoc tests

pairs(emmeans(kelp_glmm_s, spec=~Year, type="response"))

#Not working 

### SPACE ###

space_glmm_s <- glmmTMB(sp_richn ~ Habitat + (1|Sites/Replicate), data=spatial_stats, family = poisson()) 

summary(space_glmm_s)
Anova(space_glmm_s)

space_res_s <- simulateResiduals(space_glmm_s)
plot(space_res_s) #Good?

check_overdispersion(space_glmm_s) #No overdispersion (dispertion ratio = 0.537)

#Post-hoc tests

pairs(emmeans(space_glmm_s, spec=~Habitat, type="response"))

#### This one is working perfectly

### NO KELP - FRic ###

nokelp_glmm_fric <- glmmTMB(fric ~ Year + (1|Sites:Replicate), data=nokelp_stats, family = beta_family())

summary(nokelp_glmm_fric)
Anova(nokelp_glmm_fric)

nokelp_res_fric <- simulateResiduals(nokelp_glmm_fric)
plot(nokelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fric, spec=~Year, type="response"))

#Not working

### KELP - FRic ###

kelp_glmm_fric <- glmmTMB(fric ~ Year + (1|Sites:Replicate) ,data=kelp_stats, family = beta_family())

summary(kelp_glmm_fric)
Anova(kelp_glmm_fric)

kelp_res_fric <- simulateResiduals(kelp_glmm_fric)
plot(kelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fric, spec=~Year, type="response"))

#Not working

### SPATIAL - FRic ###

spatial_glmm_fric <- glmmTMB(fric ~ Habitat + (1|Sites:Replicate), data=spatial_stats, family = beta_family())

summary(spatial_glmm_fric)
Anova(spatial_glmm_fric)

spatial_res_fric <- simulateResiduals(spatial_glmm_fric)
plot(spatial_res_fric) # Good?

#Post-hoc tests

pairs(emmeans(spatial_glmm_fric, spec=~Habitat, type="response"))


### NO KELP - FDis ###

nokelp_glmm_fdis <- glmmTMB(fdis ~ Year + (1|Sites:Replicate), data=nokelp_stats, family = beta_family())

summary(nokelp_glmm_fdis)
Anova(nokelp_glmm_fdis)

nokelp_res_fdis <- simulateResiduals(nokelp_glmm_fdis)
plot(nokelp_res_fdis) # Not good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fdis, spec=~Year, type="response"))

### Not working

### KELP - FDis ###

kelp_glmm_fdis <- glmmTMB(fdis ~ Year + (1|Sites:Replicate) ,data=kelp_stats, family = beta_family())

summary(kelp_glmm_fdis)
Anova(kelp_glmm_fdis)

kelp_res_fdis <- simulateResiduals(kelp_glmm_fdis)
plot(kelp_res_fdis) # Not good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fdis, spec=~Year, type="response"))

#Not working


### SPATIAL - Fdis ###

spatial_glmm_fdis <- glmmTMB(fdis ~ Habitat + (1|Sites:Replicate) ,data=spatial_stats, family = beta_family())

summary(spatial_glmm_fdis)
Anova(spatial_glmm_fdis)

spatial_res_fdis <- simulateResiduals(spatial_glmm_fdis)
plot(spatial_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(spatial_glmm_fdis, spec=~Habitat, type="response")) # Not significant


### NO KELP - FIde 1 ###

nokelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Sites), data=nokelp_stats, family = gaussian()) #Not sure of the family

summary(nokelp_glmm_fide1)
Anova(nokelp_glmm_fide1)

nokelp_res_fide1 <- simulateResiduals(nokelp_glmm_fide1)
plot(nokelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide1, spec=~Year, type="response"))

#Not working

### KELP - FIde 1 ###

kelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide1)
Anova(kelp_glmm_fide1)

kelp_res_fide1 <- simulateResiduals(kelp_glmm_fide1)
plot(kelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide1, spec=~Year, type="response"))

## Not working


### SPATIAL - FIde 1 ###

spatial_glmm_fide1 <- glmmTMB(fide_PC1 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide1)
Anova(spatial_glmm_fide1)

spatial_res_fide1 <- simulateResiduals(spatial_glmm_fide1)
plot(spatial_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide1, spec=~Habitat, type="response")) # Not significant


### NO KELP - FIde 2 ###

nokelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Sites), data=nokelp_stats, family = gaussian()) 

summary(nokelp_glmm_fide2)
Anova(nokelp_glmm_fide2)

nokelp_res_fide2 <- simulateResiduals(nokelp_glmm_fide2)
plot(nokelp_res_fide2) # Not good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide2, spec=~Year, type="response"))

#Not working

### KELP - FIde 2 ###

kelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide2)
Anova(kelp_glmm_fide2)

kelp_res_fide2 <- simulateResiduals(kelp_glmm_fide2)
plot(kelp_res_fide2) # Not good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide2, spec=~Year, type="response"))

## Not working

### SPATIAL - FIde 2 ###

spatial_glmm_fide2 <- glmmTMB(fide_PC2 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide2)
Anova(spatial_glmm_fide2)

spatial_res_fide2 <- simulateResiduals(spatial_glmm_fide2)
plot(spatial_res_fide2) # Good

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide2, spec=~Habitat, type="response")) # Not significant

### NO KELP - FIde 3 ###

nokelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Sites:Replicate), data=nokelp_stats, family = gaussian()) #Not sure of the family

#Not working
summary(nokelp_glmm_fide3)
Anova(nokelp_glmm_fide3)

nokelp_res_fide3 <- simulateResiduals(nokelp_glmm_fide3)
plot(nokelp_res_fide3) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide3, spec=~Year, type="response"))


### KELP - FIde 3 ###

kelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide3)
Anova(kelp_glmm_fide3) # Significant

kelp_res_fide3 <- simulateResiduals(kelp_glmm_fide3)
plot(kelp_res_fide3) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide3, spec=~Year, type="response"))

#Not working

### SPATIAL - FIde 3 ###

spatial_glmm_fide3 <- glmmTMB(fide_PC3 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide3)
Anova(spatial_glmm_fide3)

spatial_res_fide3 <- simulateResiduals(spatial_glmm_fide3)
plot(spatial_res_fide3) # Good


##################################### Try linear model

### NO KELP ###

nokelp_lm_s <- lmer(sp_richn ~ Year + (1|Sites) ,data=nokelp_stats)

summary(nokelp_lm_s)
Anova(nokelp_lm_s)

mod.res <- simulateResiduals(nokelp_lm_s)
plot(mod.res) #Good

#Post-hoc tests

emmeans(nokelp_lm_s, pairwise~Year)

#Not working

### KELP ###

kelp_lm_s <- lm(sp_richn ~ Year,data=kelp_stats)

summary(kelp_lm_s)
Anova(kelp_lm_s)

kelp_res_s <- simulateResiduals(kelp_lm_s)
plot(kelp_res_s) #Good

#Post-hoc tests

pairs(emmeans(kelp_lm_s, spec=~Year, type="response"))

#Not working 

### SPACE ###

space_lm_s <- lm(sp_richn ~ Habitat, data=spatial_stats) 

summary(space_lm_s)
Anova(space_lm_s)

space_res_s <- simulateResiduals(space_lm_s)
plot(space_res_s) #Good? Check for Levene

#Post-hoc tests

pairs(emmeans(space_lm_s, spec=~Habitat, type="response"))

#### This one is working perfectly

### NO KELP - FRic ###

nokelp_lm_fric <- lm(fric ~ Year, data=nokelp_stats)

summary(nokelp_lm_fric)
Anova(nokelp_lm_fric)

nokelp_res_fric <- simulateResiduals(nokelp_lm_fric)
plot(nokelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(nokelp_lm_fric, spec=~Year, type="response"))

#Not working

### KELP - FRic ###

kelp_lm_fric <- lm(fric ~ Year ,data=kelp_stats)

summary(kelp_lm_fric)
Anova(kelp_lm_fric)

kelp_res_fric <- simulateResiduals(kelp_lm_fric)
plot(kelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(kelp_lm_fric, spec=~Year, type="response"))

#Not working

### SPATIAL - FRic ###

spatial_lm_fric <- lm(fric ~ Habitat , data=spatial_stats)

summary(spatial_lm_fric)
Anova(spatial_lm_fric)

spatial_res_fric <- simulateResiduals(spatial_lm_fric)
plot(spatial_res_fric) # Good

#Post-hoc tests

pairs(emmeans(spatial_lm_fric, spec=~Habitat, type="response"))


### NO KELP - FDis ###

nokelp_lm_fdis <- lm(fdis ~ Year, data=nokelp_stats)

summary(nokelp_lm_fdis)
Anova(nokelp_lm_fdis)

nokelp_res_fdis <- simulateResiduals(nokelp_lm_fdis)
plot(nokelp_res_fdis) # Not good

#Post-hoc tests

pairs(emmeans(nokelp_lm_fdis, spec=~Year, type="response"))

### Not working

### KELP - FDis ###

kelp_lm_fdis <-lm(fdis ~ Year,data=kelp_stats)

summary(kelp_lm_fdis)
Anova(kelp_lm_fdis)

kelp_res_fdis <- simulateResiduals(kelp_lm_fdis)
plot(kelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(kelp_lm_fdis, spec=~Year, type="response"))

#Not working


### SPATIAL - Fdis ###

spatial_lm_fdis <- lm(fdis ~ Habitat ,data=spatial_stats)

summary(spatial_lm_fdis)
Anova(spatial_lm_fdis)

spatial_res_fdis <- simulateResiduals(spatial_lm_fdis)
plot(spatial_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(spatial_lm_fdis, spec=~Habitat, type="response")) # Not significant


### Keep trying from here

### NO KELP - FIde 1 ###

nokelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Sites), data=nokelp_stats, family = gaussian()) #Not sure of the family

summary(nokelp_glmm_fide1)
Anova(nokelp_glmm_fide1)

nokelp_res_fide1 <- simulateResiduals(nokelp_glmm_fide1)
plot(nokelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide1, spec=~Year, type="response"))

#Not working

### KELP - FIde 1 ###

kelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide1)
Anova(kelp_glmm_fide1)

kelp_res_fide1 <- simulateResiduals(kelp_glmm_fide1)
plot(kelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide1, spec=~Year, type="response"))

## Not working


### SPATIAL - FIde 1 ###

spatial_glmm_fide1 <- glmmTMB(fide_PC1 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide1)
Anova(spatial_glmm_fide1)

spatial_res_fide1 <- simulateResiduals(spatial_glmm_fide1)
plot(spatial_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide1, spec=~Habitat, type="response")) # Not significant


### NO KELP - FIde 2 ###

nokelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Sites), data=nokelp_stats, family = gaussian()) 

summary(nokelp_glmm_fide2)
Anova(nokelp_glmm_fide2)

nokelp_res_fide2 <- simulateResiduals(nokelp_glmm_fide2)
plot(nokelp_res_fide2) # Not good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide2, spec=~Year, type="response"))

#Not working

### KELP - FIde 2 ###

kelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide2)
Anova(kelp_glmm_fide2)

kelp_res_fide2 <- simulateResiduals(kelp_glmm_fide2)
plot(kelp_res_fide2) # Not good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide2, spec=~Year, type="response"))

## Not working

### SPATIAL - FIde 2 ###

spatial_glmm_fide2 <- glmmTMB(fide_PC2 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide2)
Anova(spatial_glmm_fide2)

spatial_res_fide2 <- simulateResiduals(spatial_glmm_fide2)
plot(spatial_res_fide2) # Good

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide2, spec=~Habitat, type="response")) # Not significant

### NO KELP - FIde 3 ###

nokelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Sites:Replicate), data=nokelp_stats, family = gaussian()) #Not sure of the family

#Not working
summary(nokelp_glmm_fide3)
Anova(nokelp_glmm_fide3)

nokelp_res_fide3 <- simulateResiduals(nokelp_glmm_fide3)
plot(nokelp_res_fide3) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide3, spec=~Year, type="response"))


### KELP - FIde 3 ###

kelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Sites) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide3)
Anova(kelp_glmm_fide3) # Significant

kelp_res_fide3 <- simulateResiduals(kelp_glmm_fide3)
plot(kelp_res_fide3) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide3, spec=~Year, type="response"))

#Not working

### SPATIAL - FIde 3 ###

spatial_glmm_fide3 <- glmmTMB(fide_PC3 ~ Habitat + (1|Sites) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide3)
Anova(spatial_glmm_fide3)

spatial_res_fide3 <- simulateResiduals(spatial_glmm_fide3)
plot(spatial_res_fide3) # Good



##### Distance-based redundancy analysis with jaccard dissimilarity for shift in FIde


dbrda_nokelp  <- capscale(shift3D_nokelp_stats ~ Year, data = nokelp_metadata)
anova(dbrda_nokelp, by = "terms", permutations = 9999)

# Error in dfun(X, distance) : 
#   missing values are not allowed with argument 'na.rm = FALSE'



