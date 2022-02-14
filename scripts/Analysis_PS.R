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
library(DHARMa)
library(emmeans)
library(vegan)
library(here)

load(here::here("outputs", "temporal_alpha_nokelp.RData") )
load(here::here("outputs", "temporal_alpha_kelp.RData") )
load(here::here("outputs", "spatial_alpha.RData") )
load(here::here("outputs", "shift3D_nokelp_eucl_1.RData") )
load(here::here("outputs", "shift3D_kelp_eucl_1 .RData") )
load(here::here("outputs", "shift3D_space_eucl_1 .RData") )
load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "nokelp_metadata.RData") )
load(here::here("data", "spatial_metadata.RData") )
load(here::here("outputs", "temporal_alpha_nokelp.RData") )


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

#Distances for using dbRDA

fide3_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC1,     fide_PC2, fide_PC3) 

shift3D_nokelp_eucl <- dist(fide3_nokelp, method = "euclidean")

df_shift3D_nokelp_2 <- shift3D_nokelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(x1, x2, shift3D)

shift3D_nokelp_stats <- df_shift3D_nokelp_2 %>% 
  spread(x2, shift3D) %>% 
  column_to_rownames("x1") %>% 
  as.matrix()

nokelp_metadata <- nokelp_metadata %>% 
  `rownames<-`(.$Code)

df_shift3D_kelp_2 <- shift3D_kelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(x1, x2, shift3D)

shift3D_kelp_stats <- df_shift3D_kelp_2 %>% 
  spread(x2, shift3D) %>% 
  column_to_rownames("x1")

shift3D_space_stats <- shift3D_space_eucl_1 %>% 
  select(x1, x2, shift3D) %>% 
  spread(x2, shift3D) %>% 
  column_to_rownames("x1")


# Generalized Linear Mixed Modesl (GLMM) for species richnes, FRic, FDis and FIde

### NO KELP ###

nokelp_glmm_s <- glmmTMB(sp_richn ~ Year + (1|Site) ,data=nokelp_stats, family=poisson())

summary(nokelp_glmm_s)
Anova(nokelp_glmm_s)

mod.res <- simulateResiduals(nokelp_glmm_s)
plot(mod.res) #Not good

nokelp_glmm_s <- glmmTMB(log(sp_richn) ~ Year + (1|Site) ,data=nokelp_stats, family=poisson())

#Still bad with transformation, what is the correct one?


### KELP ###

kelp_glmm_s <- glmmTMB(sp_richn ~ Year + (1|Site) ,data=kelp_stats, family=poisson())

summary(kelp_glmm_s)
Anova(kelp_glmm_s)

kelp_res_s <- simulateResiduals(kelp_glmm_s)
plot(kelp_res_s) #Deavition significant

#Post-hoc tests

pairs(emmeans(kelp_glmm_s, spec=~Year, type="response"))


### SPACE ###

space_glmm_s <- glmmTMB(sp_richn ~ Habitat + (1|Site), data=spatial_stats, family = poisson()) 

summary(space_glmm_s)
Anova(space_glmm_s)

space_res_s <- simulateResiduals(space_glmm_s)
plot(space_res_s)
#Should I use the data from transects?

#Post-hoc tests

pairs(emmeans(space_glmm_s, spec=~Habitat, type="response"))



### NO KELP - FRic ###

nokelp_glmm_fric <- glmmTMB(fric ~ Year + (1|Site), data=nokelp_stats, family = beta_family())

summary(nokelp_glmm_fric)
Anova(nokelp_glmm_fric)

nokelp_res_fric <- simulateResiduals(nokelp_glmm_fric)
plot(nokelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fric, spec=~Year, type="response"))

### KELP - FRic ###

kelp_glmm_fric <- glmmTMB(fric ~ Year + (1|Site) ,data=kelp_stats, family = beta_family())

summary(kelp_glmm_fric)
Anova(kelp_glmm_fric)

kelp_res_fric <- simulateResiduals(kelp_glmm_fric)
plot(kelp_res_fric) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fric, spec=~Year, type="response"))


### SPATIAL - FRic ###

spatial_glmm_fric <- glmmTMB(fric ~ Habitat + (1|Site) ,data=spatial_stats, family = beta_family())

summary(spatial_glmm_fric)
Anova(spatial_glmm_fric)

spatial_res_fric <- simulateResiduals(spatial_glmm_fric)
plot(spatial_res_fric) # Good?

#Post-hoc tests

pairs(emmeans(spatial_glmm_fric, spec=~Habitat, type="response"))

### NO KELP - FDis ###

nokelp_glmm_fdis <- glmmTMB(fdis ~ Year + (1|Site), data=nokelp_stats, family = beta_family())

summary(nokelp_glmm_fdis)
Anova(nokelp_glmm_fdis)

nokelp_res_fdis <- simulateResiduals(nokelp_glmm_fdis)
plot(nokelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fdis, spec=~Year, type="response"))

### KELP - FDis ###

kelp_glmm_fdis <- glmmTMB(fdis ~ Year + (1|Site) ,data=kelp_stats, family = beta_family())

summary(kelp_glmm_fdis)
Anova(kelp_glmm_fdis)

kelp_res_fdis <- simulateResiduals(kelp_glmm_fdis)
plot(kelp_res_fdis) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fdis, spec=~Year, type="response"))

## How it can be that glmm is significative, but then the variables are not significative in the post-hoc?


### SPATIAL - Fdis ###

spatial_glmm_fdis <- glmmTMB(fdis ~ Habitat + (1|Site) ,data=spatial_stats, family = beta_family())

summary(spatial_glmm_fdis)
Anova(spatial_glmm_fdis)

spatial_res_fdis <- simulateResiduals(spatial_glmm_fdis)
plot(spatial_res_fdis) # Good?

#Post-hoc tests

pairs(emmeans(spatial_glmm_fdis, spec=~Habitat, type="response")) # Not significant


### NO KELP - FIde 1 ###

nokelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site), data=nokelp_stats, family = gaussian()) #Not sure of the family

summary(nokelp_glmm_fide1)
Anova(nokelp_glmm_fide1)

nokelp_res_fide1 <- simulateResiduals(nokelp_glmm_fide1)
plot(nokelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide1, spec=~Year, type="response"))

### KELP - FIde 1 ###

kelp_glmm_fide1 <- glmmTMB(fide_PC1 ~ Year + (1|Site) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide1)
Anova(kelp_glmm_fide1)

kelp_res_fide1 <- simulateResiduals(kelp_glmm_fide1)
plot(kelp_res_fide1) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide1, spec=~Year, type="response"))

## How it can be that glmm is significative, but then the variables are not significative in the post-hoc?


### SPATIAL - FIde 1 ###

spatial_glmm_fide1 <- glmmTMB(fide_PC1 ~ Habitat + (1|Site) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide1)
Anova(spatial_glmm_fide1)

spatial_res_fide1 <- simulateResiduals(spatial_glmm_fide1)
plot(spatial_res_fide1) # Good?

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide1, spec=~Habitat, type="response")) # Not significant




### NO KELP - FIde 2 ###

nokelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site), data=nokelp_stats, family = gaussian()) #Not sure of the family

summary(nokelp_glmm_fide2)
Anova(nokelp_glmm_fide2)

nokelp_res_fide2 <- simulateResiduals(nokelp_glmm_fide2)
plot(nokelp_res_fide2) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide2, spec=~Year, type="response"))

### KELP - FIde 2 ###

kelp_glmm_fide2 <- glmmTMB(fide_PC2 ~ Year + (1|Site) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide2)
Anova(kelp_glmm_fide2)

kelp_res_fide2 <- simulateResiduals(kelp_glmm_fide2)
plot(kelp_res_fide2) # Good

#Post-hoc tests

pairs(emmeans(kelp_glmm_fide2, spec=~Year, type="response"))

## How it can be that glmm is significative, but then the variables are not significative in the post-hoc?


### SPATIAL - FIde 2 ###

spatial_glmm_fide2 <- glmmTMB(fide_PC2 ~ Habitat + (1|Site) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide2)
Anova(spatial_glmm_fide2)

spatial_res_fide2 <- simulateResiduals(spatial_glmm_fide2)
plot(spatial_res_fide2) # Good?

#Post-hoc tests

pairs(emmeans(spatial_glmm_fide2, spec=~Habitat, type="response")) # Not significant




### NO KELP - FIde 3 ###

nokelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Site), data=nokelp_stats, family = gaussian()) #Not sure of the family

summary(nokelp_glmm_fide3)
Anova(nokelp_glmm_fide3)

nokelp_res_fide3 <- simulateResiduals(nokelp_glmm_fide3)
plot(nokelp_res_fide3) # Good

#Post-hoc tests

pairs(emmeans(nokelp_glmm_fide3, spec=~Year, type="response"))

## How it can be that glmm is significative, but then the variables are not significative in the post-hoc?

### KELP - FIde 3 ###

kelp_glmm_fide3 <- glmmTMB(fide_PC3 ~ Year + (1|Site) ,data=kelp_stats, family = gaussian())

summary(kelp_glmm_fide3)
Anova(kelp_glmm_fide3) # Not significant

kelp_res_fide3 <- simulateResiduals(kelp_glmm_fide3)
plot(kelp_res_fide3) # Good

### SPATIAL - FIde 3 ###

spatial_glmm_fide3 <- glmmTMB(fide_PC3 ~ Habitat + (1|Site) ,data=spatial_stats, family = gaussian())

summary(spatial_glmm_fide3)
Anova(spatial_glmm_fide3)

spatial_res_fide3 <- simulateResiduals(spatial_glmm_fide3)
plot(spatial_res_fide3) # Good?



##### Distance-based redundancy analysis with jaccard dissimilarity for shift in FIde


dbrda_nokelp  <- capscale(shift3D_nokelp_stats ~ Year, data = nokelp_metadata)
anova(dbrda_nokelp, by = "terms", permutations = 9999)

# Error in dfun(X, distance) : 
#   missing values are not allowed with argument 'na.rm = FALSE'



