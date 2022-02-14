##GLMM

FBeta_Hill_temporal <- glmmTMB(q1 ~ Year + (1|Site), data = fx_beta_stats, family = beta_family())

summary(FBeta_Hill_temporal)

mod.res <- simulateResiduals(FBeta_Hill_temporal)
plot(mod.res) #Good!

## Pairwise comparisons

emmeans(FBeta_Hill_temporal, pairwise ~ Year)



### No kelp

FBeta_Hill_temporal_nk <- glmmTMB(q1 ~ Year + (1|Site), data = fx_beta_stats_nokelp, family = beta_family())

summary(FBeta_Hill_temporal_nk)

mod.res_nk <- simulateResiduals(FBeta_Hill_temporal_nk)
plot(mod.res_nk) #Check this one - not ok 


## Pairwise comparisons

emmeans(FBeta_Hill_temporal_nk, pairwise ~ Year)


#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_year_kelp$residuals)

# QQ-plot

qqPlot(s_aov_year_kelp$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_year_kelp)

####[PS] Assumptions met!

TukeyHSD(s_aov_year_kelp)  ## check


s_aov_year_no_kelp <- aov(sp_richn~year ,data=fd_temporal_no_kelp)

summary(s_aov_year_no_kelp)

#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_year_no_kelp$residuals)

# QQ-plot

qqPlot(s_aov_year_no_kelp$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_year_no_kelp)

####[PS] Assumptions met!

TukeyHSD(s_aov_year_no_kelp)  ## check

#Significative - p < 2.2e-16

#Check assumptions

qqnorm(residuals(s_aov_year))