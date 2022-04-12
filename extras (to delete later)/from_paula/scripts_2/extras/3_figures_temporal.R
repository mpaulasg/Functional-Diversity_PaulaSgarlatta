##GRAPH

#Richness (sites that use to have kelp) - TEMPORAL DATA

s_temp_kelp <- meta_data_kelp %>% # the names of the new data frame and the data frame to be summarised
  group_by(Year) %>%   # the grouping variable
  summarise(mean = mean(Richness),  # calculates the mean of each group
            sd = sd(Richness), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(Richness)/sqrt(n())) # calculates the standard error of each group


s_temporal_kelp <- ggplot(s_temp_kelp, aes(x=Year, y=mean))+
  geom_point()+ geom_smooth(color="mediumseagreen", method = "lm")

mytheme <-theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                axis.text = element_text(size = (14), colour = "black"), axis.title = element_text(size= (18)),
                legend.key = element_rect(fill = "white"))

print(s_temporal_kelp + mytheme +labs(y="Species richness", x="Years"))


#Functional Richness

fric_temp_kelp <- fd_temporal_kelp %>% # the names of the new data frame and the data frame to be summarised
  group_by(year) %>%   # the grouping variable
  summarise(fric_mean = mean(fric),  # calculates the mean of each group
            sd = sd(fric), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(fric)/sqrt(n())) # calculates the standard error of each group


fric_temporal_kelp <- ggplot(fric_temp_kelp, aes(x=year, y=fric_mean))+
  geom_point()+ geom_smooth(color="firebrick3", method = "lm")

mytheme <-theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                axis.text = element_text(size = (14), colour = "black"), axis.title = element_text(size= (18)),
                legend.key = element_rect(fill = "white"))

print(fric_temporal_kelp + mytheme +labs(y="Functional richness", x="Years"))

#Functional dispersion

fdis_temp_kelp <- fd_temporal_kelp %>% # the names of the new data frame and the data frame to be summarised
  group_by(year) %>%   # the grouping variable
  summarise(fdis_mean = mean(fdis),  # calculates the mean of each group
            sd = sd(fdis), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(fdis)/sqrt(n())) # calculates the standard error of each group


fdis_temporal <- ggplot(fdis_temp_kelp, aes(x=year, y=fdis_mean))+
  geom_point()+ stat_smooth(color="lightsalmon", formula = y ~ x + I(x^2),
                            method = "lm")

mytheme <-theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                axis.text = element_text(size = (14), colour = "black"), axis.title = element_text(size= (16)),
                legend.key = element_rect(fill = "white"))

print(fdis_temporal + mytheme +labs(y="Functional dispersion", x="Years"))

#Functional Specialisation

fspe_temp_kelp <- fd_temporal_kelp %>% # the names of the new data frame and the data frame to be summarised
  group_by(year) %>%   # the grouping variable
  summarise(fspe_mean = mean(fspe),  # calculates the mean of each group
            sd = sd(fspe), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(fspe)/sqrt(n())) # calculates the standard error of each group


fspe_temporal <- ggplot(fspe_temp_kelp, aes(x=year, y=fspe_mean))+
  geom_point()+ stat_smooth(color="lightsalmon", formula = y ~ x + I(x^2),
                            method = "lm")

mytheme <-theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                axis.text = element_text(size = (14), colour = "black"), axis.title = element_text(size= (16)),
                legend.key = element_rect(fill = "white"))

print(fspe_temporal + mytheme +labs(y="Functional Specialisation", x="Years"))

