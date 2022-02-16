rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

kelp_data <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") )
kelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv"))
traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_kelp.csv"))

# merging metadata and data 
kelp <- left_join(kelp_metadata, kelp_data, by="Code" )

# computing mean maxN per year ----
kelp_maxN <- kelp %>% 
  pivot_longer(cols= contains("_"), names_to="Species") %>% 
mutate(value=as.numeric(value)) 

# merging traits and data 
kelp_traits <- left_join(traits, kelp_maxN, by="Species" )

######### Aggregation 

kelp_agg <- kelp_traits %>% 
  group_by(Year, Agg) %>%
  summarize(
    n = n(),
    agg_mean = mean(value),
    agg_sd = sd(value)
  ) 

# color code for agg

agg <- c(Solitary= "seagreen4", Group= "seagreen2", Pair="lightsalmon1")


agg_plot <- ggplot(kelp_agg, aes (x= Year,y = agg_mean, color=Agg)) +
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=agg) + 
  scale_fill_manual(values=agg) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(agg_plot, file=here::here("outputs/", "agg_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Size 

kelp_size <- kelp_traits %>% 
  group_by(Year, Size) %>%
  summarize(
    n = n(),
    size_mean = mean(value),
    size_sd = sd(value)
  ) 

# color code for size

size <- c(S1= "green", S2= "red", S3="blue", S4="orange", S5="yellow", S6="black")


size_plot <- ggplot(kelp_size, aes (x= Year,y = size_mean, color=Size)) +
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=size) + 
  scale_fill_manual(values=size) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(size_plot, file=here::here("outputs/", "size_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Position 

kelp_position <- kelp_traits %>% 
  group_by(Year, Position) %>%
  summarize(
    n = n(),
    pos_mean = mean(value),
    pos_sd = sd(value)
  ) 

# color code for Position

pos <- c(BenthoP= "green", Benthic= "red", Pelagic="blue")


position_plot <- ggplot(kelp_position, aes (x= Year,y = pos_mean, color=Position)) +
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=pos) + 
  scale_fill_manual(values=pos) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(position_plot, file=here::here("outputs/", "position_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Diet

kelp_diet <- kelp_traits %>% 
  group_by(Year, Diet) %>%
  summarize(
    n = n(),
    diet_mean = mean(value),
    diet_sd = sd(value)
  ) 

## saving as csv ----
write.csv(kelp_diet, file=here::here("data", "raw_data", "kelp_diet.csv"), row.names = FALSE )


diet_plot <- ggplot(kelp_diet, aes (x= Year,y = diet_mean, color=Diet)) +
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  #scale_color_manual(values=diet) + 
  #scale_fill_manual(values=diet) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(diet_plot, file=here::here("outputs/", "diet_plot.png"),
       height = 22, width = 25, unit = "cm" )







########### Use this one for traits - CHECK


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

kelp_data <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv") )
kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv"))
traits <- read.csv(here::here("from_paula", "kelp_traits.csv"))

# merging metadata and data 
kelp <- left_join(kelp_metadata, kelp_data, by=c("Code"="Site"))

# computing mean maxN per year ----
kelp_maxN <- kelp %>% 
  pivot_longer(cols= contains("_"), names_to="Species") %>% 
  mutate(value=as.numeric(value)) 

# merging traits and data 
kelp_traits <- left_join(traits, kelp_maxN, by="Species" )

######### Aggregation 

kelp_agg <- kelp_traits %>% 
  group_by(Year, Agg, Code, Site) %>%
  summarize(sum.MaxN=sum(value))
  

# color code for agg

agg <- c(Solitary= "seagreen4", Group= "seagreen2", Pair="lightsalmon1")


agg_plot <- ggplot(kelp_agg, aes (x= Year,y = sum.MaxN, color=Agg)) +
  stat_summary(geom = "point", show.legend = F)+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=agg) + 
  scale_fill_manual(values=agg) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(agg_plot, file=here::here("outputs/", "agg_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Size 

kelp_size <- kelp_traits %>% 
  group_by(Year, Size, Code, Site) %>%
  summarize(sum.MaxN=sum(value))

# color code for size

size <- c(S1= "green", S2= "red", S3="blue", S4="orange", S5="yellow", S6="black")


size_plot <- ggplot(kelp_size, aes (x= Year,y = sum.MaxN, color=Size)) +
  stat_summary(geom = "point", show.legend = F)+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=size) + 
  scale_fill_manual(values=size) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(size_plot, file=here::here("outputs/", "size_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Position 

kelp_position <- kelp_traits %>% 
  group_by(Year, Position, Code, Site) %>%
  summarize(sum.MaxN=sum(value))
    
   # color code for Position

pos <- c(BenthoP= "green", Benthic= "red", Pelagic="blue")


position_plot <- ggplot(kelp_position, aes (x= Year,y = sum.MaxN, color=Position)) +
  stat_summary(geom = "point", show.legend = F)+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=pos) + 
  scale_fill_manual(values=pos) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(position_plot, file=here::here("outputs/", "position_plot.png"),
       height = 22, width = 25, unit = "cm" )

######### Diet

kelp_diet <- kelp_traits %>% 
  group_by(Year, Diet, Site, Code) %>%
  summarize( sum.MaxN=sum(value))
   
    ## saving as csv ----
write.csv(kelp_diet, file=here::here("data", "raw_data", "kelp_diet.csv"), row.names = FALSE )


diet_plot <- ggplot(kelp_diet, aes (x= Year,y = sum.MaxN, color=Diet)) +
  stat_summary(geom = "point", show.legend = F)+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  #scale_color_manual(values=diet) + 
  #scale_fill_manual(values=diet) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(diet_plot, file=here::here("outputs/", "diet_plot.png"),
       height = 22, width = 25, unit = "cm" )


########## Kmax

kelp_kmax <- kelp_traits %>% 
  mutate(Kmax_plot=case_when(
    Kmax <0.1 ~ "<0.1",
    Kmax >0.1 & Kmax <0.4 ~ "0.1-0.4",
    Kmax >0.4 ~ ">0.4"
  )) %>% 
  group_by(Year, Kmax_plot, Site, Code) %>%
    summarize( sum.MaxN=sum(value))



Kmax_plot <- ggplot(kelp_kmax, aes (x= Year,y = sum.MaxN, color=Kmax_plot)) +
  stat_summary(geom = "point", show.legend = F)+
  geom_smooth(method = "lm")+
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  #scale_color_manual(values=diet) + 
  #scale_fill_manual(values=diet) + 
  labs(x="", y="Total maxN") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

ggsave(Kmax_plot, file=here::here("outputs/", "kmax_plot.png"),
       height = 22, width = 25, unit = "cm" )


