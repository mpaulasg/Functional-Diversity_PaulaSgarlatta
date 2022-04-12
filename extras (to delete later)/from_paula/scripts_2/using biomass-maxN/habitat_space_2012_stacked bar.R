################################################################################
##
## Script for analysing habitat data from the cross-shelf gradient
##
## 
## Code by Paula Sgarlatta
##
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# Load data

habitat <- read.csv(here::here("data", "using biomass-maxN", "habitat_solitaries_2012.csv"))

habitat_toplot <- habitat %>% 
  pivot_longer(cols= 4:10, names_to="Group") %>%
  filter(Site!="Flat_top", Site!="Look_at_me_now") %>% 
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Habitat,Transect, Group) %>%
  summarize(total=sum(value)) %>% 
  group_by(Site, Transect) %>% 
  mutate(percent_cover=(total*100/125)) %>% #125 points per transect
  select(-total)



####(FIGURE X-3) ----



habitat_summary_habitat <-habitat_toplot%>%
  group_by(Habitat, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            n = n(), mean = mean(percent_cover), se = sd(percent_cover/sqrt(n))) %>% 
  mutate(se) %>%
  group_by(Habitat) %>%
  arrange(desc(Group)) %>%
  mutate(
    pos = cumsum(Percent_cover_habitat),
    upper = pos + se/2,
    lower = pos - se/2
  ) %>%
  ungroup()

habitat_v3 <- ggplot(habitat_summary_habitat, aes(x=Habitat, y=Percent_cover_habitat, fill=Group)) +
  geom_bar(stat="identity", width = 0.7) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width=0.09, position = "identity")+
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14), color = "black"), axis.title = element_text(size= (18)),
        legend.key = element_rect(fill = "white"))+
  labs(y="Percent cover (%)", x="") + scale_fill_discrete(name = "", 
                                                          labels = c("Coral", "Ecklonia radiata", "Macroalgae", "Other invertebrates",
                                                                     "Rock & sand", "Sponges & tunicates", "Turf & CCA"))

## Error bars not working

ggsave(habitat_v3, file=here::here("outputs/", "habitat_stack.png"),
       height = 16, width = 24, unit = "cm" )