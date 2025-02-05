##Calculating mean and SD of stream morphology characteristics for each stream

library(tidyverse)
library(lubridate)
library(ggplot2)


##Import csv file containing cleaned stream conditions data  
data <- read.csv("data/physical_stream_conditions.csv", na.strings=c("","NA"))

data$newdate <-strptime(as.character(data$Date), "%d/%m/%Y")
data$txtdate <- format(data$newdate, "%y-%m-%d")

##stream morph

sm <- data %>%
  select(newdate, txtdate, Site_Code, Water_depth_cm, Wetted_width_cm, Hydraulic_head_cm) %>%
  filter(!is.na(Water_depth_cm)) %>%
  filter(!is.na(Wetted_width_cm)) %>%
  filter(Hydraulic_head_cm != 100) ##this filters out the one sample from AT that is definitely an error in recording 

sm$Wetted_width_cm <- as.numeric(sm$Wetted_width_cm)

##Boxplots 
ggplot(sm, aes(y = Water_depth_cm, x = Site_Code)) +
  geom_boxplot()

ggplot(sm, aes(y = Wetted_width_cm, x = Site_Code)) +
  geom_boxplot()

ggplot(sm, aes(y = Hydraulic_head_cm, x = Site_Code)) +
  geom_boxplot()


##Table 2: calculate mean and std dev 
sm_summary <- sm %>%
  group_by(Site_Code) %>%
  filter(!is.na(Wetted_width_cm)) %>%
  summarise_at(vars(Water_depth_cm, Wetted_width_cm, Hydraulic_head_cm), list(mean = mean, SD = sd))


