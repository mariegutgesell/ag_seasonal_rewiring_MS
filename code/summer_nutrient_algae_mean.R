###Calculate mean and SD of summer nutrient concentrations and algae cover
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)

##Calculate mean and SD of summer nutrient data
n_p <- read.csv("data/summer_nutrient_data.csv")
colnames(n_p)[2] <- "Date"

n_p$Date <- as.Date(n_p$Date, format = "%d-%b-%y")

n_p <- n_p %>%
  mutate(Site_Type = case_when(
    startsWith(Site_Code, "H") ~ "Low Impact",
    startsWith(Site_Code, "A") ~ "High Impact",
    startsWith(Site_Code, "E") ~ "Mid Impact",
  ))

##remove may -- so only summer values
n_p_summer <- n_p %>%
  filter(Date != "2018-05-10")

n_p_avg <- n_p_summer %>%
  group_by(Site_Type) %>%
  summarise_at(vars(TN_mg_per_L, TP_mg_per_L), list(mean = mean, sd = sd))


##Looking at seasonal % algae cover
df_a <- read.csv("data/summer_algae_cover.csv")


##calculate mean summer algal cover
df_a_summer <- df_a %>%
  filter(Month %in% c("June", "July", "August")) %>%
  filter(!is.na(`Algae.Coverage.Category.Midpoint`)) %>%
  group_by(Site_Code) %>%
  summarise_at(vars(`Algae.Coverage.Category.Midpoint`), list(mean = mean))
 

