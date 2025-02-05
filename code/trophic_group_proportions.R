##Stomach Content Data for Seasonal Stream Food Webs -- Creating summary figure

##Load Packages -- 
library(tidyverse)
library(ggplot2)
library(ggpubr)

##1) Read in fish data and filter to creek chub w/ stomach content data
fish_meta_cc <- read.csv("data/fish_baseline_data.csv") %>%
  filter(Species_Code == "CC") %>%
  filter(!is.na(Stomach_Contents)) %>%
  select(Season, SiteCode, Species_Code, Sample_ID, Stomach_Contents)



##2) Add Trophic Group depending on stomach contents 
##Stomach content groups:
##Empty = empty
##Not taken = unknown
##Undetermined = unknown
##Algae = Algae
##Detritus = Detritus
##Insects = Insect (also: beetles, worms,)
##Fish = Fish
##Algae + insect = Omnivore (also: plant matter + insect)
##not including parasites as part of stomach contents
##brown goo? = Unknown

stomach_content_types <- as.data.frame(unique(fish_meta_cc$Stomach_Contents))

fish_meta_cc <- fish_meta_cc %>%
  mutate(stomach_content_group = case_when(
    startsWith(Stomach_Contents, "algae, insect") ~ "Omnivore",
    startsWith(Stomach_Contents, "insects, algae") ~ "Omnivore",
    startsWith(Stomach_Contents, "algae/bug") ~ "Omnivore",
    startsWith(Stomach_Contents, "plant matter/insect") ~ "Omnivore",
    startsWith(Stomach_Contents, "gastropod/algae") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus/insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus, insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus, fish, insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus/bug") ~ "Omnivore",
    startsWith(Stomach_Contents, "worms") ~ "Predator",
    grepl("empty", Stomach_Contents) ~ "Empty",
    grepl("not taken", Stomach_Contents) ~ "Unknown",
    grepl("unid", Stomach_Contents) ~ "Unknown",
    grepl("Unid", Stomach_Contents) ~ "Unknown",
    grepl("worm/algae", Stomach_Contents) ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus") ~ "Primary Consumer (detritus)",
    startsWith(Stomach_Contents, "insect") ~ "Predator",
    startsWith(Stomach_Contents, "insects, algae") ~ "Omnivore",
    grepl("beetle", Stomach_Contents) ~ "Predator", 
    grepl("insects/worm", Stomach_Contents) ~ "Predator",
    startsWith(Stomach_Contents, "algae") ~ "Primary Consumer (algae)",
    grepl("UNK", Stomach_Contents) ~"Unknown",
    grepl("unknown", Stomach_Contents) ~"Unknown",
    grepl("damselfly", Stomach_Contents) ~"Predator",
    grepl("chironomid", Stomach_Contents) ~"Predator",
    grepl("brown goo", Stomach_Contents) ~ "Unknown"
  ))


##3) Calculate # of fish per group and total sample size (Table S4)
fish_per_cat <- fish_meta_cc %>%
  group_by(SiteCode, Season, stomach_content_group) %>%
  summarise(count = n())

sample_size <- fish_meta_cc %>%
  group_by(SiteCode, Season) %>%
  count()

##4) Calculate proportion of fish in each trophic groups
stomach_content_prop_function <- function(x){
  df_prep <- x
  n <- nrow(df_prep)
  df_prep_count <- df_prep %>%
    group_by(stomach_content_group) %>%
    summarise(count = n())
  prop <- df_prep_count$count/n
  prop <- as.data.frame(cbind(df_prep_count, prop)) 
  prop <- merge(df_prep, prop)
}

stomach_prop <- split(fish_meta_cc, paste0(fish_meta_cc$Season, fish_meta_cc$SiteCode)) %>%
  map(stomach_content_prop_function) %>%
  bind_rows() %>%
  select(Season, SiteCode, stomach_content_group, count, prop) %>%
  unique() %>%
  mutate(Site_Type = case_when(
    startsWith(SiteCode, "H") ~ "Low Impact",
    startsWith(SiteCode, "A") ~ "High Impact", 
    startsWith(SiteCode, "E") ~ "Mid Impact"
  ) )

##5) Calculate mean, SD and SE of omnivores and predators 
##calculate mean, sd, se of omnivores
stomach_prop_omni <- stomach_prop %>%
  ungroup() %>%
  filter(stomach_content_group == "Omnivore") %>%
  select(Site_Type, Season, stomach_content_group, prop) %>%
  group_by(Site_Type, stomach_content_group) %>%
  summarise_at(vars(prop), list(prop_mean = mean, prop_sd = sd)) %>%
  mutate(prop_se = prop_sd/3) ##sample size is 3 because have 3 seasons

##calculate mean, sd, se of predators
stomach_prop_insect <- stomach_prop %>%
  filter(stomach_content_group == "Predator") %>%
  select(Site_Type, Season, stomach_content_group, prop) %>%
  group_by(Site_Type, stomach_content_group) %>%
  summarise_at(vars(prop), list(prop_mean = mean, prop_sd = sd)) %>%
  mutate(prop_se = prop_sd/3) ##sample size is 3 because have 3 seasons
##calculate magnitude of difference between low and high impact proportion of predators
mag_diff_pred <- 0.6500000/0.2623188



##6) Create figure 3  ------------
##Plot seasonal trophic group proportions (Figure 3A)
##Prep dataframe to format properly
stomach_prop$Season <- ordered(stomach_prop$Season, 
                               levels = c("Spring", "Summer", "Fall"))

stomach_prop$Site_Type <- ordered(stomach_prop$Site_Type, 
                                  levels = c("Low Impact", "Mid Impact", "High Impact"))

stomach_prop$stomach_content_group <- ordered(stomach_prop$stomach_content_group, 
                                              levels = c("Primary Consumer (algae)", "Primary Consumer (detritus)", "Empty", "Predator", "Omnivore", "Unknown"))

stomach_content_prop_plot <- ggplot(stomach_prop, aes(x = Season, y = prop, fill = stomach_content_group)) +
  geom_col(position = "stack", colour = "black")  +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2")) +
  facet_wrap(~Site_Type)+
  ylab("Proportion of Individuals") +
  labs(fill = "Trophic Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 
stomach_content_prop_plot

##Plot mean % omnivores and predators 
omni_prop_plot <- ggplot(stomach_prop_omni, aes(x = Site_Type, y = prop_mean, fill = stomach_content_group )) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = prop_mean - prop_se, ymax = prop_mean + prop_se), width = 0.3, position = position_dodge(width = 0.9)) +
  theme_classic() +
  scale_fill_manual(values = c("#F0E442")) +
  ylab("Mean Proportion of Individuals") +
  labs(fill = "Trophic Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 
omni_prop_plot

insect_prop_plot <- ggplot(stomach_prop_insect, aes(x = Site_Type, y = prop_mean, fill = stomach_content_group)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = prop_mean - prop_se, ymax = prop_mean + prop_se), width = 0.3, position = position_dodge(width = 0.9)) +
  theme_classic() +
  scale_fill_manual(values = c("#009E73")) +
  ylab("Mean Proportion of Individuals") +
  labs(fill = "Trophic Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 
insect_prop_plot

##Create final figure 
stomach_content_fig <- ggarrange(stomach_content_prop_plot, 
                                 ggarrange(omni_prop_plot, insect_prop_plot, nrow = 2, labels = c("b)", "c)"), font.label = list(colour = "black", size = 12, family = "Times New Roman")),
                           ncol = 2,  labels = "a)",font.label = list(colour = "black", size = 12, family = "Times New Roman"))
stomach_content_fig

ggsave(stomach_content_fig, filename = "figures/Figure_3.png", device = "png")

##7)ANOVA and tukey ------------
##omnivores: 
omni <- stomach_prop %>%
  filter(stomach_content_group == "Omnivore") %>%
  ungroup() %>%
  select(SiteCode, Season, stomach_content_group, prop) %>%
  complete(SiteCode, Season, fill = list(prop=0))
o_aov <- aov(prop ~ SiteCode, omni)
summary(o_aov)
TukeyHSD(o_aov)

#predator:
pred <- stomach_prop %>%
  filter(stomach_content_group == "Predator") %>%
  ungroup() %>%
  select(SiteCode, Season, stomach_content_group, prop) %>%
  complete(SiteCode, Season,  fill = list(prop=0))
p_aov <- aov(prop ~ SiteCode, pred)
summary(p_aov)
TukeyHSD(p_aov)


