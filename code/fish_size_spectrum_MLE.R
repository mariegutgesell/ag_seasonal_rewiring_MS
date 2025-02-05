##size spectrum of fish caught electrofishing 

library(tidyverse)
library(ggplot)
library(readxl)
library(ggbreak)
library(ggpubr)
library(sizeSpectra)


##read in fish data
df <- read.csv("data/fish_baseline_data.csv") %>%
  filter(!is.na(Weight_per_individual_g)) ##removes baseline isotope values and fish from minnow traps (not part of standardized e-fishing sampling)


##MLE estimate of size spectrum slope with bounded power law distribution ------------
##Code for mle_func adapted from Edwards et al., 2017 
mle_func <- function (data = data, season = season) 
{
  dataForYear = dplyr::filter(data, Season == season)
  valCounts = dplyr::select(dataForYear, bodyMass, Number)
  valCounts = ungroup(valCounts)
  xmin = min(valCounts$bodyMass)
  xmax = max(valCounts$bodyMass)
  num.bins = 8
  
  sumNumber = sum(valCounts$Number)
  MLE.valCounts = dplyr::select(valCounts, bodyMass, Number)
  MLE.valCounts = dplyr::summarise(dplyr::group_by(valCounts, 
                                                   bodyMass), Count = sum(Number))
  MLE.valCounts = dplyr::arrange(MLE.valCounts, desc(bodyMass))
  sumCounts = sum(MLE.valCounts$Count)
  if (abs(sumCounts - sumNumber) > 0.001) {
    stop("Check sumCounts in eightMethods.count()")
  }
  MLE.K = dim(MLE.valCounts)[1]
  if (MLE.valCounts[1, "Count"] == 0 | MLE.valCounts[MLE.K, 
                                                     "Count"] == 0) {
    stop("Need first and last counts to be zero in\n          eightMethods.count() for MLE method")
  }
  MLE.xmin = min(MLE.valCounts$bodyMass)
  MLE.xmax = max(MLE.valCounts$bodyMass)
  MLE.sumCntLogMass = sum(MLE.valCounts$Count * log(MLE.valCounts$bodyMass))
  PL.bMLE.counts = 1/(log(MLE.xmin) - MLE.sumCntLogMass/sumCounts) - 
    1
  PLB.minLL = nlm(negLL.PLB.counts, p = PL.bMLE.counts, x = MLE.valCounts$bodyMass, 
                  c = MLE.valCounts$Count, K = MLE.K, xmin = MLE.xmin, 
                  xmax = MLE.xmax, sumclogx = MLE.sumCntLogMass)
  PLB.bMLE = PLB.minLL$estimate
  PLB.minNegLL = PLB.minLL$minimum
  bvec = seq(PLB.bMLE - 5, PLB.bMLE + 5, 1e-05)
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB.counts(bvec[i], x = MLE.valCounts$bodyMass, 
                                     c = MLE.valCounts$Count, K = MLE.K, xmin = MLE.xmin, 
                                     xmax = MLE.xmax, sumclogx = MLE.sumCntLogMass)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == 
      max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  PLB.MLE.stdErr = mean(abs((PLB.MLE.bConf - PLB.bMLE)/1.96))
  mle_results = data.frame(Season = season, Method = as.factor("MLE"), b = PLB.bMLE, confMin = PLB.MLE.bConf[1], 
                                                                        confMax = PLB.MLE.bConf[2], stdErr = PLB.MLE.stdErr, 
                                                                        row.names = NULL) 

  MLE.valCounts = dplyr::mutate(MLE.valCounts, cumSum = cumsum(Count))
  MLE.valCounts = dplyr::mutate(MLE.valCounts, cumProp = cumSum/sumCounts)
  if (MLE.valCounts[dim(MLE.valCounts)[1], "cumProp"] != 1) {
    stop("Check MLE.valcounts in eightMethods.count()")
  }
  MLE.sim = dplyr::tbl_df(data.frame(cumPropSim = seq(as.numeric(MLE.valCounts[1, 
                                                                               "cumProp"]), 1, length = ceiling(sumCounts))))
  MLE.sim = dplyr::mutate(MLE.sim, bodyMassSim = MLE.valCounts[findInterval(cumPropSim, 
                                                                            MLE.valCounts$cumProp), ]$bodyMass)
  return(mle_results)
}

##Seasonal size spectrum slope calculation ----------
hc1 <- df %>%
  filter(SiteCode == "HC") %>%
  rename(bodyMass = "Weight_per_individual_g") %>%
  select(Season, bodyMass, Number_of_Individuals) %>%
  group_by(Season,bodyMass) %>%
  count() %>%
  rename(Number = "n")

hc_summer <- mle_func(data = hc1, season = "Summer") %>%
  mutate(SiteCode = "HC")
hc_spring <- mle_func(data = hc1, season = "Spring") %>%
  mutate(SiteCode = "HC")
hc_fall <- mle_func(data = hc1, season = "Fall") %>%
  mutate(SiteCode = "HC")##come back to this one

ep1 <- df %>%
  filter(SiteCode == "EP4") %>%
  rename(bodyMass = "Weight_per_individual_g") %>%
  select(Season, bodyMass, Number_of_Individuals) %>%
  group_by(Season,bodyMass) %>%
  count() %>%
  rename(Number = "n")

ep_summer <- mle_func(data = ep1, season = "Summer")%>%
  mutate(SiteCode = "EP4")
ep_spring <- mle_func(data = ep1, season = "Spring")%>%
  mutate(SiteCode = "EP4")
ep_fall <- mle_func(data = ep1, season = "Fall")%>%
  mutate(SiteCode = "EP4")

 
at1 <- df %>%
  filter(SiteCode == "AT") %>%
  rename(bodyMass = "Weight_per_individual_g") %>%
  select(Season, bodyMass, Number_of_Individuals) %>%
  group_by(Season,bodyMass) %>%
  count() %>%
  rename(Number = "n")

at1[nrow(at1) + 1,] = list("Spring" ,0.75, 1)

at_summer <- mle_func(data = at1, season = "Summer")%>%
  mutate(SiteCode = "AT")
at_spring <- mle_func(data = at1, season = "Spring")%>%
  mutate(SiteCode = "AT")
at_fall <- mle_func(data = at1, season = "Fall")%>%
  mutate(SiteCode = "AT")


mle_seasonal <- rbind(hc_spring, hc_summer, hc_fall, ep_spring, ep_summer, ep_fall, at_spring, at_summer, at_fall)


##MLE function for all seasons combined -----------
df_mle<- df %>%
  select(SiteCode, Weight_per_individual_g, Number_of_Individuals) %>% 
  rename(Season = "SiteCode", bodyMass = "Weight_per_individual_g") %>% ##just renaming sitecode to be able to use same function 
  select(Season, bodyMass, Number_of_Individuals) %>%
  group_by(Season,bodyMass) %>%
  count() %>%
  rename(Number = "n")

hc_mle <- mle_func(data = df_mle, season = "HC") %>%
  rename(SiteCode = "Season")
ep_mle <- mle_func(data = df_mle, season = "EP4") %>%
  rename(SiteCode = "Season")
at_mle <- mle_func(data = df_mle, season = "AT") %>%
  rename(SiteCode = "Season")

mle_all <- rbind(hc_mle, ep_mle, at_mle)


##Figure 4: Plotting slopes and frequency distribution for all seasons combined --------------
##Plotting Slopes 
mle_all <- mle_all %>%
mutate(site_type = case_when(
  startsWith(SiteCode, "HC") ~ "Low Impact",
  startsWith(SiteCode, "EP4") ~ "Mid Impact",
  startsWith(SiteCode, "AT") ~ "High Impact",
))

mle_all$site_type <- ordered(mle_all$site_type,
                             levels = c("Low Impact", "Mid Impact", "High Impact"))

slopes_all <- ggplot(mle_all, aes(x = site_type, y = b, color = site_type))+
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = confMin, ymax = confMax), width = 0.1, position = position_dodge(width = 0.9))+
  scale_color_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Size Spectrum Slope (b)") +
  xlab("Site Type")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
slopes_all
  

###frequency distribution 
df <- df %>%
  mutate(site_type = case_when(
    startsWith(SiteCode, "HC") ~ "Low Impact",
    startsWith(SiteCode, "EP4") ~ "Mid Impact",
    startsWith(SiteCode, "AT") ~ "High Impact",
  ))

df$site_type <- ordered(df$site_type,
                        levels = c("Low Impact", "Mid Impact", "High Impact"))

fd_all <- df %>%
  ggplot(aes(x = Weight_per_individual_g, group = site_type, fill = site_type)) +
  scale_x_log10()+
  geom_density(alpha = 0.4) +
  scale_fill_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  xlab("Log Body Size (g)")+
  ylab("Frequency Density") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
fd_all

##combine plots for Figure 4
fig_4 <- ggarrange(fd_all, slopes_all, ncol = 2, labels = c("a)", "b)"), font.label = list(colour = "black", size = 12, family = "Times New Roman"))
fig_4

ggsave(fig_4, filename = "figures/Figure_4.png", device = "png")

##Figure S4: Plotting seasonal slopes and frequency distributions ------------
##plotting slopes 
mle_seasonal <- mle_seasonal %>%
  mutate(site_type = case_when(
    startsWith(SiteCode, "HC") ~ "Low Impact",
    startsWith(SiteCode, "EP4") ~ "Mid Impact",
    startsWith(SiteCode, "AT") ~ "High Impact",
  ))

mle_seasonal$site_type <- ordered(mle_seasonal$site_type,
                             levels = c("Low Impact", "Mid Impact", "High Impact"))


slopes_spring <- mle_seasonal %>%
  filter(Season == "Spring") %>%
  ggplot(aes(x = site_type, y = b, color = site_type))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = confMin, ymax = confMax), width = 0.1, position = position_dodge(width = 0.9))+
  scale_color_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Size Spectrum Slope (b)") +
  xlab("Site Type")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
slopes_spring

slopes_summer <- mle_seasonal %>%
  filter(Season == "Summer") %>%
  ggplot(aes(x = site_type, y = b, color = site_type))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = confMin, ymax = confMax), width = 0.1, position = position_dodge(width = 0.9))+
  scale_color_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Size Spectrum Slope (b)") +
  xlab("Site Type")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
slopes_summer

slopes_fall <- mle_seasonal %>%
  filter(Season == "Fall") %>%
  ggplot(aes(x = site_type, y = b, color = site_type))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = confMin, ymax = confMax), width = 0.1, position = position_dodge(width = 0.9))+
  scale_color_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Size Spectrum Slope (b)") +
  xlab("Site Type")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
slopes_fall





##seasonal frequency distributions 
df <- df %>%
  mutate(site_type = case_when(
    startsWith(SiteCode, "HC") ~ "Low Impact",
    startsWith(SiteCode, "EP4") ~ "Mid Impact",
    startsWith(SiteCode, "AT") ~ "High Impact",
  ))

df$site_type <- ordered(df$site_type,
                        levels = c("Low Impact", "Mid Impact", "High Impact"))

fd_spring<-df %>%
  filter(Season == "Spring") %>%
  ggplot(aes(x = Weight_per_individual_g, group = site_type, fill = site_type)) +
  scale_x_log10()+
  geom_density(alpha = 0.4) +
  scale_fill_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  xlab("Body Size (g)")+
  ylab("Frequency Density") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
fd_spring

fd_summer<-df %>%
  filter(Season == "Summer") %>%
  ggplot(aes(x = Weight_per_individual_g, group = site_type, fill = site_type)) +
  scale_x_log10()+
  geom_density(alpha = 0.4) +
  scale_fill_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  xlab("Body Size (g)")+
  ylab("Frequency Density") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_blank(), axis.title.x = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
fd_summer

fd_fall<-df %>%
  filter(Season == "Fall") %>%
  ggplot(aes(x = Weight_per_individual_g, group = site_type, fill = site_type)) +
  scale_x_log10()+
  geom_density(alpha = 0.4) +
  scale_fill_manual(values =c("#999999", "#E69F00", "#56B4E9")) +
  xlab("Log Body Size (g)")+
  ylab("Frequency Density") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
fd_fall



##combine plots for Figure S4
fig_S4 <- ggarrange(fd_spring, slopes_spring, fd_summer, slopes_summer, fd_fall, slopes_fall, ncol = 2, nrow = 3, labels = c("a)", "b)", "c)", "d)", "e)", "f)"), font.label = list(colour = "black", size = 12, family = "Times New Roman"))
fig_S4

ggsave(fig_S4, filename = "figures/Figure_S4.png", device = "png")



