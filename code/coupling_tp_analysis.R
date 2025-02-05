##Code: Seasonal Stream Food Web Stable Isotope Bi-plots - exploring data across site
##By: Marie Gutgesell
##Date: May 4, 2022

library(car)
library(tidyverse)
##1)Read in data and select only isotope data 
SI_cc_al_dt <- read.csv("data/fish_baseline_data.csv") %>%
  filter(!is.na(d15N))




##2) Conduct d2H correction for fish  ---------------
##Equations For Creek Chub: 
#### w2 = 1 - (1-w1)^2 
####dHadjusted = [dHraw - (w2*dHwater)]/1-w2

##Equations from Page et al., 2017
##dH20 estimation taken as mean of dH20 values of two sites that span our study sites from Gibson et al., 2020 - dHwater = -62.7
##w1 estimates from Solomon et al., 2009 - use w = 0.17 (avg. for all organisms) 
##Calculate w2
w2 = 1 - (1-0.17)^2

creekchub <- SI_cc_al_dt %>%
  filter(Common_Name == "Creek Chub") %>%
  select(Common_Name, Sample_ID, d2H)

output_dHadjust_cc <- matrix(nrow = 86, ncol = 1)
dH_adjust_cc <- for(i in 1:nrow(creekchub)){
  output_dHadjust_cc[i,] <- (creekchub$d2H[i] - (0.31 * (-62.7)))/(1-0.31)
}
output_dHadjust_cc
creekchub_adjusted <- cbind(creekchub, output_dHadjust_cc)

##Join adjusted values back to main df
SI_cc_al_dt_adj <- left_join(SI_cc_al_dt, creekchub_adjusted, by = c("Common_Name", "Sample_ID", "d2H")) 

##3) Remove baseline outliers (see baseline_outliers.R code to see how outliers were determined)
SI_cc_al_dt_adj <- SI_cc_al_dt_adj %>%
  filter(Sample_ID != "AT1-18-AL-01") %>%
  filter(Sample_ID != "AT1-18-DT-01")

##4) Calculate trophic position and % terrestrial energy use 
##Calculating TP & dN using bound %TC b/w 0-1 --------
si_mixing_model_pred_bound <- function(x){
  df_prep <- x
  AB <- df_prep %>%
    filter(Common_Name == "Algae")
  AB_mean_dH <- mean(AB$d2H)
  AB_mean_dN <- mean(AB$d15N)
  TB <- df_prep %>%
    filter(Common_Name == "Detritus")
  TB_mean_dH <- mean(TB$d2H)
  TB_mean_dN <- mean(TB$d15N)
  df_pred <- df_prep %>%
    filter(Common_Name =="Creek Chub")
  for(i in 1:nrow(df_pred)){
    TC_pred = ((df_pred$output_dHadjust_cc - AB_mean_dH)/(TB_mean_dH - AB_mean_dH))
    TC_pred[TC_pred>=1] <- 0.999
    TC_pred[TC_pred<=0] <- 0.001
    TP_pred = 1 + (df_pred$d15N - ((AB_mean_dN*(1 - TC_pred)) + (TB_mean_dN*TC_pred)))/3.4
    dN_pred = df_pred$d15N - (((1-TC_pred)*AB_mean_dN)+(TC_pred*TB_mean_dN))
    TC_TP_df_pred <- cbind(TC_pred, TP_pred, dN_pred)
  }
  TC_TP_df_pred <- as.data.frame(TC_TP_df_pred)
  TC_TP_df_pred$Common_Name <- df_pred$Common_Name
  TC_TP_df_pred$SiteCode <- df_pred$SiteCode
  TC_TP_df_pred$Season <- df_pred$Season
  TC_TP_df_pred$Sample_ID <- df_pred$Sample_ID
  TC_TP_df_pred
}

seasonal_df_split_pred_tp_bound <- split(SI_cc_al_dt_adj, paste0(SI_cc_al_dt_adj$SiteCode, SI_cc_al_dt_adj$Season)) %>%
  map(si_mixing_model_pred_bound) %>%
  bind_rows()


##5) Logit Transformation and Remove Outliers 
tp_dn <- left_join(seasonal_df_split_pred_tp_bound, SI_cc_al_dt_adj, by = c("Season", "SiteCode", "Common_Name", "Sample_ID")) %>%
  mutate(TC_pred_logit = logit(TC_pred))
##Function to remove outliers
##removing outliers
#' Replace outliers
#'
#' Replace outliers with NA. Outliers are defined as values that fall outside plus or minus
#' 1.5 * IQR.
#'
#' @return Numeric vector of same length.
#' @param x Numeric vector to replace.
#' @param na.rm Should NA values be replaced beforehand?
#'
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}


seasonal_df_split_pred_tp_bound_noout <- tp_dn %>% 
  group_by(SiteCode, Season) %>% 
  mutate(TC_pred_logit = remove_outliers(TC_pred_logit)) %>% 
  mutate(TP_pred = remove_outliers(TP_pred)) %>%
  ungroup() %>% 
  filter(!is.na(TC_pred_logit)) %>%
  filter(!is.na(TP_pred))



##6) Plot boxplots for % TEU and TP ---------------
seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Site_Type = case_when(
    startsWith(Sample_ID, "HC") ~ "Low Impact",
    startsWith(Sample_ID, "EP4") ~ "Mid Impact",
    startsWith(Sample_ID, "AT") ~ "High Impact"
  ) )

seasonal_df_split_pred_tp_bound_noout$Site_Type <- ordered(seasonal_df_split_pred_tp_bound_noout$Site_Type, 
                                                           levels = c("Low Impact", "Mid Impact", "High Impact"))

seasonal_df_split_pred_tp_bound_noout$Season <- ordered(seasonal_df_split_pred_tp_bound_noout$Season,
                                                        levels = c("Spring", "Summer", "Fall"))

coupling_boxplot <- ggplot(seasonal_df_split_pred_tp_bound_noout, aes(x = Season, y = TC_pred_logit, fill = Site_Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  ylab("% Terrestrial Energy Use (logit)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman"))+
  guides(fill = guide_legend(title="Site"))
coupling_boxplot


TP_boxplot <- ggplot(seasonal_df_split_pred_tp_bound_noout, aes(x = Season, y = TP_pred, fill = Site_Type)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Trophic Position") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman"))+ 
  guides(fill = guide_legend(title="Site"))

TP_boxplot



##ANOVA of % coupling and trophic position ----------------
coupling_anova <- aov(TC_pred_logit ~ Site_Type + Season, seasonal_df_split_pred_tp_bound_noout)
summary(coupling_anova)

coupling_tukey <- TukeyHSD(coupling_anova)
coupling_tukey_table <- coupling_tukey[["Site_Type:Season"]]

tp_anova <- aov(TP_pred ~ Site_Type + Season, seasonal_df_split_pred_tp_bound_noout)
summary(tp_anova)
tp_tukey <- TukeyHSD(tp_anova)
tp_tukey_table <- tp_tukey[["Site_Type:Season"]]




