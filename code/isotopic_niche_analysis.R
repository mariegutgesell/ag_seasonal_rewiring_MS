##Seasonal Stream Food Web Isotope Analysis - Isotopic Niche Width using nicheROVER on trophic position and % coupling metrics

##Using nicheRover

source("code/coupling_tp_analysis.R")
rm(list = ls()[!ls() %in% c("seasonal_df_split_pred_tp_bound_noout")])

library(nicheROVER)

###1) Format data for use in nicheRover (TP and %TC as columns, each row is an observation)
##Create column that combines site and season
seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Group = case_when(
    startsWith(Sample_ID, "HC1") ~ "HC_spring",
    startsWith(Sample_ID, "HC2") ~ "HC_summer", 
    startsWith(Sample_ID, "HC3") ~ "HC_fall",
    startsWith(Sample_ID, "EP41") ~ "EP4_spring",
    startsWith(Sample_ID, "EP42") ~ "EP4_summer", 
    startsWith(Sample_ID, "EP43") ~ "EP4_fall",
    startsWith(Sample_ID, "AT1") ~ "AT_spring",
    startsWith(Sample_ID, "AT2") ~ "AT_summer", 
    startsWith(Sample_ID, "AT3") ~ "AT_fall"
  ) )

data <- seasonal_df_split_pred_tp_bound_noout %>%
  select(Group, TC_pred_logit, TP_pred) 

data$Group <- ordered(data$Group,
                                    levels = c("HC_spring", "HC_summer", "HC_fall", "EP4_spring", "EP4_summer", "EP4_fall", "AT_spring", "AT_summer", "AT_fall"))

aggregate(data[2:3], data[1], mean) ##calculates mean TP and %TC for each site:season

##Generate the posterior distributions of u (mean) and V (variance) for isotope values of each species with the default prior
# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e4
system.time({
  fish.par <- tapply(1:nrow(data), data$Group,
                     function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))
})

##Plot mean and variance parameters
clrs <- c("black", "red", "blue", "orange", "green", "purple", "yellow", "grey", "pink") # colors for each species

# mu1 (%TC), mu2 (TP), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

##plot all mu and sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)



##Generate 10000 draws from posterior distributions 
nsamples <- 10000
fish.par <- tapply(1:nrow(data), data$Group,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))


##Calculate probability of overlap for 95% and 40% (SEA) niche areas 
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, 0.40))


#The mean overlap metrics calculated across iteratations for both niche sizes (Table )
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100 
mean_overlap <- round(over.mean, 2)
mean_overlap<- as.data.frame(mean_overlap)
#write.csv(mean_overlap, "outputs/mean_niche_overlap.csv")


##determine credible interval of overlap (Table S7 and S8)
over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
cred_overlap_95 <- round(over.cred[,,,1]) # display alpha = .95 niche region
cred_overlap_95 <- as.data.frame(cred_overlap_95)

cred_overlap_40 <- round(over.cred[,,,2])
cred_overlap_40 <- as.data.frame(cred_overlap_40)


# Generate estimates of SEA niche area (note: results not presented in MS)
p.ell <- pchisq(1,2) ##SEA 
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = p.ell)
})

rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

dev.off()

boxplot(fish.size, col = clrs,
        ylab = "SEA Niche Area", xlab = "Group")


##Create plot ellipse plot of SEA (Figure 2c) -----------------------
seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Site_Type = case_when(
    startsWith(Group, "HC") ~ "Low Impact",
    startsWith(Group, "EP4") ~ "Mid Impact", 
    startsWith(Group, "AT") ~ "High Impact"
  ))


data_2 <- seasonal_df_split_pred_tp_bound_noout %>%
  select(Site_Type, Season, TC_pred_logit, TP_pred)

data_2$Site_Type <- ordered(data_2$Site_Type,
                          levels = c("Low Impact", "Mid Impact", "High Impact"))

data_2$Season <- ordered(data_2$Season, 
                         levels = c("Spring", "Summer", "Fall"))

first.plot <- ggplot(data = data_2, 
                     aes(x = TC_pred_logit, 
                         y = TP_pred)) + 
  geom_point(aes(color = Site_Type, shape = Season), size = 5) +
  ylab("Trophic Position") +
  xlab("% Terrestrial Energy Use (logit)") + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) + 
  theme_classic()
first.plot


#p.ell <- pchisq(1,2) ##creates predictive ellipse that contains approx. p % of data
p.ell <- pchisq(1,2)
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(Season, Site_Type), 
                   fill = Site_Type,
                   color = Site_Type), 
               alpha = 0.25, 
               level = 0.40,
               type = "norm",
               geom = "polygon") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman")) 

ellipse.plot


