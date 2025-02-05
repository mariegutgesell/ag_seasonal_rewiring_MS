##Code: Seasonal Stream Food Web Stable Isotope Analysis -- Determining Baseline Outliers
##By: Marie Gutgesell
##Date: May 4, 2022

##1)Read in data and select only isotope data 
SI_all <- read.csv("data/fish_baseline_data.csv") %>%
  filter(!is.na(d15N))

##Plot scatter plot of dH and dN of algae and detritus 

baselines <- SI_all %>%
  filter(Common_Name == "Algae" | Common_Name == "Detritus")

baselines_plot <- ggplot(baselines, aes(x = d2H, y = d15N, colour = Common_Name)) +
  geom_point()

baselines_plot

##Determining outliers using Mahalanobis distance -  see follow link for code source and rational (https://towardsdatascience.com/mahalonobis-distance-and-outlier-detection-in-r-cb9c37576d7d)
###Algae
algae <- baselines %>%
  filter(Common_Name == "Algae") %>%
  select(d15N, d2H)

##Find center point
algae_center <- colMeans(algae)

##Find covariance matrix
algae_cov <- cov(algae)

##Ellipse radius from Chi-Square distribution
algae_rad = qchisq(p = 0.95, df = ncol(algae))

##Square root of Chi-Square value
algae_rad <- sqrt(algae_rad)

##Find ellipse coordinates
algae_ellipse <- car::ellipse(center = algae_center, shape = algae_cov, radius = algae_rad, draw = FALSE)

##Ellipse coordinates names should be same w/ algae dataset
algae_ellipse <- as.data.frame(algae_ellipse)
colnames(algae_ellipse) <- colnames(algae)

##Create scatter plot
algae_plot <- ggplot(algae, aes(x = d2H,  y = d15N)) +
  geom_point(size = 2) +
  geom_polygon(data = algae_ellipse, fill = "orange", colour = "orange", alpha = 0.5) +
  geom_point(aes(algae_center[2], algae_center[1]), size = 5, color = "blue") +
  ggtitle("Algae Outlier Test") +
  theme_classic()

algae_plot

##Find Mahalanobis distance
algae_distances <- mahalanobis(x = algae, center = algae_center, cov = algae_cov)
##Cutoff value for distances from Chi-Square dis
##with p = 0.95 and df = 2 which is ncol(algae) - i.e, 2 variables
cutoff <- qchisq(p = 0.95, df = ncol(algae))

##Display observation whose distance is greater than cutoff value (i.e., the outliers)
algae[algae_distances > cutoff,]
##Outlier sample code: AT1-18-AL-01


###Detritus
det <- baselines %>%
  filter(Common_Name == "Detritus") %>%
  select(d15N, d2H)

##Find center point
det_center <- colMeans(det)

##Find covariance matrix
det_cov <- cov(det)

##Ellipse radius from Chi-Square distribution
det_rad = qchisq(p = 0.95, df = ncol(det))

##Square root of Chi-Square value
det_rad <- sqrt(det_rad)

##Find ellipse coordinates
det_ellipse <- car::ellipse(center = det_center, shape = det_cov, radius = det_rad, draw = FALSE)

##Ellipse coordinates names should be same w/ algae dataset
det_ellipse <- as.data.frame(det_ellipse)
colnames(det_ellipse) <- colnames(det)

##Create scatter plot
det_plot <- ggplot(det, aes(x = d2H,  y = d15N)) +
  geom_point(size = 2) +
  geom_polygon(data = det_ellipse, fill = "orange", colour = "orange", alpha = 0.5) +
  geom_point(aes(det_center[2], det_center[1]), size = 5, color = "blue") +
  ggtitle("Detritus Outlier Test") +
  theme_classic()

det_plot

##Find Mahalanobis distance
det_distances <- mahalanobis(x = det, center = det_center, cov = det_cov)
##Cutoff value for distances from Chi-Square dis
##with p = 0.95 and df = 2 which is ncol(algae) - i.e, 2 variables
cutoff <- qchisq(p = 0.95, df = ncol(det))

##Display observation whose distance is greater than cutoff value (i.e., the outliers)
det[det_distances > cutoff,]
##Outlier sample code: AT1-18-DT-01
##Combine plots side by side
library(gridExtra)
grid.arrange(algae_plot, det_plot, ncol = 2)
