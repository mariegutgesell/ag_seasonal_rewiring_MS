##FACET Isotope figure -- joining all 3 isotope results into 1 figure

library(ggpubr)

source("code/isotopic_niche_analysis.R")
source("code/coupling_tp_analysis.R")


ellipse.plot
coupling_boxplot
TP_boxplot

##Final plot 
##Create figure legend for 3 sites 
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Low Impact", "Mid Impact", "High Impact"),
  LegendShape = c("Spring", "Summer", "Fall")
)

data$LegendData <- factor(data$LegendData, levels = c("Low Impact", "Mid Impact", "High Impact"))
data$LegendShape <- factor(data$LegendShape, levels = c("Spring", "Summer", "Fall"))
gplot <- ggplot(data, aes(Xdata, Ydata, color = LegendData, shape = LegendShape)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_shape_manual(values = c( 1, 2, 0)) +
  theme_classic()+
  guides( shape = guide_legend(title = "Season"), colour = guide_legend(title = "Site")) +
  theme(legend.title = element_text(family = "Times New Roman", size = 12, face = "bold"), legend.text = element_text(size = 10, family = "Times New Roman"), legend.position = "bottom")
gplot

leg_fig <- get_legend(gplot)


##remove legends
coupling_boxplot_noleg <- coupling_boxplot +
  theme(legend.position =  "none")
coupling_boxplot_noleg

tp_boxplot_noleg <- TP_boxplot +
  theme(legend.position =  "none")
tp_boxplot_noleg

ellipse_plot_noleg <- ellipse.plot +
  theme(legend.position =  "none")
ellipse_plot_noleg

##make final figure
library("cowplot")

iso_plot_a <- ggdraw() +
  draw_plot(coupling_boxplot_noleg, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(tp_boxplot_noleg, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(ellipse_plot_noleg, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("a)", "b)", "c)"), size = 15, family = "Times New Roman",
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

iso_plot_b <- plot_grid(iso_plot_a, leg_fig, ncol = 1, rel_heights = c(1, .1))
iso_plot_b

ggsave(iso_plot_b, filename = "figures/Figure_2.png", device = "png")

