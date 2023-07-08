# code for plotting gpd density in terms of the quotient x/sigma, where sigma
# is the scale of the GPD distribution

# import necessary libraries
library(latex2exp)
library(ggplot2)

# vector with points for calculating density
x <- seq(1, 6.5 , length.out = 1e4)

# specify 3 shape parameters to showcase the different possibilities of it being
# positive, negative or zero
shape_0 <- quaketools::dgpd(x, shape = 0)
shape_positive <- quaketools::dgpd(x, shape = 0.5)
shape_negative <- quaketools::dgpd(x, shape = -0.2)

# dataframe with data for plotting
df <- data.frame( x = rep(x, 3), y = c(shape_0, shape_positive, shape_negative),
                  z = rep(c('shape_0', 'shape_positive', 'shape_negative'), each = length(x)))
df$z <- factor(df$z, levels = c('shape_negative', 'shape_0', 'shape_positive'))

# plot
gpd_density <- ggplot(df, aes(x = x, y = y, color = z))+
  geom_line(linewidth = 1)+
  theme_classic()+
  xlab(unname(TeX("$x/\\sigma_u $")))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(family = 'sans', size=13),
        axis.text.x = element_text(family = 'sans', size=11),
        axis.text.y = element_text(family = 'sans', size=11),
        title = element_text(family = 'sans', size = 8),
        legend.title = element_blank(),
        legend.text = element_text(family = 'sans', size = 11),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        legend.position = c(0.8, 0.8),
        legend.text.align = 0)+
  scale_color_discrete(labels = unname(TeX(c("$  \\xi = -0.2$","$\\xi = 0$", "$\\xi = 0.5$"))))+
  expand_limits(y = 0) +
  coord_cartesian(expand = FALSE, clip = "off")

ggsave(here::here('01-week', 'outputs', 'figures', 'gpd_density.png'),
       plot = gpd_density,
       units = "px",
       width = 997,
       height = 683,
       dpi = 300)
