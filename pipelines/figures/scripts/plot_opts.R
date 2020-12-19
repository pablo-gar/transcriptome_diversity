# PLOT
p.label_size <- 3.2
p.text_small <- 8
p.text_mid <- 10

p.point_size_baby <- 0.2
p.point_size_tiny <- 0.5
p.point_size_small <- 1
p.point_size_medium <- 2


# COLORS
p.positive_color <- '#e0cd51'
p.negative_color <- '#5089ec'
p.positive_color_green <- '#517e24'
p.negative_color_red <- '#791c9a'
p.bar_gray_fill <- 'grey70'
p.bar_gray_colour <- 'black'
p.abline_colour <- 'cadetblue4'


p.color_dark_blue <- '#4068ad'

# Theme
theme_pars <- list(text=element_text(family='sans'),
                   axis.text=element_text(size=p.text_small),
                   axis.title=element_text(face='bold', size=p.text_mid),
                   strip.text=element_text(face='bold.italic', size=p.text_small),
                   strip.background=element_blank(),
                   legend.text=element_text(size=p.text_small),
                   legend.title=element_text(size=p.text_mid)
                   )

guide_pars <- list(color = guide_legend(override.aes = list(size = 0.2)))
