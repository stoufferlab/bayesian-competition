require(ggpubr)
require(ggplot2)

theme_alba<- theme(#legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = c(0.850, 0.65),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  text = element_text(size = 10),
  panel.background = element_rect(fill = "white",
                                  colour = NA),
  panel.border = element_rect(fill = NA,
                              colour = "grey20"),
  panel.grid = element_line(colour = "white"),
  panel.grid.minor = element_line(size = rel(0.5)),
  strip.background = element_rect(fill = "white",
                                  colour = "grey20"),
  legend.key = element_rect(fill = "white",
                            colour = NA),
  strip.text.x = element_text(size = 12, colour = "black"),
  strip.text.y = element_text(size = 12, colour = "black"))

palette_alba<-ggpubr::get_palette(palette = "YlGn",k=16)[c(7,10,13,16)]