require(ggpubr)
require(ggplot2)

theme_alba<- theme(#legend.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.position = c(0.70, 0.70),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  text = element_text(size = 10),
  panel.background = element_rect(fill = "white",
                                  colour = NA),
  panel.border = element_rect(fill = NA,
                              colour = "grey20"),
  panel.grid = element_line(colour = "grey92"),
  panel.grid.minor = element_line(size = rel(0.5)),
  strip.background = element_rect(fill = "grey85",
                                  colour = "grey20"),
  legend.key = element_rect(fill = "white",
                            colour = NA))

palette_alba<-ggpubr::get_palette(palette = "YlGn",k=16)[c(7,10,13,16)]