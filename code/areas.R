library(ggplot2)
library(ggpubr)

# figures for the paper ---------------------------------------------------

source("code/access_results.R")

sunny_results <- filter(final_results, environment=="open") %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )
woody_results <- filter(final_results, environment=="woody")  %>% filter( vero_model != "Lotka-Volterra")%>% filter( trcy_model != "Lotka-Volterra" )



col1 <- rethinking::col.alpha("mediumseagreen",1)
col2 <- rethinking::col.alpha("#0e77be",1)
col3<- rethinking::col.alpha("#f0bb3f",1)



col4 <- rethinking::col.alpha("mediumseagreen",.7)
col5 <- rethinking::col.alpha("#0e77be",.7)
col6<- rethinking::col.alpha("#f0bb3f",.7)



xx <- expression("Area in monoculture,"~ gamma)
yy <- expression("Biologically-constrained feasibility domain,"~ beta)


new_order <- c("both", "trcy", "vero")


sunny_results$outcome_mean <- factor(sunny_results$outcome_mean, levels = new_order)
woody_results$outcome_mean <- factor(woody_results$outcome_mean, levels = new_order)


sun1<-ggplot(sunny_results) +
  geom_point(
    mapping = aes(
      x = area,
      y = area_alone,
      col = as.factor(outcome),
    ),
    shape=19,
    show.legend = TRUE,
    size = 1
  )+
  scale_color_manual(
    values = c(col4, col5, col6),
    name = "Posterior prediction",
    labels = c("Coexistence",
               expression(italic("T. cyanopetala")),
               expression(italic(
                 "V. rosea"
               ))
    ))+
  facet_grid(trcy_model ~ vero_model)  +
  geom_abline(
    intercept = 0,
    slope = 0,
    linetype = "dashed",
    col = "grey50"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 12),
    legend.text.align = 0,
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin()
  ) +
  xlab(xx)+
  ylab(yy)+
  geom_abline()

open<- sun1 +
  geom_point(
    mapping = aes(
      x = proportion_mean,
      y = distance_mean,
      fill = as.factor(outcome_mean)
    ),
    col="black",
    shape = 24,
    stroke=0.5,
    size = 3
  ) +
  scale_fill_manual(
    values = c(col1,col2,col3),
    drop =FALSE,
    name = "Median prediction",
    labels = c("Coexistence",
               expression(italic("T. cyanopetala")),
               expression(italic(
                 "V. rosea"
               )))
  ) 

open <- annotate_figure(open,
                        top = text_grob( "Vellia rosea", size = 12, face = "italic"),
                        right = text_grob( "Trachymene cyanopetala", size = 12, face = "italic", rot=-90)
)


ggsave(open, filename = "../bayesian_competition_ms/sunny_results.pdf", width = 8, height = 6)
