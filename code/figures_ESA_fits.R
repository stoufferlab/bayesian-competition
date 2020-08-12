### figure  for ESA


BEV_vero <- readRDS("~/bayesian-competition/BEV_vero.RDS")
RC_vero <- readRDS("~/bayesian-competition/RC_vero.RDS")
LAW_vero <- readRDS("~/bayesian-competition/LAW_vero.RDS")


make_conditions()

df <- data.frame("conspecifics"= seq(0,80,1))


bev<-marginal_effects(LAW_vero ,effects = "conspecifics", select_points = 10) 

plot(bev, plot = FALSE)[[1]] + 
  theme_bw()+
  geom_line(color="black", size=1)+
  xlim(0,80)  +
  ylim(0,20)  +
  xlab("Number of conspecifics")+
  ylab("Number of seeds produced")  +
  geom_point(
    aes(x = conspecifics, y = totalseeds), 
    # this is the key!
    data =vero_focal , 
    color = "olivedrab",
    size = 2,
    alpha = .8,
    # This tells it to ignore the ymin and ymax settings used elsewhere
    inherit.aes = FALSE
  )

