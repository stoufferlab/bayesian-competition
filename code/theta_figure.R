
#Function to make the cloud figure. It is essentially a wrapper of the feasibility_wrapper, plus a few options regarding plots


make_figure_theta <-function(vero_model,
           vero_function,
           vero_exp,
           vero_name,
           trcy_model,
           trcy_function,
           trcy_exp,
           trcy_name,
           env,
           show_top = TRUE,
           show_left = TRUE) {
    col1 <- rethinking::col.alpha("mediumseagreen", .3)
    col2 <- rethinking::col.alpha("grey50", .3)
    
    post        <-
      feasibility_wrapper(
        vero_model = vero_model,
        vero_function = vero_function,
        vero_exp = vero_exp,
        vero_name = vero_name,
        trcy_model = trcy_model,
        trcy_function = trcy_function,
        trcy_exp = trcy_exp,
        trcy_name = trcy_name,
        env = env
      )
    coexistence <- sum(post$feasibility_results) / nrow(post)
    coexistence <- round(coexistence, digits = 2)
    
    
    
    
    theta_omega <- ggplot(post) +
      geom_point(
        mapping = aes(
          x = omega_results,
          y = theta_results,
          col = as.factor(feasibility_results)
        ),
        show.legend = FALSE
      ) +
      theme_alba +
      scale_color_manual(values = c(col2, col1)) +
      xlim(0, 0.5)  +
      ylim(0, 45)
    
    
    if (show_top & show_left) {
      pp <-
        annotate_figure(
          theta_omega,
          top = text_grob(vero_name, size = 16, face = "bold"),
          left = text_grob(
            trcy_name,
            rot = 90,
            size = 16,
            face = "bold"
          )
        )
    }
    
    if (show_top & show_left == FALSE) {
      pp <-
        annotate_figure(theta_omega, top = text_grob(vero_name, size = 16, face =
                                                       "bold"))
    }
    
    
    if (show_top == FALSE & show_left) {
      pp <-  annotate_figure(theta_omega,
                             left = text_grob(
                               trcy_name,
                               rot = 90,
                               size = 16,
                               face = "bold"
                             ))
      
    }
    return(pp)
    
}
