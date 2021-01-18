require("sp")
require("SpatialGraph")
require("tidyverse")
source("code/integration_toolbox.R")

#This function determines the boundary of feasibility for a vector of angle theta
# R_max is the maximum radius calculated for a particular combo of constraints and alpha matrix
feasibility_boundary <- function(theta,
                                 alpha, 
                                 R_max,
                                 rconstraints=NULL,
                                 Nupper=NULL){
  # for every theta value we evaluate if ri and rj (determined by theta and R) along an R sequence
  N <- 500
  R_seq <- seq(0.0001, R_max, length.out = N) %>% as.list()
  R_boundary <- lapply(R_seq,
                       function(R,
                                theta,
                                alpha,
                                rconstraints=NULL,
                                Nupper=NULL){
                         
                         ri <- R * cos(theta)
                         rj <- R * sin(theta)
                         r <- c(ri, rj)
                         
                         #is this particular magnitude feasible
                         feasible <- check_feasibility(r = r,
                                                       alpha = alpha, 
                                                       rconstraints = rconstraints,
                                                       Nupper = Nupper)
                         
                         results <- data.frame( "theta" = theta, 
                                                "R_bound"= R,
                                                "ri"=ri,
                                                "rj"=rj ,
                                                "feasible"= feasible)
                         return(results)
                         
                       } ,theta = theta,
                       alpha = alpha,
                       rconstraints = rconstraints,
                       Nupper = Nupper)
  
                     R_bounds <- do.call(rbind, R_boundary)
                     feasible_bounds <- filter(R_bounds, feasible ==1) 
                     return(feasible_bounds)
                     
}


# this function iteraties over a series of thetas to determine the shape and bounds of the feasibility domain
feasibility_shape<-function( alpha, 
                             R_max,
                             rconstraints=NULL,
                             Nupper=NULL){
  N <- 100
  thetas <- seq(0 , 2*pi, length.out = N) %>% as.list()
  
  #we apply the feasibility_boundary function to a set of thetas
  bound <- lapply(thetas, function(t,
                                    alpha,
                                    R_max, 
                                    rconstraints= rconstraints,
                                    Nupper= Nupper){
     bounded_points  <- feasibility_boundary(theta = t,
                                             alpha = alpha,
                                             R_max =  R_max,
                                             rconstraints =  rconstraints,
                                             Nupper= Nupper)
     return(bounded_points)
    
    
                                            
  }, alpha=alpha,
  R_max =  R_max,
  rconstraints =  rconstraints,
  Nupper= Nupper)
  
  shape <- do.call(rbind, bound)
   
  #we determine the boundary of the points
  coord_points <- chull(x = shape[,"ri"], y=shape[,"rj"]) 
  #we add the first point to complete the polygon
  coord_points<- c(coord_points, coord_points[1])
  shape_bounds <- shape[coord_points,]
  
  return(shape_bounds)
}


#function that returns the shortest distance from the point defined by the growth rates and 
#the bounds of the feasibility domain
#r is a vector of growth rates
#shape is the bounds of the feasibility domain, defined by the funciton feasibility_shape

shortest_distance<-function(r, shape){
  N <- nrow(shape) -1
  distances <- c()
  cc <- c()
  for( i in 1:N){
    
    p0 <- c(shape[i,]$ri, shape[i,]$rj)
    p1 <- c(shape[i+1,]$ri, shape[i+1,]$rj)
    
    line_matrix <- rbind(p0,p1) %>% as.matrix()
    r <- as.matrix(r) %>% t
    dist<- SpatialGraph::pointLineD(xy = line_matrix,
                                    xyp = r )
    dist <- dist$d
    
    data <- data.frame("distance"= dist,
                       #"distance_point1"= dist_pt1,
                       #"distance_point2"= dist_pt2,
                       "p1x"=p0[1],
                       "p1y"=p0[2],
                       "p2x"=p1[1],
                       "p2y"=p1[2])
    
    
    distances <- rbind(distances, data)
    
  }
  
  minimum_distance <- distances[which(distances$dist==min(distances$dist)),]
  # points(r[1], r[2], col="blue", pch=20)
  # lines(x = c(minimum_distance[,"p1x"], minimum_distance[,"p2x"]),
  #       y=c(minimum_distance[,"p1y"], minimum_distance[,"p2y"]), 
  #       lwd=2, 
  #       col="darkgoldenrod")
 return(min(distances$dist))
  #return(minimum_distance)
}



#wrapper of the functions defined above for it's usefull implementation
#its output is the shortest distance from the point of growth rates to the bounds of the feasibility domain
distance_from_limit <- function(alpha, 
                                R_max,
                                rconstraints=NULL,
                                Nupper=NULL,
                                r){
  shape <- feasibility_shape(alpha = alpha,
                             R_max = R_max,
                             rconstraints = rconstraints, 
                             Nupper = Nupper)
  
  distance <- shortest_distance(r = r,
                                shape = shape)
  return(distance)
}


