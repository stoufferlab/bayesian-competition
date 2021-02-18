require("sp")
require("SpatialGraph")
require("tidyverse")
require("polylabelr")


#This function determines the boundary of feasibility domainfor a vector of angle theta and a Radius = R_max
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
  #we check if there is actually an area that is detecte
  nn <- nrow(shape)
  
  if(nn==0){
    shape[1,"ri"] <- 0
    shape[1,"rj"]<- 0
    shape[2,"ri"] <- 0
    shape[2,"rj"]<- 0
    return(shape)
  }else{
    #we determine the boundary of the points
    coord_points <- chull(x = shape[,"ri"], y=shape[,"rj"]) 
    #we add the first point to complete the polygon
    coord_points<- c(coord_points, coord_points[1])
    shape_bounds <- shape[coord_points,]
    return(shape_bounds)
  }
  
  
}


#function to get the euclidean distance between to points
calculate_distance <- function(p1,p2){
  distance <- sqrt( (p2[1]-p1[1])^2 + (p2[2] - p1[2])^2   )
  return(distance)
}

#function that returns the shortest distance from the point defined by the growth rates and 
#the bounds of the feasibility domain
#r is a vector of growth rates
#shape is the bounds of the feasibility domain, defined by the funciton feasibility_shape


shortest_distance<-function(r, 
                            shape,
                            col){
  N <- nrow(shape) -1
  distances <- c()
  cc <- c()
  #for every line defined by two points, we calculate the distance of our growth rates to it 
  for(i in 1:N) {
    # points that define the line
    p0 <- c(shape[i, ]$ri, shape[i, ]$rj)
    p1 <- c(shape[i + 1, ]$ri, shape[i + 1, ]$rj)
    
    #we make everything a matrix
    line_matrix <- rbind(p0, p1) %>% as.matrix()
    r <- as.matrix(r) %>% t
    
    #we get the shortest distance from the point to the line, which is the distance between the point and the perpendicula projections of the line
    dist <- SpatialGraph::pointLineD(xy = line_matrix,
                                     xyp = r)
    
    
    # does the perpendicular projection of the points crosses the segment or not
    cross <- dist$cross
    
    #if it does not, then the shortest distance is the distance to one of the edges of the line
    if (cross == 0) {
      dist_p0 <- calculate_distance(p1 = p0,
                                    p2 = r)
      
      dist_p1 <- calculate_distance(p1 = p1,
                                    p2 = r)
      dist_m <-  ifelse(dist_p0 > dist_p1,
                        dist_p1,
                        dist_p0)
    } else{
      #if it does, then the shortest distance is to the perpendicular projection of the line
      dist_m <- dist$d
    }
    
    
    data <- data.frame(
      "distance" = dist_m,
      #"distance_point1"= dist_pt1,
      #"distance_point2"= dist_pt2,
      "p1x" = p0[1],
      "p1y" = p0[2],
      "p2x" = p1[1],
      "p2y" = p1[2]
    )
    
    
    distances <- rbind(distances, data)
    
  }
  
  #and we return which distance is the shortest
  minimum_distance <-
    distances[which(distances$dist == min(distances$distanc)), ]
  
  # lines(x = c(minimum_distance[,"p1x"], minimum_distance[,"p2x"]),
  #       y=c(minimum_distance[,"p1y"], minimum_distance[,"p2y"]),
  #       lwd=3,
  #       col=col)
  # pch=20)
  
  #But it matters if you are inside or outside the feasibility domain
  return(minimum_distance$distance)
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
  
  #the center of the polygon
  center <- poi(x = shape$ri,
                y = shape$rj)
  
  #we get the shortest distance from the center of the polygon to an edge
  center_distance <- shortest_distance(r= c(center$x, center$y),
                                       shape = shape)
  #we get the shortest distance from our growth rates to an edge
  growth_distance <- shortest_distance(r = r,
                                       shape = shape)
  
  results <- data.frame("center_distance" =center_distance,
                        "growth_distance"= growth_distance)
  
  return(results)
}


