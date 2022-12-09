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
                         
                         # if(feasible){
                         #   points(ri,rj, pch= 20, col = rethinking::col.alpha( "dodgerblue",1 ))}
                         # }else{
                         #   points(ri,rj, pch= 20, col = rethinking::col.alpha( "firebrick", 0.1 ))
                         # }
                         
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
  #we only return the points that are inside the feasibility area
                     R_bounds <- do.call(rbind, R_boundary)
                     feasible_bounds <- filter(R_bounds, feasible ==1) 
                     return(feasible_bounds)
                     
}


# this function iteraties over a series of thetas to determine the shape and bounds of the feasibility domain


feasibility_shape<-function( alpha, 
                             R_max,
                             rconstraints=NULL,
                             Nupper=NULL,
                             N){
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
  
  return(shape)
  }
 
  
  
#function to make sure we detect an area
determine_feasibility_shape <- function(alpha, 
                                        R_max,
                                        rconstraints=NULL,
                                        Nupper=NULL){
  # first we try with wider slices to see if we detect a shape
  shape <- feasibility_shape(alpha = alpha,
                             R_max = R_max,
                             rconstraints = rconstraints,
                             Nupper = Nupper,
                             N = 1000)
  

  nn <- nrow(shape)
  lines_in_shape <- unique(shape$theta) %>% length()
  #if no area is detected or only a thin line, then we run it again with thinner slices
  if(nn==0 | lines_in_shape < 2){
    
    shape <- feasibility_shape(alpha = alpha,
                               R_max = R_max,
                               rconstraints = rconstraints,
                               Nupper = Nupper,
                               N = 2000)
    nn_thin <- nrow(shape)
    
    if(nn_thin == 0){
      #if again no area is detected, we return a data frame with no rows, which can be dealt with later in the function to determine the distance
      return(shape)
      
    }else{
      #we determine the boundary of the points
      coord_points <- chull(x = shape[,"ri"], y=shape[,"rj"]) 
      #we add the first point to complete the polygon
      coord_points<- c(coord_points, coord_points[1])
      shape_bounds <- shape[coord_points,]
      return(shape_bounds)
    }
    
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
#shape is the bounds of the feasibility domain, defined by the funciton determine_feasibility_shape


shortest_distance<-function(r, 
                            shape,
                            feasibility){

  N <- nrow(shape) -1
  data_dist <- c()
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
      type <- "point"
      
    } else{
      #if it does, then the shortest distance is to the perpendicular projection of the line
      dist_m <- dist$d
      type <- "line"
    }
    
    
    data <- data.frame(
      "distance_from_edge" = dist_m,
      #"distance_point1"= dist_pt1,
      #"distance_point2"= dist_pt2,
      "p1x" = p0[1],
      "p1y" = p0[2],
      "p2x" = p1[1],
      "p2y" = p1[2],
      "type"= type
    )
    
    
    data_dist <- rbind(data_dist, data)
    
  }
  
  #and we return which distance is the shortest
   minimum_distance <- data_dist[which(data_dist$distance_from_edge == min(data_dist$distance_from_edge)), ]
   #but it has to be ONLY ONe
   minimum_distance <- unique(minimum_distance$distance_from_edge)
   
  # 
    # lines(x = c(minimum_distance[,"p1x"], minimum_distance[,"p2x"]),
    #       y=c(minimum_distance[,"p1y"], minimum_distance[,"p2y"]),
    #       lwd=3,
    #       col=2)

  #But it matters if you are inside or outside the feasibility domain, with a sign
   if (feasibility){
     return(minimum_distance)
   }else{
     return( - minimum_distance)
   }
}


#wrapper of the functions defined above for it's usefull implementation
#its output is the shortest distance from the point of growth rates to the bounds of the feasibility domain
distance_from_limit <- function(alpha, 
                                R_max,
                                rconstraints=NULL,
                                Nupper=NULL,
                                r,
                                feasibility){
  
  feas <- ifelse(feasibility,1,-1)
  shape <- determine_feasibility_shape(alpha = alpha,
                             R_max = R_max,
                             rconstraints = rconstraints, 
                             Nupper = Nupper)
  nn <- nrow(shape)
  lines_in_shape <- unique(shape$theta) %>% length()
  #if there is no detectable shape
  if(nn==0 | lines_in_shape < 2){
    print(1)
    distance <- calculate_distance( p1 =  r,
                                    p2 = c(0,0))
    distance <- distance*feas
    results <- data.frame("center_distance" = distance,
                          "growth_distance"= 0)
    
    #return(results)
  }else{
    #if the feasibility domain is only a line, and not a volume, then we can not detect the minimum distance to a boundary
      col1 <- rethinking::col.alpha("black", alpha = 0.3)
      lines(shape$ri, shape$rj, col=col1)
      #the center of the polygon
      center <- poi(x = shape$ri,
                    y = shape$rj)
      
      #we get the shortest distance from the center of the polygon to an edge, it is always feasible
      center_distance <- shortest_distance(r= c(center$x, center$y),
                                           shape = shape,
                                           feasibility = 1)
      #we get the shortest distance from our growth rates to an edge
      growth_distance <- shortest_distance(r = r,
                                           shape = shape,
                                           feasibility = feasibility)
      
      results <- data.frame("center_distance" =center_distance,
                            "growth_distance"= growth_distance)
      
     # return(results)
      
    }
  
}


