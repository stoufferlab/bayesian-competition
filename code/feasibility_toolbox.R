require("sp")
require("SpatialGraph")
require("tidyverse")
require("polylabelr")

#FUNCTIONS TO GET THE SHAPE AND AREA OF THE FEASIBILITY DOMAIN, AS WELL AS THE SHORTEST DISTANCE FROM AN EDGE OF THE GROWTH RATES


#function to ge the invese of a matrix, just to save time
inverse_matrix<-function(alpha){
  
  deteminant_alpha <-
    (alpha[1, 1] * alpha[2, 2]) - (alpha[2, 1] * alpha[1, 2])
  
  inverse_det  <- 1 / deteminant_alpha
  
  
  adjugate <- matrix(NA, nrow = 2, ncol = 2)
  adjugate[1, 1] <- alpha[2, 2]
  adjugate[2, 2] <- alpha[1, 1]
  
  adjugate[1, 2] <- -alpha[1, 2]
  adjugate[2, 1] <- -alpha[2, 1]
  
  
  ii <- inverse_det * adjugate
  
  return(ii)
  
  
}

#funciton to calculate the abundances of two species given the inverse of an interaction matrix and species growth rates
calculate_abundances <- function(r, inv_alpha){
  N1 <- (inv_alpha[1, 1] * r[1]) + (inv_alpha[1, 2] * r[2])
  N2 <- (inv_alpha[2, 1] * r[1]) + (inv_alpha[2, 2] * r[2])
  #check their feasibility
  N <- c(N1, N2)
  return(N)
}

#function to check if a point is inside the biological boundaries imposed by the model
check_r_boundaries <- function(r, rconstraints= NULL){
  #returns TRUE or FALSE
  if (is.null(rconstraints)) {
    return(TRUE)
  }else{
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
    return(all(r_good))
  }
  
}

# function to check if a set of abundances are below the maximum abundances expected per species
check_N_boundaries <- function(N, Nupper = NULL){
  #returns TRUE or FALSE
  if (is.null(Nupper)) {
    return(TRUE)
  }else{
    N_good <- (N <= Nupper) 
    return(all(N_good))
  }
}

#function to see if the point we generated is inside a radius
check_radius_boundaries <- function(r, R_max){
  #returns TRUE or FALSE
  radius <- sqrt(sum(r^2))
  
  return( radius <= R_max)
  
}

#function that for a given point, checks all of the growth rate boundaries, maximum abundances and feasibility.
check_point <- function(r,R_max,inv_alpha,rconstraints=NULL,Nupper=NULL){
  #returns NA if point is outside boundary
  #FALSE if it is inside boundary but infeasible
  #TRUE if it is inside boundary and feasible
  
  if(!check_radius_boundaries(r = r,
                         R_max = R_max)){
    return(NA)
  }
  
  
  
  if(!check_r_boundaries(r = r,
                         rconstraints = rconstraints)){
    return(NA)
  }
  
  #solve fo abundances
  N  <- calculate_abundances(r = r,
                            inv_alpha = inv_alpha)
  
  if(!check_N_boundaries(N = N,
                         Nupper = Nupper)){
    return(NA)
  }
  
  N_feasible <- (N > 0)
  N_feasible <- all(N_feasible)
  
  return(N_feasible)
}

#Function to generate a random point given the constraints
generate_point<- function(R_max){
  return(runif(2, -R_max, R_max))
}

#Function that generates many random points and checks our constraints and feasibility, by monte carlo sampling. 
integrate_area <- function(R_max,
                           alpha,
                           rconstraints=NULL,
                           Nupper=NULL,
                           n_samples){
  #returns the proportion of the area that is feasible, inside the boundaries and
  # the coordinates of the points that are feasible
  inv_alpha <- inverse_matrix(alpha)

  
  samples<- 0
  total <- 0
  coordinates_ri <- c()
  coordinates_rj <- c()
  coordinates_unfeasible_ri <- c()
  coordinates_unfeasible_rj <- c()
  while(samples < n_samples){
    #generate a random point inside a Radius defined by R_max
    point <- generate_point(R_max= R_max)
    #we check if it is within the boundaries, and if it is feasible
    result <- check_point(r = point,
                          R_max = R_max,
                          inv_alpha = inv_alpha,
                          rconstraints = rconstraints,
                          Nupper = Nupper)
    if(is.na(result)){
      next
    }
    
  #  col1<-ifelse(result,"blue","red")
  #  if the point is inside the boundary and feasible, we save its coordenates
     if(result){
       coordinates_ri <-c(coordinates_ri, point[1])
       coordinates_rj <-c(coordinates_rj, point[2])
     }else{
       coordinates_unfeasible_ri <- c(coordinates_unfeasible_ri, point[1])
       coordinates_unfeasible_rj <- c(coordinates_unfeasible_rj, point[2])
     }
    
   # points(point[1], point[2],col=col1, pch=20)
    
    samples <- samples + 1
    total <- total + result
  }
  
   ri_rj <- data.frame("ri"= coordinates_ri, "rj"= coordinates_rj)
   unfeasible <- data.frame("ri"= coordinates_unfeasible_ri, "rj"= coordinates_unfeasible_rj)
   proportion <- total/n_samples
   
  return(list("proportion"= proportion,
              "coords"= ri_rj,
              "unfeasible"= unfeasible))
  
}


#function to get the shape of the feasibility domain takes in a data frame of coordinates
determine_boundary_shape <- function(shape){
  nn <- nrow(shape)
  if(nn <= 1){
    shape_bounds <- data.frame("ri"=0, "rj"=0)
    return(list("bounds"=shape_bounds,
                "area"= 0))
      
  }else{
  
    #returns the coordinates of the limits
    # the area of the feasibility domain
    #we determine the boundary of the points, assuming it is a convex hull
    coord_points <- grDevices::chull(x = shape[,"ri"], y=shape[,"rj"]) 
    #we add the first point to complete the polygon
    coord_points<- c(coord_points, coord_points[1])
    shape_bounds <- shape[coord_points,]
    
    area_polygon <- sp::Polygon(coords = shape_bounds,
                                hole = FALSE)
    return(list("bounds"=shape_bounds,
                "area"= area_polygon@area))
  }
  

}



calculate_convex <-function(shape,
                            unfeasible){
  prop <- sp::point.in.polygon(point.x = unfeasible$ri,
                               point.y = unfeasible$rj,
                               pol.x = shape$ri,
                               pol.y = shape$rj)
return(  mean(prop))
  
}

#function to get the euclidean distance between to points
calculate_distance <- function(p1,p2){
  distance <- sqrt( (p2[1]-p1[1])^2 + (p2[2] - p1[2])^2   )
  return(distance)
}

#function that returns the shortest distance from the point defined by the growth rates and 
#the bounds of the feasibility domain
shortest_distance<-function(r, 
                            shape,
                            feasibility){
  #returns a minimum distance from the growth rates to a boundary of the feasibility domain
  # positive if it is inside the feasible
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


#Function to calculate the distance from limit fom the center of a feasibility domain and from the growth rates
distance_from_limit <- function(r,
                                shape,
                                feasibility){
  nn <- nrow(shape)
  #if there is no detectable shape
  if(nn <= 4){
    print("no feasibility domain detected")
    #and we detect the distance from our growth rates to the point 0,0
    distance <- calculate_distance( p1 =  r,
                                    p2 = c(0,0))
    distance <- -distance 
    results <- data.frame("center_distance" = distance,
                          "growth_distance"= 0)

  }else{

  
    
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

  }
  return(results)
  
}






