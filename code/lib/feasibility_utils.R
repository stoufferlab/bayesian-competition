
# we use functions from these packages below
require("grDevices")
require("polylabelr")
require("sp")
require("SpatialGraph")

#FUNCTIONS TO GET THE SHAPE AND AREA OF THE FEASIBILITY DOMAIN, AS WELL AS THE SHORTEST DISTANCE FROM AN EDGE OF THE GROWTH RATES

# "fast" function to get the inverse of a 2x2 matrix
inverse_matrix<-function(alpha){
  
  determinant_alpha <- (alpha[1, 1] * alpha[2, 2]) - (alpha[2, 1] * alpha[1, 2])
  # inverse_det  <- 1 / determinant_alpha

  adjugate <- matrix(NA, nrow = 2, ncol = 2)
  adjugate[1, 1] <- alpha[2, 2]
  adjugate[2, 2] <- alpha[1, 1]
  adjugate[1, 2] <- -alpha[1, 2]
  adjugate[2, 1] <- -alpha[2, 1]
  
  ialpha <- adjugate / determinant_alpha
  
  return(ialpha)
}

# function to calculate the equilibrium abundances of two species given the inverse of an interaction matrix and species growth rates
calculate_abundances <- function(r, inv_alpha){
  return(t(inv_alpha %*% r))
}

#function to check if a vector of growth rates is inside the biological boundaries imposed by the model
#returns TRUE or FALSE
within_r_boundaries <- function(r, rconstraints = NULL){
  if (is.null(rconstraints)) {
    return(TRUE)
  }else{
    r_good <- (r >= rconstraints$lower) & (r <= rconstraints$upper)
    return(all(r_good))
  }
}

#function to check if a set of abundances are below the maximum abundances expected per species
#returns TRUE or FALSE
within_N_boundaries <- function(r, inv_alpha, gN_max = NULL){
  if (is.null(gN_max)) {
    return(TRUE)
  }else{
    # solve for equilibrium abundances
    Nequil <- calculate_abundances(r = r, inv_alpha = inv_alpha)
    return(all(Nequil <= gN_max))
  }
}

# given interaction coefficients and abundance constraints
# determine the maximum permissable R value
maximum_radius <-function(alpha, gN_max){
  # separate the maximum values
  giNi_max <- gN_max[1]
  gjNj_max <- gN_max[2]

  # giNi = 0 and gjNj = 0
  R1 <- 0

  # giNi = giNi_max and gjNj = 0
  R2 <- sqrt((alpha[1,1]*giNi_max)^2 + (alpha[2,1]*giNi_max)^2)

  # giNi = 0 and gjNj = gjNj_max
  R3 <- sqrt((alpha[1,2]*gjNj_max)^2 + (alpha[2,2]*gjNj_max)^2)

  # giNi = giNi_max & gjNj = gjNj_max
  R4 <- sqrt((alpha[1,1]*giNi_max + alpha[1,2]*gjNj_max)^2 + (alpha[2,1]*giNi_max + alpha[2,2]*gjNj_max)^2)
  
  # 0 < giNi < giNi_max & gjNj = 0
  R5 <- 0

  # giNi = 0 & 0 < gjNj < gjNj_max
  R6 <- 0

  # 0 < Ni < giNi_max & gjNj = gjNj_max
  Ni <- -gjNj_max*(alpha[1,1]*alpha[1,2] + alpha[2,1]*alpha[2,2])/(alpha[1,1]^2 + alpha[2,1]^2)
  R7 <- sqrt((alpha[1,1]*Ni + alpha[1,2]*gjNj_max)^2 + (alpha[2,1]*Ni + alpha[2,2]*gjNj_max)^2)

  # giNi = giNi_max & 0 < gjNj < gjNj_max
  gjNj <- -giNi_max*(alpha[1,1]*alpha[1,2] + alpha[2,1]*alpha[2,2])/(alpha[1,2]^2 + alpha[2,2]^2)
  R8 <- sqrt((alpha[1,1]*giNi_max + alpha[1,2]*gjNj)^2 + (alpha[2,1]*giNi_max + alpha[2,2]*gjNj)^2)
  
  # 0 < giNi < giNi_max & 0 < gjNj < gjNj_max
  R9 <- 0

  # set of all values
  Rvalues <- c(R1, R2, R3, R4, R5, R6, R7, R8, R9)
  
  # return the maximum R
  return(max(Rvalues))
}


#function to see if the vector of growth rates is inside a radius
#returns TRUE or FALSE
within_radius_boundaries <- function(r, R_max = 1){
  radius <- sqrt(sum(r^2))
  return(radius <= R_max)
}

#function to see if the vector of growth rates and inverse alpha matrix leads to feasible equilibrium
#returns TRUE or FALSE
is_feasible <- function(r, inv_alpha){
  Nequil <- calculate_abundances(r = r, inv_alpha = inv_alpha)
  return(all(Nequil > 0))
}

# # TODO: need a function that determines which species competitively excludes which other species
# # returns TRUE or FALSE
# can_invade <- function(r, alpha){
#   Nequil_resident <- r / diag(alpha)
#   Nequil <- calculate_abundances(r = r, inv_alpha = inv_alpha)
#   return(all(Nequil > 0))
# }

#function that for a given vector of growth rates checks all of:
# (i) the growth rate boundaries
# (ii) maximum abundances
# (iii) equilibrium feasibility
#returns:
# NA if point is outside boundary
# FALSE if it is inside boundary but infeasible
# TRUE if it is inside boundary and feasible
check_point <- function(r,inv_alpha,R_max,rconstraints=NULL,gN_max=NULL){
  feasibility <- is_feasible(r = r, inv_alpha = inv_alpha)
  if(!feasibility)
    return(feasibility)
  else{
    faults <- NULL
    if(!within_radius_boundaries(r = r, R_max = R_max))
      faults <- c(faults, "RMAX")

    if(!within_r_boundaries(r = r, rconstraints = rconstraints))
      faults <- c(faults, "MODEL")

    if(!within_N_boundaries(r = r, inv_alpha = inv_alpha, gN_max = gN_max))
      faults <- c(faults, "ABUNDANCE")

    if(is.null(faults)){
      return(feasibility)
    }else{
      return(paste0(faults, collapse="_"))
    }
  }
}

#Function to generate a random point within a circle of radius R_max
generate_point <- function(R_max){
  a <- runif(1) * 2 * pi
  r = R_max * sqrt(runif(1))
  return(c(r*cos(a), r*sin(a)))
}

# Function that generates many random points and checks our constraints and feasibility by monte carlo sampling
# desired feasible is the number of points we want in the feasibility region
# max_samples is the max number of samples before it gives up
integrate_area <- function(
    alpha,
    R_max,
    rconstraints=NULL,
    gN_max=NULL,
    desired_feasible=1e3,
    max_samples=1e6
){
  # returns:
  # the coordinates of the points that are feasible and inside the constraints
  # the coordinates of the points that are not
  # the proportion of the points that are feasible and inside the constraints

  inv_alpha <- inverse_matrix(alpha)
  
  samples<- 0
  total <- 0
  coordinates_ri <- c()
  coordinates_rj <- c()
  coordinates_unfeasible_ri <- c()
  coordinates_unfeasible_rj <- c()
  
  while(total < desired_feasible & samples < max_samples){
    #generate a random point inside a radius defined by R_max
    point <- generate_point(R_max = R_max)
    #we check if it is within the boundaries, and if it is feasible
    result <- check_point(
      r = point,
      inv_alpha = inv_alpha,
      R_max = R_max,
      rconstraints = rconstraints,
      gN_max = gN_max
    )
    if(!is.logical(result)){
      next
    }
    
    #  if the point is inside the boundary and feasible, we save its coordinates
    if(result){
      coordinates_ri <-c(coordinates_ri, point[1])
      coordinates_rj <-c(coordinates_rj, point[2])
    }else{
      coordinates_unfeasible_ri <- c(coordinates_unfeasible_ri, point[1])
      coordinates_unfeasible_rj <- c(coordinates_unfeasible_rj, point[2])
    }
    
    samples <- samples + 1
    total <- total + result
  }
  
  ri_rj <- data.frame("ri"= coordinates_ri, "rj"= coordinates_rj)
  unfeasible <- data.frame("ri"= coordinates_unfeasible_ri, "rj"= coordinates_unfeasible_rj)
  proportion <- total/samples
  
  return(list("coords"= ri_rj, "unfeasible"=unfeasible, proportion = proportion))
}


# function to convert a set of points into their convex hull and the area thereof
# takes
# shape: a set of x,y coordinates as input
# returns
# bounds: the x,y coordinates of the limits of the convex hull
# area: the area of the convex hull
determine_boundary_shape <- function(shape){
  if(nrow(shape) <= 4){
    shape_bounds <- data.frame("ri"=0, "rj"=0)
    shape_area <- 0
  }else{
    # determine the boundary of the points assuming it is a convex hull
    coord_points <- grDevices::chull(x = shape[,"ri"], y=shape[,"rj"]) 
    #we add the first point to complete the polygon
    coord_points<- c(coord_points, coord_points[1])
    # and use the points to select only boundary points (in order)
    shape_bounds <- shape[coord_points,]
    
    shape_area <- sp::Polygon(
      coords = shape_bounds,
      hole = FALSE
    )@area
  }
  return(list(bounds=shape_bounds,area=shape_area))
}

# determine if any of a sample of points is inside the feasiblity domain
# this helps establish the validity of treating the domain as a convex hull
convexity_check <-function(
  pts,
  shape
){
  pips <- sp::point.in.polygon(
    point.x = pts$ri,
    point.y = pts$rj,
    pol.x = shape$ri,
    pol.y = shape$rj
  )
  return(sum(pips == 0) / length(pips))
}

#function to get the euclidean distance between to points
calculate_distance <- function(p1,p2){
  distance <- sqrt( (p2[1]-p1[1])^2 + (p2[2] - p1[2])^2   )
  return(distance)
}

#function that returns the shortest distance from the point defined by the growth rates and 
#the bounds of the feasibility domain
shortest_distance<-function(
  r, 
  shape
){
  #returns a minimum distance from the growth rates to a boundary of the feasibility domain

  # make r a matrix to use a function below
  r <- t(as.matrix(r))

  #for every line defined by two points, we calculate the distance of our growth rates to it 
  data_dist <- c()
  for(i in seq.int(nrow(shape)-1)){
    # points that define the line
    p0 <- c(shape[i, ]$ri, shape[i, ]$rj)
    p1 <- c(shape[i + 1, ]$ri, shape[i + 1, ]$rj)
    
    # convert the line into a matrix
    line_matrix <- as.matrix(rbind(p0, p1))
    
    #we get the shortest distance from the point to the line, which is the distance between the point and the perpendicular projections of the line
    dist <- SpatialGraph::pointLineD(
      xy = line_matrix,
      xyp = r
    )
    # does the perpendicular projection of the points cross the segment or not?
    cross <- dist$cross
    #if it does not, then the shortest distance is the distance to one of the ends of the line
    if (cross == 0) {
      dist_p0 <- calculate_distance(p1 = p0,
                                    p2 = r)
      
      dist_p1 <- calculate_distance(p1 = p1,
                                    p2 = r)
      dist_m <-  min(dist_p0, dist_p1)
      type <- "point"
    } else{
      #if it does, then the shortest distance is to the perpendicular projection of the line
      dist_m <- dist$d
      type <- "line"
    }
    
    data <- data.frame(
      "distance_from_edge" = dist_m,
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
  #but it has to be ONLY ONE
  minimum_distance <- unique(minimum_distance$distance_from_edge)

  return(minimum_distance)
}


# calculate the distance from (i) the center of a feasibility domain and (ii) from the growth rates to the edge of the feasbility domain
distance_from_limit <- function(
  r,
  shape
){
  #if there is no detectable shape
  if(nrow(shape) <= 4){
    warning("no feasibility domain detected")
    #and we detect the distance from our growth rates to the point 0,0
    distance <- calculate_distance(
      p1 =  r,
      p2 = c(0,0)
    )
    results <- data.frame("center_distance" = distance,
                          "growth_distance" = 0, # DEBUG: shouldn't this be "distance"?
                          "detection" = FALSE)

  }else{
    #the center of the polygon
    center <- polylabelr::poi(
      x = shape$ri,
      y = shape$rj
    )
    
    #we get the shortest distance from the center of the polygon to an edge, it is always feasible
    center_distance <- shortest_distance(
      r = c(center$x, center$y),
      shape = shape
    )
    
    #we get the shortest distance from our growth rates to an edge
    growth_distance <- shortest_distance(
      r = r,
      shape = shape
    )
    
    results <- data.frame(
      "center_distance" = center_distance,
      "growth_distance" = growth_distance,
      "detection" = TRUE
    )

  }
  return(results)
  
}

# functions to get the area of species without competition
get_boundary_r <-function(intraspecific_competition,
                          gN_max,
                          lower,
                          upper){
  max_growth_rate <- intraspecific_competition * gN_max
  
  #if there is facilitation, then the area alone should include negative growth rates, 
  #but still be limited by model constraints
  if(intraspecific_competition < 0){
    bounds <- max(max_growth_rate, lower)
  }else{
    bounds <- min(max_growth_rate, upper)
  }
  
  return(bounds)
}

# size of growth rate parameter space that allows the two species to survive if they don't interact with each other
area_species_alone <- function(alpha,
                               gN_max,
                               rconstraints){
  # min/max of ri that allows monoculture equilibrium
  ri_bound <- get_boundary_r(intraspecific_competition = alpha[1,1],
                             gN_max = gN_max[1],
                             lower = rconstraints$lower[1],
                             upper = rconstraints$upper[1])
  
  # min/max of rj that allows monoculture equilibrium
  rj_bound <- get_boundary_r(intraspecific_competition = alpha[2,2],
                             gN_max = gN_max[2],
                             lower = rconstraints$lower[2],
                             upper = rconstraints$upper[2])
  
  # monoculture domain is always a rectangle
  area_alone <- abs(ri_bound * rj_bound)
  
  return(area_alone)
}

# given a vector of growth rates, an alpha matrix, and a variety of constraints
# perform a series of checks regarding the size of the feasibility domain, etc
structural_stability_wrapper <- function(
  r,
  alpha,
  R_max,
  gN_max,
  rconstraints,
  desired_feasible = 2000,
  max_samples = 2E5
){

  #we do a mcmc integration of the feasible area
  integration<- integrate_area(
    alpha = alpha,
    R_max = R_max,
    rconstraints = rconstraints,
    gN_max = gN_max,
    desired_feasible = desired_feasible,
    max_samples = max_samples
  )

  # use the feasible points to define the convex hull of the feasibility domain
  shape <- determine_boundary_shape(shape = integration$coords)

  # these points define the limits of the convex hull
  bounds <- shape$bounds
  # we also calculate the area of the convex hull
  area <- shape$area
  
  #with the bounds we can then get the distance from the limit of our growth rates
  distances <- distance_from_limit(
    r = r,
    shape = bounds
  )
  center_to_edge <- distances$center_distance
  observed_to_edge <- distances$growth_distance

  #calculate the proportion of unfeasible points  things outside the convex hull
  domain_convex <- convexity_check(
    pts = integration$unfeasible,
    shape = bounds
  )

  # check if the inferred growth rates correspond to a biologically constrained feasible equilibrium
  biological_feasibility <- check_point(
    r = r,
    R_max = R_max,
    inv_alpha = inverse_matrix(alpha),
    rconstraints = rconstraints,
    gN_max = gN_max
  )

  # calculate the size of the region that allows both to exist in monoculture
  area_alone <- area_species_alone(
    alpha = alpha,
    gN_max = gN_max,
    rconstraints = rconstraints
  )

  # calculate the monoculture equilibria
  Nmonoc <- r / diag(alpha)
  Nmonoc[Nmonoc<0] <- NA

  # calculate the location of the multispecies equilibrium
  Nequil <- calculate_abundances(r, inverse_matrix(alpha))
  feasibility <- all(Nequil > 0)

  # DEBUG: determine which species "wins" when outside the two-species feasibility domain
  if(!feasibility){
    if(Nequil[1]>0 & Nequil[2]<0)
      coexist <- "i"
    else if(Nequil[1]<0 & Nequil[2]>0)
      coexist <- "j"
    else if(Nequil[1]<0 & Nequil[2]<0)
      coexist <- NA
  }else
    coexist <- "BOTH"
  # # TODO: the above needs to be checked more rigorously
  # coexist <- NA

  # stitch all results together in a dataframe
  results <- data.frame(
    area_feasible = area,
    area_alone = area_alone,
    proportion = area / area_alone,
    center_to_edge = center_to_edge,
    observed_to_edge = observed_to_edge,
    observed_to_edge_signed = ifelse(biological_feasibility==TRUE, observed_to_edge, -observed_to_edge),
    domain_detected = distances$detection,
    domain_convex = domain_convex,
    giNi_mono = Nmonoc[1],
    gjNj_mono = Nmonoc[2],
    giNi_equil = Nequil[1],
    gjNj_equil = Nequil[2],
    coexist = coexist,
    feasibility = feasibility,
    biological_feasibility = biological_feasibility
  )

  return(results)
}
