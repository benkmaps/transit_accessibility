# Function that computes  accessibility scores for one GTFS time period according to the 12 specified impedance functions
calculate_accessibility = function(ttm) {
  
  mgaus = data.frame(t = c(25,50,75),
                     w = 0.5)
  mgaus = mgaus %>% mutate(beta = -(t^2)/log(w))
  
  nexp = data.frame(t = c(15,30,45),
                    w = c(0.5))
  nexp = nexp %>% mutate(beta = -log(w)/t)
  
  
  accessibility = ttm %>% 
    mutate(toId = as.numeric(toId)) %>% 
    mutate(fromId = as.numeric(fromId)) %>% 
    filter(travel_time < 120) %>% 
    inner_join(dest_centroids, by = c("toId" = "id")) %>% 
    
    mutate(dest_CUMR30 = dest_count * ifelse(travel_time <= 30, 1, 0)) %>% 
    mutate(dest_CUMR60 = dest_count * ifelse(travel_time <= 60, 1, 0)) %>% 
    mutate(dest_CUMR90 = dest_count * ifelse(travel_time <= 90, 1, 0)) %>% 
    
    mutate(dest_CUML30 = dest_count * (30-travel_time) / 30 * ifelse(travel_time <= 30, 1, 0)) %>% 
    mutate(dest_CUML60 = dest_count * (60-travel_time) / 60 * ifelse(travel_time <= 60, 1, 0)) %>% 
    mutate(dest_CUML90 = dest_count * (90-travel_time) / 90 * ifelse(travel_time <= 90, 1, 0)) %>% 
    
    mutate(dest_NEXP_1 = dest_count * exp(-nexp[1,"beta"] * travel_time)) %>% 
    mutate(dest_NEXP_2 = dest_count * exp(-nexp[2,"beta"] * travel_time)) %>% 
    mutate(dest_NEXP_3 = dest_count * exp(-nexp[3,"beta"] * travel_time)) %>% 
    
    mutate(dest_MGAUS_1 = dest_count * exp(-1 * travel_time^2 / mgaus[1,"beta"])) %>% 
    mutate(dest_MGAUS_2 = dest_count * exp(-1 * travel_time^2 / mgaus[2,"beta"])) %>% 
    mutate(dest_MGAUS_3 = dest_count * exp(-1 * travel_time^2 / mgaus[3,"beta"])) %>% 
    
    dplyr::select(c(fromId, dest_CUMR30:dest_MGAUS_3)) %>% 
    group_by(fromId) %>% 
    summarise_all(sum) %>% 
    
    left_join(origin_grid, by = c("fromId" = "id")) %>% 
    st_as_sf()
  
  return(accessibility)
}

# A function that computes the difference between two accessibility scores by impedance function
calculate_accessibility_change = function(after, before) {
  change = inner_join(st_drop_geometry(after),
                      st_drop_geometry(before),
                      by = "fromId",
                      suffix = c("_aft","_bef")) %>% 
    mutate(dest_CUMR30_delta = dest_CUMR30_aft - dest_CUMR30_bef) %>% 
    mutate(dest_CUMR60_delta = dest_CUMR60_aft - dest_CUMR60_bef) %>% 
    mutate(dest_CUMR90_delta = dest_CUMR90_aft - dest_CUMR90_bef) %>% 
    
    mutate(dest_CUML30_delta = dest_CUML30_aft - dest_CUML30_bef) %>% 
    mutate(dest_CUML60_delta = dest_CUML60_aft - dest_CUML60_bef) %>% 
    mutate(dest_CUML90_delta = dest_CUML90_aft - dest_CUML90_bef) %>% 
    
    mutate(dest_NEXP_1_delta = dest_NEXP_1_aft - dest_NEXP_1_bef) %>% 
    mutate(dest_NEXP_2_delta = dest_NEXP_2_aft - dest_NEXP_2_bef) %>% 
    mutate(dest_NEXP_3_delta = dest_NEXP_3_aft - dest_NEXP_3_bef) %>% 
    
    mutate(dest_MGAUS_1_delta = dest_MGAUS_1_aft - dest_MGAUS_1_bef) %>% 
    mutate(dest_MGAUS_2_delta = dest_MGAUS_2_aft - dest_MGAUS_2_bef) %>% 
    mutate(dest_MGAUS_3_delta = dest_MGAUS_3_aft - dest_MGAUS_3_bef) %>% 
    
    inner_join(grid, by = c("fromId" = "id")) %>% 
    st_as_sf()
  
  return(change)
}

# A function that calculates travel time matrices for a before and after GTFS time period
# and uses calculate_accessibility and calculate_accessibility_change to create a table with
# accessibility before, after, and delta values for each impedance function
compare_accessibility = function(r5r_core,
                                 origins,
                                 destinations,
                                 datetimes,
                                 time_window)
{
  # Routing Inputs
  mode = c("WALK","TRANSIT")
  max_walk_dist = 1600
  max_trip_duration = 120
  
  departure_datetime = c(as.POSIXct(datetimes[1],
                                    format = "%d-%m-%Y %H:%M:%S"),
                         as.POSIXct(datetimes[2],
                                    format = "%d-%m-%Y %H:%M:%S"))
  
  # Travel Time Matrices
  
  ttm_after = travel_time_matrix(r5r_core,
                                 origins = origins,
                                 destinations = destinations,
                                 mode = mode,
                                 departure_datetime = departure_datetime[1],
                                 max_walk_dist = max_walk_dist,
                                 max_trip_duration = max_trip_duration,
                                 time_window = time_window,
                                 verbose = F)
  
  ttm_before = travel_time_matrix(r5r_core,
                                  origins = origins,
                                  destinations = destinations,
                                  mode = mode,
                                  departure_datetime = departure_datetime[2],
                                  max_walk_dist = max_walk_dist,
                                  max_trip_duration = max_trip_duration,
                                  time_window = time_window,
                                  verbose = F)
  
  # Calculate Accessibility
  
  accessibility_grid_after = calculate_accessibility(ttm_after)
  accessibility_grid_before = calculate_accessibility(ttm_before)
  
  accessibility_grid_change = calculate_accessibility_change(accessibility_grid_after,
                                                             accessibility_grid_before)
  
  return(accessibility_grid_change) 
}

# Standardizes values by dividing by the standard deviation,
# WITHOUT centering on the mean
standardize_0 = function(x) {
  sd = sd(x)
  x = x/sd
  return(x)
}

# Standardized the accessibility grid created from compare_accessibility using standarzide_0
standardized_grid = function(grid) {
  grid_strd = grid %>% 
    st_drop_geometry() %>% 
    dplyr::select(dest_CUMR30_aft:dest_MGAUS_3_delta) %>% 
    sapply(standardize_0) %>% 
    as.data.frame() %>% 
    cbind(grid[,"fromId"]) %>% 
    rowwise(fromId) %>% 
    mutate(sd = sd(c_across(dest_CUMR30_aft:dest_MGAUS_3_delta))) %>% 
    st_as_sf()
  
  return(grid_strd)
}
