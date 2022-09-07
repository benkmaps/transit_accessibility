options(java.parameters = "-Xmx8G")

library(r5r)
library(sf)
library(data.table)
library(ggplot2)
library(akima)
library(dplyr)
library(bit64)
library(ggmap)
library(tidyverse)
library(tmap)
library(spdep)
library(ineq)
library(corrplot)

#setwd("C:\\Users\\bklar3\\Documents\\Transit_Accessibility\\transit_accessibility")
data_path <- "data_van\\otp" # Change to your data folder path

r5r_core <- setup_r5(data_path, verbose = FALSE)


# Set up destinations table
test_area = read_sf("data_van\\stops_buffer.shp") %>% st_transform(32610)
land = read_sf("data_van\\gcma000b11a_e.shp") %>% st_transform(32610) %>% 
  filter(CMANAME == "Vancouver")
rb_sf = read_sf("data_van\\rapidbus.shp")

test_area = st_intersection(test_area, land %>% filter(CMANAME == "Vancouver")) %>% 
  dplyr::select(Id)

#test_area = read_sf("data_r5r\\AdminBoundary.shp")[,"ShortName"] %>% 
#  st_transform(4326)

dest_pts = read_sf("data_van\\core_poi_van_shp.shp") %>% st_transform(32610)
dest_pts = dest_pts %>%
  st_intersection(test_area) %>% 
  rowid_to_column("dest_id") %>% 
  dplyr::select(dest_id)

# Create Grid

cellsize = 300
crs = 32610

grid = st_make_grid(test_area,
                    cellsize = cellsize,
                    what = "polygons",
                    square = FALSE
) %>% 
  st_as_sf() %>% 
  st_intersection(test_area) %>% 
  filter(st_geometry_type(.) %in% c("POLYGON","MULTIPOLYGON")) %>% 
  rowid_to_column("id")

dest_grid = grid %>%
  st_join(dest_pts) %>% 
  na.omit() %>% 
  group_by(id) %>% 
  summarise(dest_count = n()) %>% 
  st_transform(4326)

origin_grid = grid %>% 
  #st_join(test_area) %>% 
  #na.omit() %>%
  st_intersection(test_area) %>% 
  dplyr::select(id,x) %>% 
  st_transform(4326)

dest_centroids = st_centroid(dest_grid)
dest_centroids$lat = st_coordinates(dest_centroids)[,"Y"]
dest_centroids$lon = st_coordinates(dest_centroids)[,"X"]

origin_centroids = st_centroid(origin_grid)
origin_centroids$lat = st_coordinates(origin_centroids)[,"Y"]
origin_centroids$lon = st_coordinates(origin_centroids)[,"X"]

#Impedance Function Names
imp_fctn = paste0("dest_",
                  c("CUMR30",
                    "CUMR60",
                    "CUMR90",
                    "CUML30",
                    "CUML60",
                    "CUML90",
                    "NEXP_1",
                    "NEXP_2",
                    "NEXP_3",
                    "MGAUS_1",
                    "MGAUS_2",
                    "MGAUS_3"))

# Accessibility Calculation Functions

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
    
    #dplyr::select(!c((travel_time_p025:dest_count), lat:lon, x)) %>% 
    dplyr::select(c(fromId, dest_CUMR30:dest_MGAUS_3)) %>% 
    group_by(fromId) %>% 
    summarise_all(sum) %>% 
    
    left_join(origin_grid, by = c("fromId" = "id")) %>% 
    st_as_sf()
  
  return(accessibility)
}

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

# moransI_gini = function (accessibility_grid) {
#   
#   moransI_Gini = data.frame(imp_fctn = imp_fctn,
#                        I_bef = NA,
#                        I_aft = NA,
#                        I_delta = NA,
#                        Gini_bef = NA,
#                        Gini_aft = NA)
#   
#   for(n in 1:nrow(moransI_Gini)) {
#     
#     moransI_Gini[n,"I_bef"] = moran.test(accessibility_grid[,paste0(imp_fctn[n],"_bef")] %>% st_drop_geometry() %>% pull(),
#                        nb2listw(poly2nb(accessibility_grid),
#                                 zero.policy = TRUE),
#                        zero.policy = TRUE)[["estimate"]][["Moran I statistic"]]
#     
#     moransI_Gini[n,"I_aft"] = moran.test(accessibility_grid[,paste0(imp_fctn[n],"_aft")] %>% st_drop_geometry() %>% pull(),
#                        nb2listw(poly2nb(accessibility_grid),
#                                 zero.policy = TRUE),
#                        zero.policy = TRUE)[["estimate"]][["Moran I statistic"]]
#     
#     moransI_Gini[n,"I_delta"] = moran.test(accessibility_grid[,paste0(imp_fctn[n],"_delta")] %>% st_drop_geometry() %>% pull(),
#                                          nb2listw(poly2nb(accessibility_grid),
#                                                   zero.policy = TRUE),
#                                          zero.policy = TRUE)[["estimate"]][["Moran I statistic"]]
#     
#     
#     moransI_Gini[n,"Gini_bef"] = Gini(accessibility_grid[,paste0(imp_fctn[n],"_bef")] %>% st_drop_geometry() %>% pull())
#     moransI_Gini[n,"Gini_aft"] = Gini(accessibility_grid[,paste0(imp_fctn[n],"_aft")] %>% st_drop_geometry() %>% pull())
#     
#   }
#   return(moransI_Gini)
#   
# }

# lisa_map = function(x, v, var_name) {
#   
#   global = moran.test(v,
#                       nb2listw(poly2nb(x),
#                                zero.policy = TRUE),
#                       zero.policy = TRUE)[["estimate"]][["Moran I statistic"]]
#   
#   local = localmoran(v,
#                      nb2listw(poly2nb(x),
#                               zero.policy = TRUE),
#                      zero.policy = TRUE)
#   
#   quadrant <- vector(mode="numeric",length=nrow(local))
#   
#   # centers the variable of interest around its mean
#   m.qualification <- v - mean(v)     
#   
#   # centers the local Moran's around the mean
#   m.local <- local[,1] - mean(local[,1])    
#   
#   # significance threshold
#   signif <- 0.1 
#   
#   # builds a data quadrant
#   quadrant[m.qualification >0 & m.local>0] <- 4  
#   quadrant[m.qualification <0 & m.local<0] <- 1      
#   quadrant[m.qualification <0 & m.local>0] <- 2
#   quadrant[m.qualification >0 & m.local<0] <- 3
#   quadrant[local[,5]>signif] <- 0
#   
#   quadrant = factor(quadrant,
#                     levels = c(0, 1, 2, 3, 4),
#                     labels = c("Insignificant",
#                                "Low-Low",
#                                "Low-High",
#                                "High-Low",
#                                "High-High"),
#                     ordered = F)
#   
#   x = cbind(x, quadrant)
#   
#   colors <- c("grey90",
#               "blue",
#               rgb(0,0,1,alpha=0.4),
#               rgb(1,0,0,alpha=0.4),
#               "red")
#   
#   tm_shape(x) +
#     tm_fill(col = "quadrant",
#             palette = colors,
#             title = "LISA Clusters") +
#     tm_layout(title = paste0(sub("dest_","",var_name),"\n","I = ",round(global,2)),
#               title.position = c("right","bottom")) +
#     tm_shape(rb_sf) +
#     tm_lines("green4",
#              lwd = 3)
#   
#   
# }
# 
# create_lisa_maps = function(accessibility_grid, path) {
#   
#   for(x in imp_fctn) {
#     
#     map = lisa_map(accessibility_grid,
#                    accessibility_grid[paste0(x,"_delta")] %>% st_drop_geometry() %>% pull(),
#                    x)
#     
#     tmap_save(map, paste0(path,period_name,"_",x,".png"))
#     
#   }
# }

# lisa_map = function(x, v) {
#   global = moran.test(v,
#                       nb2listw(poly2nb(x),
#                                zero.policy = TRUE),
#                       zero.policy = TRUE)[["estimate"]][["Moran I statistic"]]
#   
#   local = localmoran(v,
#                      nb2listw(poly2nb(x),
#                               zero.policy = TRUE),
#                      zero.policy = TRUE)
#   
#   quadrant <- vector(mode="numeric",length=nrow(local))
#   
#   # centers the variable of interest around its mean
#   m.qualification <- v - mean(v)     
#   
#   # centers the local Moran's around the mean
#   m.local <- local[,1] - mean(local[,1])    
#   
#   # significance threshold
#   signif <- 0.1 
#   
#   # builds a data quadrant
#   quadrant[m.qualification >0 & m.local>0] <- 4  
#   quadrant[m.qualification <0 & m.local<0] <- 1      
#   quadrant[m.qualification <0 & m.local>0] <- 2
#   quadrant[m.qualification >0 & m.local<0] <- 3
#   quadrant[local[,5]>signif] <- 0
#   
#   quadrant = factor(quadrant,
#                     levels = c(0, 1, 2, 3, 4),
#                     labels = c("Insignificant",
#                                "Low-Low",
#                                "Low-High",
#                                "High-Low",
#                                "High-High"),
#                     ordered = F)
#   
#   x = cbind(x, quadrant)
#   
#   return(list(lisa = x, global = global))
# }
# 
# create_lisa_maps = function(accessibility_grid, path, var_name) {
# 
#   lisa_quads = data.frame(fromId = c(),
#                           quadrent = c(),
#                           fctn = c())
#   
#   global_Is = list()
#   
#   for(x in imp_fctn) {
#     
#     lisa = lisa_map(accessibility_grid,
#                    accessibility_grid[paste0(x,"_delta")] %>% st_drop_geometry() %>% pull())
#     
#     lisa_quads = rbind(lisa_quads,cbind(lisa$lisa[,c("fromId","quadrant")],fctn = sub("dest_","",x)))
#     
#     global_Is[[x]] = lisa$global
#     
#   }
#   
#   colors <- c("grey90",
#               "blue",
#               rgb(0,0,1,alpha=0.4),
#               rgb(1,0,0,alpha=0.4),
#               "red")
#   
#   map = tm_shape(lisa_quads) +
#     tm_fill(col = "quadrant",
#             palette = colors,
#             title = "LISA Clusters") +
#     tm_facets(by = "fctn") +
#     tm_shape(rb_sf) +
#     tm_lines("green4",
#              lwd = 3,
#              alpha = 0.6,
#              title.col = "RapidBus Lines",
#              legend.col.show = T)
#   
#   return(map)
# 
# }
# 
# create_combined_lisa_map = function(accessibility_grid) {
#   
#   maps = list()
#   
#   for(x in imp_fctn) {
#     
#     map = lisa_map(accessibility_grid,
#                    accessibility_grid[paste0(x,"_delta")] %>% st_drop_geometry() %>% pull(),
#                    x)
#     
#     maps[[imp_fctn]] = map
#     
#   }
#   
#   return(maps)
# }
# 
standardize_0 = function(x) {
  sd = sd(x)
  x = x/sd
  return(x)
}

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

create_accessibility_maps = function(grid_strd) {
  
  grid_strd %>% 
    dplyr::select(dest_CUMR30_delta : dest_MGAUS_3_delta) %>% 
    rename_with(substr, dest_CUMR30_delta : dest_MGAUS_3_delta, start = 6, stop = 12) %>% 
    rename(CUMR30 = CUMR30_,
           CUMR60 = CUMR60_,
           CUMR90 = CUMR90_,
           CUML30 = CUML30_,
           CUML60 = CUML60_,
           CUML90 = CUML90_,
           NEXP_1 = NEXP_1_,
           NEXP_2 = NEXP_2_,
           NEXP_3 = NEXP_3_) %>% 
    dplyr::select(CUMR30, CUMR60, CUMR90,
                   NEXP_1, NEXP_2, NEXP_3,
                   MGAUS_1, MGAUS_2, MGAUS_3,
                   CUML30, CUML60, CUML90) %>%
    
    pivot_longer(CUMR30:CUML90, names_to = "fctn", values_to = "value") %>% 
    st_as_sf() %>% 
    
    mutate(fctn = factor(fctn,
                         ordered = T,
                         levels = c("CUMR30", "CUMR60", "CUMR90",
                                    "NEXP_1", "NEXP_2", "NEXP_3",
                                    "MGAUS_1", "MGAUS_2", "MGAUS_3",
                                    "CUML30", "CUML60", "CUML90"))) %>% 
    
    tm_shape() +
    tm_fill(col = "value",
            title = "Standardized Accessibility Change Score",
            palette = "PRGn",
            midpoint = 0,
            breaks = c(-Inf, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, Inf),
            legend.is.portrait = F,
            legend.z = 1) +
    tm_facets(by = "fctn",
              ncol = 3,
              nrow = 4) +
      tm_shape(rb_sf) +
      tm_lines("orange",
               lwd = 1,
               legend.col.show = T,
               title.col = "RapidBus Lines",
               legend.col.z = 2,
               title.lwd = "RapidBus Lines") +
      tmap_mode("plot") +
    tm_layout(legend.outside.position = "bottom")
}


### END OF FUNCTIONS


# Routing Inputs
mode = c("WALK","TRANSIT")
max_walk_dist = 1600
max_trip_duration = 120

# Time Periods

# ----- SEPTEMBER 2019 to JANUARY 2020
#       Weekdays 7:00 to 9:00

period_name = "AMpeak"

accessibility_grid = compare_accessibility(r5r_core,
                                           origins = origin_centroids,
                                           destinations = dest_centroids,
                                           datetimes = c("21-01-2020 7:00:00","17-09-2019 7:00:00"),
                                           time_window = 120)

write_sf(accessibility_grid, paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))
accessibility_grid = read_sf(paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))

accessibility_grid = read_sf(paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))

accessibility_grid_strd = standardized_grid(accessibility_grid)

write_sf(accessibility_grid_strd, paste0("accessibility_grids\\van\\grid_strd_van_",period_name,".gpkg"))
accessibility_grid_strd = read_sf(paste0("accessibility_grids\\van\\grid_strd_van_",period_name,".gpkg"))

#write_csv(moransI_gini(accessibility_grid), paste0("exports\\moran_gini\\van\\moranGini_van_",period_name,".csv"))

# create_lisa_maps(accessibility_grid, "exports\\lisa_maps\\van\\LISA_van_")
# tmap_save(create_lisa_maps(accessibility_grid),
#           "exports\\lisa_maps\\LISA_van.png",
#           width = 8.5,
#           height = 11,
#           units = 'in')

tmap_save(create_accessibility_maps(accessibility_grid_strd),
          "exports\\maps\\accessibilityChangeStrd_van.png",
          width = 8.5,
          height = 8.5,
          units = "in")

accessibility_grid %>% 
  st_drop_geometry() %>% 
  dplyr::select(dest_CUMR30_delta:dest_MGAUS_3_delta) %>% 
  rename_with(substr, start = 6, stop = 12) %>% 
  rename(CUMR30 = CUMR30_,
              CUMR60 = CUMR60_,
              CUMR90 = CUMR90_,
              CUML30 = CUML30_,
              CUML60 = CUML60_,
              CUML90 = CUML90_,
              NEXP_1 = NEXP_1_,
              NEXP_2 = NEXP_2_,
              NEXP_3 = NEXP_3_) %>% 
  dplyr::select(CUMR30, CUMR60, CUMR90,
                NEXP_1, NEXP_2, NEXP_3,
                MGAUS_1, MGAUS_2, MGAUS_3,
                CUML30, CUML60, CUML90) %>% 
  cor() %>%
  corrplot(type = "lower",
           addCoef.col = "black")

