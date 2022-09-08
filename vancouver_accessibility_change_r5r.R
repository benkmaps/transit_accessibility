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
library(corrplot)

source("accessibility_functions.R")

# Path to data folder
data_path <- "data_van\\otp"

# Set up r5r network
r5r_core <- setup_r5(data_path, verbose = FALSE)

# Set up destinations table
test_area = read_sf("data_van\\stops_buffer.shp") %>% st_transform(32610)
land = read_sf("data_van\\gcma000b11a_e.shp") %>% st_transform(32610) %>% 
  filter(CMANAME == "Vancouver")
rb_sf = read_sf("data_van\\rapidbus.shp")

test_area = st_intersection(test_area, land %>% filter(CMANAME == "Vancouver")) %>% 
  dplyr::select(Id) # The test area is the intersection of land and the Vancouver CMA

# Set up destination points
dest_pts = read_sf("data_van\\core_poi_van_shp.shp") %>% st_transform(32610) # File of points representing POIs
dest_pts = dest_pts %>%
  st_intersection(test_area) %>% 
  rowid_to_column("dest_id") %>% 
  dplyr::select(dest_id) # Selecting only points that are within the test area

# Create Grid

cellsize = 300 # 300 unit grid cells (metres, since UTM 17N is being used)
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

# Count the number of destinations in each grid cell
dest_grid = grid %>%
  st_join(dest_pts) %>% 
  na.omit() %>% 
  group_by(id) %>% 
  summarise(dest_count = n()) %>% 
  st_transform(4326)

# Create an origin grid
origin_grid = grid %>% 
  st_intersection(test_area) %>% 
  dplyr::select(id,x) %>% 
  st_transform(4326)

# Find the centroids of the origin grid cells
dest_centroids = st_centroid(dest_grid)
dest_centroids$lat = st_coordinates(dest_centroids)[,"Y"]
dest_centroids$lon = st_coordinates(dest_centroids)[,"X"]

# Find the centroids of the destination grid cells
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

# Function for creating accessibility grid maps, specific to Vancouver
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

# Create accessibility grid
accessibility_grid = compare_accessibility(r5r_core,
                                           origins = origin_centroids,
                                           destinations = dest_centroids,
                                           datetimes = c("21-01-2020 7:00:00","17-09-2019 7:00:00"),
                                           time_window = 120)

# Save accessibility grid
write_sf(accessibility_grid, paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))
accessibility_grid = read_sf(paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))

accessibility_grid = read_sf(paste0("accessibility_grids\\van\\grid_van_",period_name,".gpkg"))

# Create standardized accessibility grid
accessibility_grid_strd = standardized_grid(accessibility_grid)

# Save standardized accessibility grid
write_sf(accessibility_grid_strd, paste0("accessibility_grids\\van\\grid_strd_van_",period_name,".gpkg"))
accessibility_grid_strd = read_sf(paste0("accessibility_grids\\van\\grid_strd_van_",period_name,".gpkg"))

# Create and save accessibility grid maps
tmap_save(create_accessibility_maps(accessibility_grid_strd),
          "exports\\maps\\accessibilityChangeStrd_van.png",
          width = 8.5,
          height = 8.5,
          units = "in")

# Correlation of accessibility values among impedance functions
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

# Significance of correlation of accessibility values among impedance functions
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
  cor.mtest(conf.level = 0.999)

