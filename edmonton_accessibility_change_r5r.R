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
data_path <- "data_edm\\otp"

# Set up r5r network
r5r_core <- setup_r5(data_path, verbose = FALSE)


# Set up destinations table
test_area = read_sf("data_edm\\edm_boundary.shp") %>% st_transform(32612)

# Read bus lines shapefile
buslines_after = read_sf("data_edm\\otp\\gtfs_edm_20210504\\shapes.shp")
buslines_before = read_sf("data_edm\\otp\\gtfs_edm_20210504\\shapes_0222.shp")

# Set up destination points
dest_pts = read_sf("data_edm\\core_poi_edm_shp.shp") %>% st_transform(32612) # File of points representing POIs
dest_pts = dest_pts %>% 
  st_intersection(test_area) %>% 
  rowid_to_column("dest_id") %>% 
  dplyr::select(dest_id) # Selecting only points that are within the test area

# Create Grid

cellsize = 300 # 300 unit grid cells (metres, since UTM 17N is being used)
crs = 32612

grid = st_make_grid(test_area,
                    cellsize = cellsize,
                    crs = crs,
                    what = "polygons",
                    square = FALSE
) %>% 
  st_as_sf() %>% 
  st_intersection(test_area) %>% 
  filter(st_geometry_type(.) %in% c("POLYGON","MULTIPOLYGON")) %>% 
  dplyr::select(x) %>% 
  rowid_to_column("id")

# Count the number of destinations in each grid cell
dest_grid = grid %>%
  st_join(dest_pts) %>% 
  na.omit() %>% 
  group_by(id) %>% 
  summarise(dest_count = n()) %>% 
  st_transform(4326)

# Ce=reate an origin grid
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

# Function for creating accessibility grid maps, specific to Edmonton
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
            legend.is.portrait = F) +
    tm_facets(by = "fctn", ncol = 3, nrow = 4) +
    
    tm_shape(buslines_after) +
    tm_lines("orange",
             lwd = 0.25,
             alpha = 0.5,
             legend.col.show = T,
             title.col = "Bus Lines - April 2021",
             legend.col.z = 2) +

    # tm_shape(buslines_before) +
    # tm_lines("blue",
    #          lwd = 0.25,
    #          alpha = 0.5,
    #          legend.col.show = T,
    #          title.col = "Bus Lines - February 2021",
    #          legend.col.z = 2) +
    
    tm_layout(legend.outside.position = "bottom") +
    
    tmap_mode("plot")
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
                             datetimes = c("11-05-2021 7:00:00","06-04-2021 7:00:00"),
                             time_window = 120)

# Save accessibility grid
write_sf(accessibility_grid, paste0("accessibility_grids\\edm\\grid_edm_",period_name,".gpkg"))
accessibility_grid = read_sf(paste0("accessibility_grids\\edm\\grid_edm_",period_name,".gpkg"))

# Create standardized accessibility grid
accessibility_grid_strd = standardized_grid(accessibility_grid)

# Save standardized accessibility grid
write_sf(accessibility_grid_strd, paste0("accessibility_grids\\edm\\grid_strd_edm_",period_name,".gpkg"))
accessibility_grid_strd = read_sf(paste0("accessibility_grids\\edm\\grid_strd_edm_",period_name,".gpkg"))

# Create and save accessibility grid maps
tmap_save(create_accessibility_maps(accessibility_grid_strd),
          "exports\\maps\\accessibilityChangeStrd_edm.png",
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

