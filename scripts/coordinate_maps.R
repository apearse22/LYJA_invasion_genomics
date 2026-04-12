##############################################

# creating map of sample locations 
# code developed by Abby Pearse, Jessie Pelosi
# last updated: 04/07/26

##############################################

### libraries

library(ggplot2)
library(dplyr)
library(ggmap)
library(maps)
library(readr)
library(gridExtra)
library(ggspatial)
library(terra)
library(sf)

######################################################### creating basemaps #####################################################

### invaded range

states_map <- map_data("state")
invaded <- states_map %>% 
  filter(region == "alabama" |
           region == "georgia" |
           region == "south carolina" |
           region == "florida" |
           region == "arkansas" |
           region == "texas" |
           region == "mississippi" |
           region == "louisiana" |
           region == "oklahoma") # no samples in oklahoma, just for aesthetic purposes

invaded_basemap <- ggplot(invaded, aes(long, lat, group=group)) +
  geom_polygon(fill='gray97', color='black') +
  coord_fixed(1.3) +
  theme_classic() +
  xlab('Latitude (\u00B0)') +
  ylab('Longitude (\u00B0)') 


### native

country_map <- map_data("world")
native <- country_map %>% 
  filter(region == "China" |
           region == "Japan" |
           region == "Vietnam" |
           region == "Philippines" |
           region == "Palau" |
           region == "Taiwan" |
           region == "Indonesia")

native_basemap <- ggplot(native, aes(long, lat, group=group)) +
  geom_polygon(fill='gray97', color='black') +
  coord_fixed(1.3) +
  theme_classic() +
  xlab('Latitude (\u00B0)') +
  ylab('Longitude (\u00B0)') 



################ playing around with doing the whole continent and subsetting by xlim and ylim

us_basemap <- ggplot(states_map, aes(long, lat, group=group)) +
  geom_polygon(fill='white', color='black') +
  coord_fixed(1.3) +
  theme_classic() +
  ylab('Latitude (\u00B0)') +
  xlab('Longitude (\u00B0)') +
  coord_sf(xlim = c(-96, -79), ylim = c(25, 36.5)) 
  #annotation_scale(location = "bl", unit_category = "metric", width_hint = 0.1, plot_unit = "km") 
  #scalebar(location = "bottomleft", dist_unit = "km", transform = TRUE,
         #  model = "WGS84")

############################################################ plotting coordinates #####################################################


### getting coordinates from supplmental information table 


si_table <- read.csv("files/popmap_updatedcoords.csv")

coordinates <- si_table %>% 
  dplyr::select(Approx..Latitude, Approx..Longitude)

#coordinates <- read.csv("files/lyja_si_coordinates.csv")

si_table <- si_table %>% 
  mutate(Approx..Latitude = as.numeric(Approx..Latitude)) %>% 
  mutate(Approx..Longitude = as.numeric(Approx..Longitude))

invaded_coordinates <- si_table %>% 
  filter(Invaded == "Y")

native_coordinates <- si_table %>% 
  filter(Invaded == "N") %>% 
  na.omit()


# CONVERT TO 'sf' OBJECTS: 
# This is the magic step that fixes your scale bar. It tells R these are 
# geographic coordinates (crs = 4326, which is WGS84).
invaded_sf <- st_as_sf(invaded_coordinates, coords = c("Approx..Longitude", "Approx..Latitude"), crs = 4326)

# Assuming native_coordinates keeps the original column names
native_sf <- st_as_sf(native_coordinates, coords = c("Approx..Longitude", "Approx..Latitude"), crs = 4326)

### plotting invaded coordinates

######################################################### 
# 2. Creating Basemaps & Plotting 
#########################################################

# Convert map boundaries directly to sf objects instead of map_data dataframes
states_map_sf <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
country_map_sf <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))

# Find the global min and max years across both datasets
global_year_limits <- range(c(invaded_sf$Collection.Year, native_sf$Collection.Year), na.rm = TRUE)

### PLOTTING INVADED RANGE ###

invaded_coordinates_mapped <- ggplot(data = states_map_sf) +
  geom_sf(fill = 'gray97', color = 'black') +
  geom_sf(data = invaded_sf, aes(color = Collection.Year), size = 3) +
  scale_color_viridis_c(limits = global_year_limits) +
  labs(color = "Collection Year", x = 'Longitude (\u00B0)', y = 'Latitude (\u00B0)') +
  theme_classic() +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-96, -79), ylim = c(25, 37), crs = 4326, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)

# Save Invaded Map
ggsave("invaded_sample_coords.pdf", plot = invaded_coordinates_mapped, width = 8, height = 5.5)


### PLOTTING NATIVE RANGE ###

# Filter the world map for native regions
native_basemap_sf <- country_map_sf %>% 
  filter(ID %in% c("China", "Japan", "Vietnam", "Philippines", "Palau", "Taiwan", "Indonesia", "Mongolia", 
                   "Laos", "Myanmar", "India", "Thailand", "Burma", "Cambodia", "Malaysia","Bangladesh", 
                   "North Korea", "South Korea", "Singapore", "Brunei", "Papua New Guinea"))

native_coordinates_mapped <- ggplot(data = native_basemap_sf) +
  geom_sf(fill = 'gray97', color = 'black') +
  geom_sf(data = native_sf, aes(color = Collection.Year), size = 3) +
  scale_color_viridis_c(limits = global_year_limits) +
  labs(color = "Collection Year", x = 'Longitude (\u00B0)', y = 'Latitude (\u00B0)') +
  theme_classic() +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(95, 143), ylim = c(-11, 40), crs = 4326, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(1, "in"),
                         style = north_arrow_fancy_orienteering)

# Save Native Map
ggsave("native_samples_coords.pdf", plot = native_coordinates_mapped, width = 6.5, height = 7)
