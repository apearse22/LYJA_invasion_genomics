##############################################

# creating map of sample locations 
# code developed by Abby Pearse, Jessie Pelosi
# last updated: 2/2/26

##############################################

### libraries

library(ggplot2)
library(dplyr)
library(ggmap)
library(maps)
library(readr)
library(gridExtra)

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


############################################################ plotting coordinates #####################################################


### getting coordinates from supplmental information table 

coordinates <- read.csv("data/lyja_si_coordinates.csv")
coordinates <- coordinates %>% 
  mutate(Approximate.Latitude = as.numeric(Approximate.Latitude)) %>% 
  mutate(Approximate.Longitude = as.numeric(Approximate.Longitude))

invaded_coordinates <- coordinates %>% 
  filter(Invaded == "Y")
native_coordinates <- coordinates %>% 
  filter(Invaded == "N")


### plotting invaded coordinates

invaded_coordinates_mapped <- invaded_basemap +
  geom_point(invaded_coordinates, mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude, color = Collection.Year), 
             group="Abbreviation", size= 3) +
  scale_color_viridis_c(limits = range(coordinates$Collection.Year)) +
  labs(color = "Collection Year") 
  #geom_segment(aes(x=-78, y=30, xend=-78, yend=32), arrow=arrow(), size=1) +
  #annotate("text", x = -77.5, y = 32, label = "N", color = 'black', size = 5)
  # does not plot a 1940 Florida Sample


### plotting native coordinates

native_coordinates_mapped <- native_basemap +
  geom_point(native_coordinates, mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude, color = Collection.Year), 
             group="Abbreviation", size=3) +
  scale_color_viridis_c(limits = range(coordinates$Collection.Year)) +
  labs(color = "Collection Year") # does not plot 7 samples
 

### stacking plots 

stacked.maps <- ggarrange(invaded_coordinates_mapped, native_coordinates_mapped, ncol = 2, common.legend = TRUE, legend = "right")


  