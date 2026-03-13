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
           region == "oklahoma" |
           region == "north carolina" |
           region == "tennessee" |
           region == "virginia") # no samples in oklahoma, just for aesthetic purposes

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

coordinates <- read.csv("LYJA_invasion_genomics/files/lyja_si_coordinates.csv")
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

stacked.maps <- ggarrange(invaded_coordinates_mapped, native_coordinates_mapped, ncol = 1, common.legend = TRUE, legend = "right")


## Adding points from GBIF to look at the progression of the invasion over time

library(rgbif)

gbif_occs <- occ_data(scientificName = "Lygodium japonicum", limit = 10000, continent = "north_america", year = '1900, 1950')
gbif_occs.df <- gbif_occs$data
gbif_occs_earlyInv <- gbif_occs.df %>% 
  filter(decimalLongitude < 0)

gbif_occs <- occ_data(scientificName = "Lygodium japonicum", limit = 10000, continent = "north_america", year = '1950, 1975')
gbif_occs.df <- gbif_occs$data
gbif_occs_MidInv <- gbif_occs.df %>% 
  filter(decimalLongitude < 0)

gbif_occs <- occ_data(scientificName = "Lygodium japonicum", limit = 10000, continent = "north_america", year = '1975, 2000')
gbif_occs.df <- gbif_occs$data
gbif_occs_lateMidInv <- gbif_occs.df %>% 
  filter(decimalLongitude < 0)

gbif_occs <- occ_data(scientificName = "Lygodium japonicum", limit = 10000, continent = "north_america", year = '2000, 2025')
gbif_occs.df <- gbif_occs$data
gbif_occs_lateInv <- gbif_occs.df %>% 
  filter(decimalLongitude < 0)
  
# early invasion 
early_invasion_map <- invaded_basemap +
  geom_point(gbif_occs_earlyInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  xlim(-110, -75) + ylim(25, 38) + ggtitle("1900-1950")

early_invasion_map

# early mid-invasion 
Mid_invasion_map <- invaded_basemap +
  geom_point(gbif_occs_earlyInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_earlyMidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_MidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  xlim(-110, -75) + ylim(25, 38) + ggtitle("1950-1975")

Mid_invasion_map

# Late mid invasion 
lateMid_invasion_map <- invaded_basemap +
  geom_point(gbif_occs_earlyInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_earlyMidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_MidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_lateMidInv, mapping = aes(x =decimalLongitude, y = decimalLatitude),
             group="Abbreviation", size= 2, color = "gray32") +
  xlim(-110, -75) + ylim(25, 38) + ggtitle("1975-2000")

lateMid_invasion_map

# Late invasion 
late_invasion_map <- invaded_basemap +
  geom_point(gbif_occs_earlyInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_earlyMidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_MidInv, mapping=aes(x=decimalLongitude, y=decimalLatitude), 
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_lateMidInv, mapping = aes(x =decimalLongitude, y = decimalLatitude),
             group="Abbreviation", size= 2, color = "gray32") +
  geom_point(gbif_occs_lateInv, mapping = aes(x =decimalLongitude, y = decimalLatitude),
             group="Abbreviation", size= 2, color = "gray32") +
  xlim(-110, -75) + ylim(25, 38) + ggtitle("2000-2025")

late_invasion_map

grid.arrange(early_invasion_map, Mid_invasion_map, lateMid_invasion_map, late_invasion_map)
