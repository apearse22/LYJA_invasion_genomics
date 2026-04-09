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

#coordinates <- si_table %>% 
  dplyr::select(Approx..Latitude, Approx..Longitude)

#coordinates <- read.csv("files/lyja_si_coordinates.csv")

si_table <- si_table %>% 
  mutate(Approx..Latitude = as.numeric(Approx..Latitude)) %>% 
  mutate(Approx..Longitude = as.numeric(Approx..Longitude))

invaded_coordinates <- si_table %>% 
  filter(Invaded == "Y")

native_coordinates <- si_table %>% 
  filter(Invaded == "N")


### plotting invaded coordinates

colnames(invaded_coordinates) <- c("Ind", "Herbarium", "Voucher.Identifier", "Location", "Lat", "Long", "Invaded", "Collection.Year",
                                   "Temporal_Group", "Specimen.Comments", "total...bases", "Specimen.Link")


myCRS1 <- CRS("WGS84")
crs(invaded_coordinates$Lat, invaded_coordinates$Long) <- myCRS1

invaded_coordinates_mapped <- us_basemap +
  geom_spatial_point(invaded_coordinates, mapping=aes(x=Long, y=Lat, color = Collection.Year), 
             group="Abbreviation", size= 3, crs = "WGS84") +
  scale_color_viridis_c(limits = range(invaded_coordinates$Collection.Year)) +
  labs(color = "Collection Year") +
  annotation_scale(plot_unit = "km", ) +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  coord_sf(xlim = c(-96, -79), ylim = c(25, 36.5), crs = "WGS84") 



ggsave("invaded_sample_coords.pdf", width = 8, height = 5.5)


### plotting native coordinates

native_coordinates_mapped <- native_basemap +
  geom_point(native_coordinates, mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude, color = Collection.Year), 
             group="Abbreviation", size=3) +
  scale_color_viridis_c(limits = range(coordinates$Collection.Year)) +
  labs(color = "Collection Year") # does not plot 7 samples
 
ggsave("native_samples_coords.pdf", width = 6.5, height = 7)
