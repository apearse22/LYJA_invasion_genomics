###################################################

# Creating maps of environmental bioclim variables
# Code developed: Abby Pearse, Jessie Pelosi
# Last updated: 02/24/2026

###################################################


### Libraries

library(raster)
library(sp)
library(rnaturalearth)
library(devtools)
library(readr)
library(psych)
library(ggplot2)
library(ggfortify)
library(dplyr)

### Reading in files

lyja_si <- read.csv("data/lyja_si_coordinates.csv")

################################ Cleaning dataframes / generating coordinates ###################################

### Cleaning dataframes

lyja_inv <- lyja_si %>% 
  dplyr::select(SRR, Voucher.Identifier, Approximate.Latitude, Approximate.Longitude, Invaded) %>% 
  filter(Invaded == "Y") # creating a data frame of relevant information for invaded samples

lyja_lat_long <- lyja_inv %>% 
  dplyr::select(Approximate.Latitude, Approximate.Longitude) # creating a data frame of coordinates

lyja_lat_long$Approximate.Latitude <- as.numeric(lyja_lat_long$Approximate.Latitude)
lyja_lat_long$Approximate.Longitude <- as.numeric(lyja_lat_long$Approximate.Longitude)
lyja_lat_long <- na.omit(lyja_lat_long)


### Generating coordinates

lyja_sp <- SpatialPoints(lyja_lat_long)

coordinates(lyja_lat_long) <- ~ Approximate.Longitude + Approximate.Latitude 

myCRS1 <- CRS("WGS84")
crs(lyja_lat_long) <- myCRS1


### Creating base maps

countries <- ne_countries(scale = 10, returnclass = "sf")
states <- ne_states(country = "united states of america", returnclass = "sf") # creates state boundaries

reprojected_map <- sf::st_transform(countries, "WGS84")


############################################################# Obtaining bioclim variables ####################################


### Cleaning data

clim_list <- list.files("data/GEA_files/bioclim/wc2.1_30s_bio/", pattern=".tif", full.names = TRUE)
clim <- raster::stack(clim_list)

conditions_occ <- extract(clim, lyja_lat_long)
conditions_occ_df <- as.data.frame(conditions_occ)

colnames(conditions_occ_df) <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19",
                                 "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9")

lyja_inv <- filter(lyja_inv, SRR != "SRR29127793") # this sample does not have any coordinates associated with it
lyja_bioclim <- cbind(lyja_inv, conditions_occ_df) # creating data frame of sample, sample info, and associated bioclim variables


### Creating a pairs panel between bioclim variables (visualizes corrleation)

pairs.panels(lyja_bioclim[, 6:24], scale=TRUE, gap=0, method = "pearson")


############################################### Accessing bioclim variables #########################

# limiting bioclim variables to our georgraphic area of interest 
e_coords <- extent(-100, -75, 25, 37)
inv.clim <- crop(clim, e_coords)


bio1 <-  as.data.frame(inv.clim$wc2.1_30s_bio_1, xy = TRUE)
bio2 <-  as.data.frame(inv.clim$wc2.1_30s_bio_2, xy = TRUE)
bio3 <-  as.data.frame(inv.clim$wc2.1_30s_bio_3, xy = TRUE)
bio4 <-  as.data.frame(inv.clim$wc2.1_30s_bio_4, xy = TRUE)
bio5 <-  as.data.frame(inv.clim$wc2.1_30s_bio_5, xy = TRUE)
bio6 <-  as.data.frame(inv.clim$wc2.1_30s_bio_6, xy = TRUE)
bio7 <-  as.data.frame(inv.clim$wc2.1_30s_bio_7, xy = TRUE)
bio8 <-  as.data.frame(inv.clim$wc2.1_30s_bio_8, xy = TRUE)
bio9 <-  as.data.frame(inv.clim$wc2.1_30s_bio_9, xy = TRUE)
bio10 <-  as.data.frame(inv.clim$wc2.1_30s_bio_10, xy = TRUE)
bio11 <-  as.data.frame(inv.clim$wc2.1_30s_bio_11, xy = TRUE)
bio12 <-  as.data.frame(inv.clim$wc2.1_30s_bio_12, xy = TRUE)
bio13 <- as.data.frame(inv.clim$wc2.1_30s_bio_13, xy = TRUE)
bio14 <- as.data.frame(inv.clim$wc2.1_30s_bio_14, xy = TRUE)
bio15 <- as.data.frame(inv.clim$wc2.1_30s_bio_15, xy = TRUE)
bio16 <- as.data.frame(inv.clim$wc2.1_30s_bio_16, xy = TRUE)
bio17 <- as.data.frame(inv.clim$wc2.1_30s_bio_17, xy = TRUE)
bio18 <- as.data.frame(inv.clim$wc2.1_30s_bio_18, xy = TRUE)
bio19 <- as.data.frame(inv.clim$wc2.1_30s_bio_19, xy = TRUE)

##################################### Plotting bioclim variables in invaded range ########################################

# Here, we are plotting the gradient of environmental bioclim variables against our sample coordinates to determine the amount of variation occurring
# among different sample locations

### Bioclim #1: Annual mean temperature (C)

bio1_plot <- ggplot() +
  geom_raster(data = bio1, aes(x, y, fill = wc2.1_30s_bio_1)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Annual Mean Temperature (C)") +
  labs(fill = "Annual Mean Temperature (C)")
  
### Bioclim #2: Mean Diurnal Range

bio2_plot <- ggplot() +
  geom_raster(data = bio2, aes(x, y, fill = wc2.1_30s_bio_2)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Mean Diurnal Range (C)") +
  labs(fill = "Mean Diurnal Range (C)")
  
### Bioclim #3: Isothermality

bio3_plot <- ggplot() +
  geom_raster(data = bio3, aes(x, y, fill = wc2.1_30s_bio_3)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Isothermality") +
  labs(fill = "Isothermality")

### Bioclim #4: Temperature Seasonality

bio4_plot <- ggplot() +
  geom_raster(data = bio4, aes(x, y, fill = wc2.1_30s_bio_4)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Temperature Seasonality (C)") +
  labs(fill = "Temperature Seasonality(C)")


### Bioclim #5: Max Temperature of Warmest Month

bio5_plot <- ggplot() +
  geom_raster(data = bio5, aes(x, y, fill = wc2.1_30s_bio_5)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Max Temperature of Warmest Month (C)") +
  labs(fill = "Max Temperature of Warmest Month (C)")


### Bioclim #6: Min Temperature of Coldest Month

bio6_plot <- ggplot() +
  geom_raster(data = bio6, aes(x, y, fill = wc2.1_30s_bio_6)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Min Temperaure of Coldest Month (C)") +
  labs(fill = "Min Temperature of Coldest Month (C)")


### Bioclim #7: Temperature Annual Range

bio7_plot <- ggplot() +
  geom_raster(data = bio7, aes(x, y, fill = wc2.1_30s_bio_7)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Temperature Annual Range (C)") +
  labs(fill = "Temperature Annual Range (C)")

### Bioclim #8: Mean Temperature of Wettest Quarter

bio8_plot <- ggplot() +
  geom_raster(data = bio8, aes(x, y, fill = wc2.1_30s_bio_8)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Mean Temperature of Wettest Quarter (C)") +
  labs(fill = "Mean Temperature of Wettest Quarter (C)")


### Bioclim #9: Mean Temperature of Driest Quarter

bio9_plot <- ggplot() +
  geom_raster(data = bio9, aes(x, y, fill = wc2.1_30s_bio_9)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Mean Temperature of Driest Quarter (C)") +
  labs(fill = "Mean Temperature of Driest Quarter (C)")

### Bioclim #10: Mean Temperature of Warmest Quarter

bio10_plot <- ggplot() +
  geom_raster(data = bio10, aes(x, y, fill = wc2.1_30s_bio_10)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Mean Temperature of Warmest Quarter (C)") +
  labs(fill = "Mean Temperature of Warmest Quarter (C)")

### Bioclim #11: Mean Temperature of Coldest Quarter

bio11_plot <- ggplot() +
  geom_raster(data = bio11, aes(x, y, fill = wc2.1_30s_bio_11)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Mean Temperature of Coldest Quarter (C)") +
  labs(fill = "Mean Temperature of Coldest Quarter (C)")


### Bioclim #12: Annual Precipitation

bio12_plot <- ggplot() +
  geom_raster(data = bio12, aes(x, y, fill = wc2.1_30s_bio_12)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Annual Precipitation (mm)") +
  labs(fill = "Annual Precipitation (mm)")

### Bioclim #13: Precipitation of Wettest Month

bio13_plot <- ggplot() +
  geom_raster(data = bio13, aes(x, y, fill = wc2.1_30s_bio_13)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude)) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Wettest Month (mm)") +
  labs(fill = "Precipitation of Wettest Month (mm)")

### Bioclim #14: Precipitation of Driest Month

bio14_plot <- ggplot() +
  geom_raster(data = bio14, aes(x, y, fill = wc2.1_30s_bio_14)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Driest Month (mm)") +
  labs(fill = "Precipitation of Driest Month (mm)")
  

### Bioclim #15: Precipitation Seasonality

bio15_plot <- ggplot() +
  geom_raster(data = bio15, aes(x, y, fill = wc2.1_30s_bio_15)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation Seasonality (mm)") +
  labs(fill = "Precipitation Seasonality (mm)")

### Bioclim #16: Precipitation of Wettest Quarter

bio16_plot <- ggplot() +
  geom_raster(data = bio16, aes(x, y, fill = wc2.1_30s_bio_16)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Wettest Quarter (mm)") +
  labs(fill = "Precipitation of Wettest Quarter (mm)")


### Bioclim #17: Precipitation of Driest Quarter

bio17_plot <- ggplot() +
  geom_raster(data = bio17, aes(x, y, fill = wc2.1_30s_bio_17)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Driest Quarter (mm)") +
  labs(fill = "Precipitation of Driest Quarter (mm)")


### Bioclim #18: Precipitation of Warmest Quarter

bio18_plot <- ggplot() +
  geom_raster(data = bio18, aes(x, y, fill = wc2.1_30s_bio_18)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Warmest Quarter (mm)") +
  labs(fill = "Precipitation of Warmest Quarter (mm)")

### Bioclim #19: Precipitation of Coldest Quarter

bio19_plot <- ggplot() +
  geom_raster(data = bio19, aes(x, y, fill = wc2.1_30s_bio_19)) +
  scale_fill_distiller(palette = "Spectral", na.value = "transparent") +
  geom_sf(data=reprojected_map, color='black', size=0.05, fill = NA) + 
  geom_sf(data = states, color = 'black', size = 0.05, fill = NA) +
  geom_point(as.data.frame(lyja_lat_long@coords), mapping=aes(x=Approximate.Longitude, y=Approximate.Latitude), size = 1.75) +
  coord_sf(xlim=c(-100, -75), ylim=c(25, 37)) +
  theme_classic() +
  scale_x_continuous(limits = c(-100, -75), expand = c(0,0)) +
  scale_y_continuous(limits = c(25, 37), expand = c(0,0)) +
  xlab('Approximate Longitude') +
  ylab('Approximate Latitude') +
  ggtitle("Precipitation of Coldest Quarter (mm)") +
  labs(fill = "Precipitation of Coldest Quarter (mm)")
  

################################################## PCA of Magnitude ##########################################

# Here, we made a PCA of the magnitude of variance explained by each bioclim variable.

lyja_bioclims_only <- lyja_bioclim %>% 
  dplyr::select(-Voucher.Identifier, -Approximate.Latitude, -Approximate.Longitude, -SRR, -Invaded)

rownames(lyja_bioclims_only) <- lyja_bioclim$SRR

lyja_bioclim_norm <- scale(lyja_bioclims_only)
pca_lyja_bioclim <- stats::princomp(lyja_bioclim_norm, scale = TRUE) 


loadings <- pca_lyja_bioclim$loadings
autoplot(pca_lyja_bioclim, loadings=TRUE, loadings.label=TRUE, 
         loadings.label.colour = "forestgreen", loadings.colour = 'black') + theme_bw()


# We chose to use bioclim variables 12, 14, 2, 15, 13, and 1, as these variables explained a large amount of variance in all directions.

lyja_bioclims_sub <- lyja_bioclims_only %>% 
  dplyr::select(bio12, bio14, bio2, bio15, bio13, bio1)


#write.table(lyja_bioclims_sub, "LYJA_bioclim1-2-13-14-15", sep="\t", quote=FALSE) writing out a file

pairs.panels(lyja_bioclims_sub, scale=TRUE, gap=0, method = "pearson")

