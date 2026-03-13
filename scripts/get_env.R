###############################################################################

# Get environmental data for invaded range samples to run GEA 
# Code developed by Jessie Pelosi, Abby Pearse
# Last updated: March 2026

###############################################################################


library(dismo)
library(terra)

precip <- rast("env_data/cru_ts4.08.1901.2023.pre.dat.nc")
tmp_max <- rast("env_data/cru_ts4.08.1901.2023.tmx.dat.nc")
tmp_min <- rast("env_data/cru_ts4.08.1901.2023.tmn.dat.nc")

popmap <- read.csv("popmap.csv")

coords_year <- dplyr::select(popmap, Ind, Approx..Latitude, Approx..Longitude, Collection.Year)
coords_year$Approx..Latitude <- as.numeric(coords_year$Approx..Latitude)
coords_year$Approx..Longitude <- as.numeric(coords_year$Approx..Longitude)
coords_year$Collection.Year <- as.numeric(coords_year$Collection.Year)

samples <- na.omit(coords_year)
samples_sv <- vect(samples, geom=c("Approx..Longitude", "Approx..Latitude"), crs = "WGS84")

results_list <- list()

for(i in 1:nrow(samples_sv)) {
  
  # Get year and calculate CRU layer indices
  current_year <- as.numeric(samples_sv$Collection.Year[i])
  start_idx <- ((current_year - 1901) * 12) + 1
  end_idx   <- start_idx + 11
  
  # Extract 12 months for this specific point
  # We use a small buffer (5000m) to avoid NAs on coastal/edge pixels
  p_vals <- extract(precip[[start_idx:end_idx]], samples_sv[i], 
                    method="simple", buffer=5000, fun=mean, na.rm=TRUE)[,-1]
  n_vals <- extract(tmp_min[[start_idx:end_idx]], samples_sv[i], 
                    method="simple", buffer=5000, fun=mean, na.rm=TRUE)[,-1]
  x_vals <- extract(tmp_max[[start_idx:end_idx]], samples_sv[i], 
                    method="simple", buffer=5000, fun=mean, na.rm=TRUE)[,-1]
  
  # Convert to numeric vectors (biovars needs 12-month vectors)
  p_vec <- as.numeric(p_vals)
  n_vec <- as.numeric(n_vals)
  x_vec <- as.numeric(x_vals)
  
  # Calculate Bioclim (handles the NA check)
  if(!any(is.na(c(p_vec, n_vec, x_vec)))) {
    results_list[[i]] <- biovars(p_vec, n_vec, x_vec)
  } else {
    results_list[[i]] <- matrix(NA, nrow=1, ncol=19) # Fill with NAs if data missing
    warning(paste("Sample", i, "at year", current_year, "returned NAs. Check coordinates/land mask."))
  }
}

bioclim_mat <- do.call(rbind, results_list)
colnames(bioclim_mat) <- paste0("BIO", 1:19)
sample_env <- cbind(as.data.frame(samples_sv), bioclim_mat)

