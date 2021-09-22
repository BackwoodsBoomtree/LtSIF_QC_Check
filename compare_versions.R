
library(terra)
library(ncdf4)

new_ver_dir <- "G:/LtSIF/Lite_B10306Ar_r02/2021/04"
old_ver_dir <- "G:/LtSIF/Lite_B10206r_r02/2021/04"
meta_list   <- c("SoundingID", "latitude", "longitude",
                 "MeasurementMode", "FootprintID", "Quality_Flag",
                 "daily_correction_factor", "continuum_radiance_757nm", "continuum_radiance_771nm")
sif_list    <- c("SIF_740nm", "SIF_Daily_740nm", "SIF_Uncertainty_740nm",
                 "SIF_757nm", "SIF_Daily_757nm", "SIF_Uncertainty_757nm",
                 "SIF_771nm", "SIF_Daily_771nm", "SIF_Uncertainty_771nm",
                 "SIF_Relative_757nm", "SIF_Relative_771nm")
geom_list   <- c("SAz", "SZA",
                 "VAz", "VZA",
                 "PA", "RAz")
meteo_list  <- c("specific_humidity", "surface_pressure", "temperature_skin",
                 "temperature_two_meter", "vapor_pressure_deficit", "wind_speed")
cloud_list  <- c("cloud_flag_abp", "co2_ratio", "delta_pressure_abp",
                 "o2_ratio", "surface_albedo_abp")
land_list   <- c("IGBP_index", "sounding_land_fraction")

new_file_list <- list.files(new_ver_dir, full.names = TRUE, pattern = "*.nc4$", recursive = TRUE)
old_file_list <- list.files(old_ver_dir, full.names = TRUE, pattern = "*.nc4$", recursive = TRUE)

# Filter out seq files if they are there
new_file_list <- new_file_list[!grepl("*LtSeq", new_file_list)]
old_file_list <- old_file_list[!grepl("*LtSeq", old_file_list)]

cosd                <- function(degrees) {
  radians <- cos(degrees * pi / 180)
  return(radians)
}
sind                <- function(degrees) {
  radians <- sin(degrees * pi / 180)
  return(radians)
}
compute_phase_angle <- function(sat){
  ## Input is satellite data in data.table format with: sza, vza, saa, vaa
  ## Phase angles are set to negative values if the observational azimuth angle is bigger than the solar azimuth angle
  ## (negative phase angle means the sun is to the right of the satellite)
  necessary_names <- c("saa", "sza", "vaa", "vza")
  if (all(necessary_names %in% names(sat))) {
    vaa <- sat$vaa ## viewing azimuth angle
    vza <- sat$vza ## viewing zenith angle
    sza <- sat$sza ## solar zenith angle
    saa <- sat$saa ## solar azimuth angle
    pa  <- phase <- raa <- rep(NA, length(sza))
    phase[vaa > saa] <- -1.
    phase[vaa < saa] <-  1.
    # Relative azimuth angle
    raa <- abs(saa - vaa)
    for (i in 1:length(raa)) {
      if (raa[i] > 180) {
        raa[i] <- abs(raa[i] - 360)
      }
    }
    pa  <- acos(cosd(sza) * cosd(vza) + sind(vza) * sind(sza) * cosd(raa)) * 180. / pi
    pa  <- pa * phase
    return(pa)
  } else {
    print("!!! Necessary input is missing, function returns NULL !!!")
    return(NULL)
  }
}
build_data          <- function(input_file) {
  env <- new.env()
  df  <- nc_open(input_file) # Open file
  # Metadata
  id    <- ncvar_get(df, "Metadata/SoundingId") # "YYYYMMDDHHMMSS"
  mode  <- ncvar_get(df, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
  f_id  <- ncvar_get(df, "Metadata/FootprintId")
  # Geolocation
  lon_center    <- ncvar_get(df, "Geolocation/longitude")
  lat_center    <- ncvar_get(df, "Geolocation/latitude")
  # Land
  igbp          <- ncvar_get(df, "Science/IGBP_index")
  percent_cover <- ncvar_get(df, "Science/sounding_land_fraction")
  # Cloud
  cloud_flag_abp      <- ncvar_get(df, "Cloud/cloud_flag_abp") # 0 - \"Classified clear\", 1 - \"Classified cloudy\", 2 - \"Not classified\", all other values undefined; not used in SIF processing
  co2_ratio           <- ncvar_get(df, "Cloud/co2_ratio")
  delta_pressure_abp  <- ncvar_get(df, "Cloud/delta_pressure_abp")
  o2_ratio            <- ncvar_get(df, "Cloud/o2_ratio") 
  surface_albedo_abp  <- ncvar_get(df, "Cloud/surface_albedo_abp")
  # SIF and Flag Data
  q_flag     <- ncvar_get(df, "Quality_Flag") # 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
  sif740_D   <- ncvar_get(df, "Daily_SIF_740nm")
  sif757_D   <- ncvar_get(df, "Daily_SIF_757nm")
  sif771_D   <- ncvar_get(df, "Daily_SIF_771nm")
  sif740     <- ncvar_get(df, "SIF_740nm")
  sif757     <- ncvar_get(df, "Science/SIF_757nm")
  sif771     <- ncvar_get(df, "Science/SIF_771nm")
  sif740_U   <- ncvar_get(df, "SIF_Uncertainty_740nm")
  sif757_U   <- ncvar_get(df, "Science/SIF_Uncertainty_757nm")
  sif771_U   <- ncvar_get(df, "Science/SIF_Uncertainty_771nm")
  sif757_R   <- ncvar_get(df, "Science/SIF_Relative_757nm")
  sif771_R   <- ncvar_get(df, "Science/SIF_Relative_771nm")
  daily_corr <- ncvar_get(df, "Science/daily_correction_factor")
  # Radiance
  rad757 <- ncvar_get(df, "Science/continuum_radiance_757nm")
  rad771 <- ncvar_get(df, "Science/continuum_radiance_771nm")
  # Geometry
  sza         <- ncvar_get(df, "SZA")
  saa         <- ncvar_get(df, "SAz")
  vza         <- ncvar_get(df, "VZA")
  vaa         <- ncvar_get(df, "VAz")
  # Meteo
  humidity         <- ncvar_get(df, "Meteo/specific_humidity")
  surface_pressure <- ncvar_get(df, "Meteo/surface_pressure")
  temp_skin        <- ncvar_get(df, "Meteo/temperature_skin")
  temp_2m          <- ncvar_get(df, "Meteo/temperature_two_meter")
  vpd              <- ncvar_get(df, "Meteo/vapor_pressure_deficit")
  wind             <- ncvar_get(df, "Meteo/wind_speed")

  # Phase Angle
  pa_table <- data.frame(sza, vza, saa, vaa)  # build table
  pa       <- compute_phase_angle(pa_table)
  # Relative azimuth angle
  raa <- abs(saa - vaa)
  for (i in 1:length(raa)) {
    if (raa[i] > 180) {
      raa[i] <- abs(raa[i] - 360)
    }
  }
  nc_close(df) # Close nc file
  df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "Quality_Flag" = q_flag,
                   "longitude" = lon_center, "latitude" = lat_center, "FootprintID" = f_id,
                   "SIF_Daily_740nm" = sif740_D, "SIF_Daily_757nm" = sif757_D, "SIF_Daily_771nm" = sif771_D,
                   "SIF_740nm" = sif740, "SIF_Uncertainty_740nm" = sif740_U,
                   "SIF_757nm" = sif757, "SIF_Uncertainty_757nm" = sif757_U, "SIF_Relative_757nm" = sif757_R, "continuum_radiance_757nm" = rad757,
                   "SIF_771nm" = sif771, "SIF_Uncertainty_771nm" = sif771_U, "SIF_Relative_771nm" = sif771_R, "continuum_radiance_771nm" = rad771,
                   "daily_correction_factor" = daily_corr,
                   "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "RAz" = raa, "PA" = pa,
                   "specific_humidity" = humidity, "surface_pressure" = surface_pressure, "temperature_skin" = temp_skin,
                   "temperature_two_meter" = temp_2m, "vapor_pressure_deficit" = vpd, "wind_speed" = wind,
                   "IGBP_index" = igbp, "sounding_land_fraction" = percent_cover,
                   "cloud_flag_abp" = cloud_flag_abp, "co2_ratio" = co2_ratio, "delta_pressure_abp" = delta_pressure_abp,
                   "o2_ratio" = o2_ratio, "surface_albedo_abp" = surface_albedo_abp)
  df <- na.omit(df) # drop rows that contain an NA anywhere
  return(df)
}
match_data          <- function(df1, df2) {
  df <- subset(df1, SoundingID %in% df2$SoundingID)
  return(df)
}
s_plot              <- function(df1, df2, var){
  
  x <- old_file[names(old_file) == var]
  y <- new_file[names(new_file) == var]
  
  x <- x[,1]
  y <- y[,1]
  
  fit <- lm(y ~ x)
  cf <- round(coef(fit), 4)
  eq <- paste0("y = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")
  
  plot(y ~ x, xlab = paste0("Old ", var), ylab = paste0("New ", var),
       cex.lab = 2, cex.axis = 2, cex = 1.5, cex.main = 2, tck = 0.03)
  
  if (var != "SoundingID"){
    mtext(eq, 3, line=-2)
  }
}


new_file <- build_data(new_file_list[1])
old_file <- build_data(old_file_list[1])

new_file <- match_data(new_file, old_file)
old_file <- match_data(old_file, new_file)

new_file <- new_file[order(new_file$SoundingID),]
old_file <- old_file[order(old_file$SoundingID),]

# Plot meta data
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/meta_check.jpg", width = 900, height = 900, quality = 100)
par(mfrow = c(3, 3))

for (i in 1:length(meta_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, meta_list[i])
}
dev.off()

# Plot SIF
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/sif_check.jpg", width = 900, height = 1200, quality = 100)
par(mfrow = c(4, 3))

for (i in 1:length(sif_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, sif_list[i])
}
dev.off()

# Plot geometry
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/geom_check.jpg", width = 1200, height = 800, quality = 100)
par(mfrow = c(2, 3))

for (i in 1:length(geom_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, geom_list[i])
}
dev.off()

# Plot meteo
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/meteo_check.jpg", width = 800, height = 1200, quality = 100)
par(mfrow = c(3, 2))

for (i in 1:length(meteo_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, meteo_list[i])
}
dev.off()

# Plot cloud
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/cloud_check.jpg", width = 1200, height = 800, quality = 100)
par(mfrow = c(2, 3))

for (i in 1:length(cloud_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, cloud_list[i])
}
dev.off()

# Plot land
jpeg("C:/Russell/Git/R/LtSIF_QC_Check/figs/land_check.jpg", width = 800, height = 400, quality = 100)
par(mfrow = c(1, 2))

for (i in 1:length(land_list)) {
  op <- par(mar = c(5.1, 6.1, 4.1, 2.1))
  s_plot(old_file, new_file, land_list[i])
}
dev.off()

