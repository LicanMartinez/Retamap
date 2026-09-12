# Limpieza del enviroment ------ 
ls () 
rm (list = ls()) 
ls ()

# paquetes ----

library(terra)
library(landscapemetrics)
library(dplyr)

# Project directory
project_dir <- here::here()

# Grillas disponibles -------------------------------------------------------------------
grids_dir <- file.path(
  project_dir,
  "r_analysis",
  "15_Grids_with_predictors",
  "15.1_Empty_grids")

raster_dir <- file.path(
  project_dir,
  "r_analysis",
  "14_Final_rasters_RF2_clean_invasion_per_year"
)

#Results
results_dir <- file.path(
  project_dir,
  "r_analysis",
  "19_code_interim_results",
  "19.3_Res_CODE3"
)

figure_dir <- file.path(
  project_dir,
  "r_analysis",
  "18_Images"
)

#-------------------------------------------------------------------------------
# Grillas vacías
grid7000 <- read.csv(file.path(grids_dir,"grid_7000m_clipped.csv"))
grid7000 <- as.data.frame(grid7000)

grid3000 <- read.csv(file.path(grids_dir,"grid_3000m_clipped.csv"))
grid3000 <- as.data.frame(grid3000)

grid1000 <-read.csv(file.path(grids_dir,"grid_1000m_clipped.csv"))
grid1000 <- as.data.frame(grid1000)
colnames(grid1000)

# Stack rasters
years <- 2023

raster_suffix <- "_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k.tif"

files <- file.path(
  raster_dir,
  paste0("08_RF2_prediction_", years, raster_suffix)
)

stopifnot(all(file.exists(files)))


files
r_list <- lapply(files, rast)

raster <- rast(r_list)

# convertir WKT a vector espacial
grid_vect <- vect(
  grid1000,
  geom = "WKT",
  crs = crs(raster)
)
grid_vect$plot_id <- 1:nrow(grid_vect)

# Número de patches
np <- sample_lsm(
  landscape = raster,
  y = grid_vect,
  what = "lsm_c_np"
)

# Área invadida

ca <- sample_lsm(
  landscape = raster,
  y = grid_vect,
  what = "lsm_c_ca"
)

# Largest patch index
lpi <- sample_lsm(
  landscape = raster,
  y = grid_vect,
  what = "lsm_c_lpi"
)

# Distancia media entre patches
enn <- sample_lsm(
  landscape = raster,
  y = grid_vect,
  what = "lsm_c_enn_mn"
)

# limpiar tablas
np2 <- np %>%
  filter(class == 1) %>%   # SOLO retama
  select(plot_id, value) %>%
  rename(np = value)

ca2 <- ca %>%
  filter(class == 1) %>%
  select(plot_id, value) %>%
  rename(ca = value)

lpi2 <- lpi %>%
  filter(class == 1) %>%
  select(plot_id, value) %>%
  rename(lpi = value)

enn2 <- enn %>%
  filter(class == 1) %>%
  select(plot_id, value) %>%
  rename(enn = value)

grid_vect_df <- as.data.frame(grid_vect,  
                              geom = "WKT")
grid_vect_df$plot_id <- 1:nrow(grid_vect_df)

grid_metrics <- grid_vect_df %>%
  left_join(np2, by = "plot_id") %>%
  left_join(ca2, by = "plot_id") %>%
  left_join(lpi2, by = "plot_id") %>%
  left_join(enn2, by = "plot_id")

grid_vect$np <- np2$np
grid_vect$ca <- ca2$ca
grid_vect$lpi <- lpi2$lpi
grid_vect$enn <- enn2$enn

grid_metrics_vect <- vect(
  grid_metrics,
  geom = "geometry",
  crs = crs(grid_vect)
)

# Fin cálculo de métricas y guardado ----
writeVector(
  grid_metrics_vect,
  file.path(results_dir,"grid_metrics_vect_2023.gpkg"),
  overwrite = TRUE
)
