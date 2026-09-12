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


# Analysis directories -----
interim_results <- file.path(
  project_dir,
  "r_analysis",
  "19_code_interim_results",
  "19.3_Res_CODE3"
)

# Final shapefiles
shapefile_dir <- file.path(
  project_dir,
  "r_analysis",
  "16_Qgis_project",
  "16.1_Final_shapefiles"
)

#Results
results_dir <- file.path(
  project_dir,
  "r_analysis",
  "19_code_interim_results",
  "19.4_Res_CODE4"
)

figure_dir <- file.path(
  project_dir,
  "r_analysis",
  "18_Images"
)


# --- capas extra ---

apn <- vect(file.path(
  shapefile_dir,
  "16.1.4_APN_roi",
  "apn_selected_PatagoniaNorte.shp"
))

lagos <- vect(file.path(
  shapefile_dir,
  "16.1.7_Lakes_roi",
  "2_lagos_study_area.shp"
))

urbano <- vect(file.path(
  shapefile_dir,
  "16.1.5_Urban_roi",
  "4_poblados_study_area.shp"
))

# Carga de métricas desde archivo -----
grid_metrics_vect <- vect( file.path(interim_results,
  "grid_metrics_vect_2023.gpkg")
)

# Número de parches "np" bajo
# Área invadida "ca" bajo
# Largest patch index - "lpi" bajo
# Distancia media entre patches "enn" alta

plot(grid_metrics_vect, "lpi")
names(apn)
unique(apn$nombre)
nh <- apn[apn$nombre == "LanÃ­n", ]
grid_nh <- crop(grid_metrics_vect, nh)
lagos_nh <- crop(lagos, nh)

plot(
  grid_nh,
  "lpi"
)
cols <- c("#440154", "#21908C", "#FDE725")
q <- quantile(
  grid_nh$lpi,
  probs = c(0, 0.33, 0.66, 1),
  na.rm = TRUE
)

plot(
  grid_nh,
  "lpi",
  border = NA,
  col = NA,
  breaks = q,
  legend = FALSE
)

plot(
  lagos_nh,
  col = "lightblue",
  border = NA,
  add = TRUE
)

plot(
  grid_nh,
  "lpi",
  add = TRUE,
  border = NA,
  col = hcl.colors(5, "YlOrRd")[1:3],
  breaks = 3
)

plot(
  nh,
  add = TRUE,
  col = NA,
  lwd = 2
)


nh_buffer <- buffer(nh, width = 3000)
urbano_nh <- crop(urbano, nh_buffer)

analysis_area <- crop(grid_metrics_vect, urbano_nh)
analysis_area <- union(analysis_area, grid_nh)

plot(
  analysis_area,
  "enn"
)


df <- as.data.frame(grid_metrics_vect)

q_np  <- quantile(df$np,  c(0.33, 0.66), na.rm = TRUE)
q_ca  <- quantile(df$ca,  c(0.33, 0.66), na.rm = TRUE)
q_lpi <- quantile(df$lpi, c(0.33, 0.66), na.rm = TRUE)
q_enn <- quantile(df$enn, c(0.33, 0.66), na.rm = TRUE)

cor(
  df[, c("ca", "lpi", "np", "enn")],
  use = "pairwise.complete.obs"
)

df_inv <- df %>%
  filter(!is.nan(ca))

df <- df %>%
  mutate(
    
    erad_score =
      (ca < q_ca[1]) +
      (lpi < q_lpi[1]) +
      (np < q_np[1]),
    
    cont_score =
      (ca > q_ca[2]) +
      (lpi > q_lpi[2]) +
      (np > q_np[2]),
    
    manejo = case_when(
      is.nan(ca) ~ "No invasion",
      erad_score >= 2 ~ "Eradication",
      cont_score >= 2 ~ "Containment",
      TRUE ~ "Restoration"
    )
  )


grid_metrics_vect$manejo <- df$manejo
prop.table(table(df$manejo))

terra::writeVector(
  grid_metrics_vect,
  file.path(results_dir,"grid_index_vect_2023.gpkg"),
  overwrite = TRUE
)
