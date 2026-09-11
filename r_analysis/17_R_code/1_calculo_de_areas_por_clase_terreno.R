##########################################
# Author: Sofía Cingolani
# Script: Area quantification
##########################################
# Enviroment cleaning------ 
ls () 
rm (list = ls()) 
ls ()

# packages ---------------------

library(terra)
library(dplyr)
library(ggplot2)
library(egg)

# Graph theme -----------------------------------------------------------
theme_mio <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "plain", color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_blank()
    )
}

theme_set(theme_mio())

# 1 Configuration -------------
# Calculo de pixeles en todo el área de estudio

# Years included in the analysis
years <- 2017:2025

# threshold for determining invasion (presence/absence)
threshold <- 0.1

# main directory

main_dir <- "/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama"

# in and out directories

raster_dir <- file.path(main_dir, "04_final_rasters") # in
figure_dir <- file.path(main_dir, "6_Imagenes")       # out  
# reference coordinates system
crs_ref <- "EPSG:4326"

# 2. Raster preparation ----------------------
files <- paste0(
  "/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/",
  "08_RF2_prediction_", years, "_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k",".tif"
)

files
stack <- rast(files)
names(stack) <- years
stack<- project(stack, crs_ref)

out_file <- "/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/stack_utm_4326.tif"
writeRaster(stack, out_file, overwrite = TRUE)

stack <-rast("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/stack_utm_4326.tif")

# ---- Classify pixels as invaded and count positive pixels by year ----


pixeles_positivos <- global(
  stack > threshold,
  fun = "sum",
  na.rm = TRUE
)
out_file2 <- "/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/positive_stack_utm_4326.tif"

write.csv(pixeles_positivos, out_file2)
pixeles_positivos <-read.csv("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/positive_stack_utm_4326.tif")

  
# Organize annual pixel counts
df_pix <- data.frame(
  year = years,
  n_pixels = pixeles_positivos[,1]
)

df_pix

ggplot(df_pix,aes( x = year,y = n_pixels)) +
  geom_point() +
  scale_x_continuous(breaks=seq(2017, 2025, 1))

# Calculo de área total ----

pixel_size_m <- res(stack)[1]
pixel_area_m2 <- 75.14  # 75.14 m2 deformación del pixel en EPSG:4326 sobre latitud de zona de estudio

df_pix$area_m2 <-df_pix$n_pixels * pixel_area_m2
df_pix$area_ha <- df_pix$area_m2 / 10000


ggplot(df_pix, aes(x = year, y = area_ha)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  scale_x_continuous(breaks=seq(2017, 2025, 1)) +
  labs(
    x = "Año",
    y = "Área (ha)",
    title = "Área clasificada de retama por año"
  ) 

# 4. Land-cover classes --------------------------------------  
# Carga de clases de terreno ------

# Area classes
roi <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/3_study_area_retama.shp")

apn <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/apn_selected_PatagoniaNorte.shp")

poblados <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/4_poblados_study_area.shp")

rutas <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/2_rutas_study_area.shp")

rios <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/2_rios_study_area.shp")

lagos <- vect("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/Area de Estudio/2_lagos_study_area.shp")

roi       <- project(roi, "EPSG:4326")
apn       <- project(apn, "EPSG:4326")
rios      <- project(rios, "EPSG:4326")
rutas     <- project(rutas, "EPSG:4326")
poblados  <- project(poblados, "EPSG:4326")
roi       <- project(roi, "EPSG:4326")

# -----

# buffers
rios_buf  <- buffer(rios, width = 20)
rutas_buf <- buffer(rutas, width = 20)

# APN dentro / fuera
apn_union <- aggregate(apn)
fuera_apn <- erase(roi, apn_union)

# Urbanizaciones dentro / fuera
poblados_union <- aggregate(poblados)
fuera_poblados <- erase(roi, poblados_union)

# funcion para contar pixeles
contar_pixeles <- function(r, zona_vect, nombre_zona, year, threshold, pixel_area_m2){
  
  r_bin <- r > threshold
  
  r_crop <- mask(crop(r_bin, zona_vect), zona_vect)
  
  n_pix <- global(r_crop, "sum", na.rm = TRUE)[1,1]
  
  data.frame(
    anio = year,
    zona = nombre_zona,
    n_pix = n_pix,
    ha_retama = (n_pix * pixel_area_m2) / 10000
  )
}
  
resultados <- list()

zonas <- list(
  "Dentro APN" = apn_union,
  "Fuera APN" = fuera_apn,
  "Buffer rios 20m" = rios_buf,
  "Buffer Rutas 20m" = rutas_buf,
  "Urbanizaciones" = poblados_union,
  "Fuera de Urbanizaciones" = fuera_poblados
)

for(i in 1:nlyr(stack)){
  
  r <- stack[[i]]
  year <- years[i]
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, apn_union, "Dentro APN", year, threshold, pixel_area_m2)
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, fuera_apn, "Fuera APN", year, threshold, pixel_area_m2)
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, rios_buf, "Buffer rios 20m", year, threshold, pixel_area_m2)
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, rutas_buf, "Buffer Rutas 20m", year, threshold, pixel_area_m2)
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, poblados_union, "Urbanizaciones", year, threshold, pixel_area_m2)
  
  resultados[[length(resultados)+1]] <- contar_pixeles(r, fuera_poblados, "Fuera de Urbanizaciones", year, threshold, pixel_area_m2)
}
  
 retama_x_clase <- bind_rows(resultados)

#saveRDS(retama_x_clase, "/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/retama_x_clase.rds")

retama_x_clase <- readRDS("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/04_final_rasters/retama_x_clase.rds")

# ---
ggplot(retama_x_clase, aes(x = anio, y = ha_retama)) +
  geom_line() +
  facet_wrap(vars(zona), nrow = 2) 

retama_x_clase <- retama_x_clase %>%
  mutate(
    zona = factor(zona, levels = c(
      "Dentro APN",
      "Fuera APN",
      "Buffer rios 20m",
      "Buffer Rutas 20m",
      "Urbanizaciones",
      "Fuera de Urbanizaciones"
    ))
  )

retama_x_clase1 <- retama_x_clase %>%
  filter(zona %in% c("Dentro APN",
                     "Fuera APN") )

retama_x_clase2 <- retama_x_clase %>%
  filter(zona %in% c("Urbanizaciones",
                     "Fuera de Urbanizaciones"))

retama_x_clase3 <- retama_x_clase %>%
  filter(zona %in% c("Buffer rios 20m",
                     "Buffer Rutas 20m"))

A <- ggplot(retama_x_clase1, aes(x = anio, y = ha_retama)) +
  geom_line() +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2)) +
  ylim(0, 1200)
A
B <- ggplot(retama_x_clase2, aes(x = anio, y = ha_retama)) +
  geom_line() +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))+
  ylim(0, 900)
B

C <- ggplot(retama_x_clase3, aes(x = anio, y = ha_retama)) +
  geom_line() +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))+
  ylim(0, 250)
C
ggarrange(A, B, C, ncol = 1)


library(gstat) # variograma 
library(spmodel)
library(sp)
library(mgcv)

retama_x_clase$zona <- as.factor(retama_x_clase$zona)

# 5. GAM models ----------------------------------------------
# a: Gaussian model
a <- gam(ha_retama ~ zona + s(anio, by = zona, k = 3),
    retama_x_clase,
    family = gaussian())
summary(a)

newdata <- expand.grid(
  anio = seq(min(retama_x_clase$anio),
             max(retama_x_clase$anio),
             length.out = 100),
  zona = levels(retama_x_clase$zona)
)

pred <- predict(a,
                newdata = newdata,
                se.fit = TRUE)

newdata <- newdata %>%
  mutate(
    fit = pred$fit,
    se = pred$se.fit,
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )


ggplot(newdata,
       aes(x = anio,
           y = fit,
           color = zona,
           fill = zona)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  geom_point(data = retama_x_clase,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  geom_line(data = retama_x_clase,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 2, scales = "free") +
  
  theme_mio()

pred_total <- newdata %>%
  filter(zona %in% c("Urbanizaciones", "Fuera de Urbanizaciones")) %>%
  group_by(anio) %>%
  summarise(
    fit = sum(fit),
    se = sqrt(sum(se^2)),
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )


retama_total <- aggregate(ha_retama ~ anio,
                          retama_x_clase2,
                          sum)
ggplot() +
  
  geom_ribbon(data = pred_total,
              aes(x = anio,
                  ymin = lower,
                  ymax = upper),
              alpha = 0.2) +
  
  geom_line(data = pred_total,
            aes(x = anio,
                y = fit),
            linewidth = 1.2) +
  
  geom_point(data = retama_total,
             aes(x = anio,
                 y = ha_retama)) +
  
  geom_line(data = retama_total,
            aes(x = anio,
                y = ha_retama)) +
  
  theme_mio()

newdata_row1 <- newdata[newdata$zona == "Dentro APN" | newdata$zona == "Fuera APN", ]
newdata_row2 <- newdata[newdata$zona == "Urbanizaciones" | newdata$zona == "Fuera de Urbanizaciones", ]
newdata_row3 <- newdata[newdata$zona == "Buffer rios 20m" | newdata$zona == "Buffer Rutas 20m", ]

A <- ggplot(newdata_row1,
            aes(x = anio,
                y = fit,
                color = zona,
                fill = zona)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  geom_point(data = retama_x_clase1,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  geom_line(data = retama_x_clase1,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))+
  ylim(180,1290) +
  theme_mio()
A

B <- ggplot(newdata_row2,
            aes(x = anio,
                y = fit,
                color = zona,
                fill = zona)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  geom_point(data = retama_x_clase2,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  geom_line(data = retama_x_clase2,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))+
  ylim(300, 900)
B

C <- ggplot(newdata_row3,
            aes(x = anio,
                y = fit,
                color = zona,
                fill = zona)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  geom_point(data = retama_x_clase3,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  geom_line(data = retama_x_clase3,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all")  +
  xlab("Año") +
  ylab("Área ocupada por\n C. scoparius (ha)") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))
C
ggarrange(A, B, C, ncol = 1)

# b: Gamma model 
b <- gam(ha_retama ~ zona + s(anio, by = zona, k = 3),
         retama_x_clase,
         family = Gamma(link = log))
summary(b)

newdata <- expand.grid(
  anio = seq(min(retama_x_clase$anio),
             max(retama_x_clase$anio),
             length.out = 100),
  zona = levels(retama_x_clase$zona)
)

pred <- predict(b,
                newdata = newdata,
                se.fit = TRUE)

newdata <- newdata %>%
  mutate(
    fit = exp(pred$fit),
    lower = exp(pred$fit - 1.96 * pred$se.fit),
    upper = exp(pred$fit + 1.96 * pred$se.fit)
  )


ggplot(newdata,
       aes(x = anio,
           y = fit,
           color = zona,
           fill = zona)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  geom_point(data = retama_x_clase,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE) +
  geom_line(data = retama_x_clase,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 2, scales = "free") +
  
  theme_mio()

pred_total <- newdata %>%
  filter(zona %in% c("Urbanizaciones", "Fuera de Urbanizaciones")) %>%
  group_by(anio) %>%
  summarise(
    fit = sum(fit),
    lower_t = sum(lower),
    upper_t = sum(upper)
  )

retama_total <- aggregate(ha_retama ~ anio,
                          retama_x_clase2,
                          sum)
df_pix

# 7. Figures -------------------------------------------------
# Figura 1: Total invaded area ----- 
fig1 <- ggplot() +
  geom_ribbon(data = pred_total,
              aes(x = anio,
                  ymin = lower_t,
                  ymax = upper_t),
              alpha = 0.2, fill = "#CDAA7D") +
  geom_line(data = pred_total,
            aes(x = anio,
                y = fit),
            linewidth = 1.2, color = "#CDAA7D") +
  geom_point(data = df_pix,
             aes(x = year,
                 y = area_ha)) +
  geom_line(data = df_pix,
            aes(x = year,
                y = area_ha)) +
  xlab("Year") +
  ylab(expression(atop("Area occupied by", 
                  italic("C. scoparius") ~ " (ha)"))) +
  scale_x_continuous(breaks = seq(2017, 2025, 2)) +
  theme_mio() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))
fig1

ggsave("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/6_Imagenes/1_retama_temporal_trends_total.png",
       fig1,
       width = 8,
       height = 4,
       dpi = 300,
       bg = "white")


newdata_row1 <- newdata[newdata$zona == "Dentro APN" | newdata$zona == "Fuera APN", ]
newdata_row2 <- newdata[newdata$zona == "Urbanizaciones" | newdata$zona == "Fuera de Urbanizaciones", ]
newdata_row3 <- newdata[newdata$zona == "Buffer rios 20m" | newdata$zona == "Buffer Rutas 20m", ]


newdata_row1$zona <- factor(newdata_row1$zona,
                               levels = c("Dentro APN", "Fuera APN",
                                          "Urbanizaciones", "Fuera de Urbanizaciones",
                                          "Buffer rios 20m", "Buffer Rutas 20m"),
                               labels = c("National Parks", "Outside National Parks",
                                          "Urban centers", "Outside urban centers",
                                          "Rivers 20m buffer", "Roads 20m buffer"))
newdata_row2$zona <- factor(newdata_row2$zona,
                               levels = c("Dentro APN", "Fuera APN",
                                          "Urbanizaciones", "Fuera de Urbanizaciones",
                                          "Buffer rios 20m", "Buffer Rutas 20m"),
                               labels = c("National Parks", "Outside National Parks",
                                          "Urban centers", "Outside urban centers",
                                          "Rivers 20m buffer", "Roads 20m buffer"))
newdata_row3$zona <- factor(newdata_row3$zona,
                               levels = c("Dentro APN", "Fuera APN",
                                          "Urbanizaciones", "Fuera de Urbanizaciones",
                                          "Buffer rios 20m", "Buffer Rutas 20m"),
                               labels = c("National Parks", "Outside National Parks",
                                          "Urban centers", "Outside urban centers",
                                          "Rivers 20m buffer", "Roads 20m buffer"))

retama_x_clase1$zona <- factor(retama_x_clase1$zona,
                          levels = c("Dentro APN", "Fuera APN",
                                     "Urbanizaciones", "Fuera de Urbanizaciones",
                                     "Buffer rios 20m", "Buffer Rutas 20m"),
                          labels = c("National Parks", "Outside National Parks",
                                     "Urban centers", "Outside urban centers",
                                     "Rivers 20m buffer", "Roads 20m buffer"))
retama_x_clase2$zona <- factor(retama_x_clase2$zona,
                          levels = c("Dentro APN", "Fuera APN",
                                     "Urbanizaciones", "Fuera de Urbanizaciones",
                                     "Buffer rios 20m", "Buffer Rutas 20m"),
                          labels = c("National Parks", "Outside National Parks",
                                     "Urban centers", "Outside urban centers",
                                     "Rivers 20m buffer", "Roads 20m buffer"))
retama_x_clase3$zona <- factor(retama_x_clase3$zona,
                          levels = c("Dentro APN", "Fuera APN",
                                     "Urbanizaciones", "Fuera de Urbanizaciones",
                                     "Buffer rios 20m", "Buffer Rutas 20m"),
                          labels = c("National Parks", "Outside National Parks",
                                     "Urban centers", "Outside urban centers",
                                     "Rivers 20m buffer", "Roads 20m buffer"))

A <- ggplot(newdata_row1,
            aes(x = anio,
                y = fit), 
            ) +
  
  geom_line(linewidth = 0.5, color = "#CDAA7D") +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              fill = "#CDAA7D") +
  
  geom_point(data = retama_x_clase1,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE, 
             size = 1) +
  geom_line(data = retama_x_clase1,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab("") +
  scale_x_continuous(breaks = seq(2017, 2025, 2))+
  #ylim(180,1290) +
  theme_mio() +
  theme(
    axis.title.x = element_blank()
  )
A

B <- ggplot(newdata_row2,
            aes(x = anio,
                y = fit)) +
  
  geom_line(linewidth = 0.5, color = "#CDAA7D") +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              fill = "#CDAA7D") +
  
  geom_point(data = retama_x_clase2,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE,size = 1) +
  geom_line(data = retama_x_clase2,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all") +
  xlab("Año") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_x_continuous(breaks = seq(2017, 2025, 2)) +
  #ylim(300, 900) +
  theme_mio () +
  theme(
    axis.title.x = element_blank()
  )
B

C <- ggplot(newdata_row3,
            aes(x = anio,
                y = fit)) +
  
  geom_line(linewidth = 0.5, color = "#CDAA7D") +
  
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              fill = "#CDAA7D") +
  
  geom_point(data = retama_x_clase3,
             aes(y = ha_retama, x = anio),
             inherit.aes = FALSE, size = 1) +
  geom_line(data = retama_x_clase3,
            aes(y = ha_retama, x = anio),
            inherit.aes = FALSE) +
  facet_wrap(vars(zona), nrow = 1, axes = "all")  +
  xlab("Year") +
  ylab("") +
  scale_x_continuous(breaks = seq(2017, 2025, 2)) +
  theme_mio() 
C

# Figura 2: Invaded area by land-cover class ------
fig2 <- ggarrange(A +
                    theme(axis.text = element_text(size = 15),
                          axis.title = element_text(size = 15),
                          strip.text = element_text(size = 15)), B+
                            theme(axis.text = element_text(size = 15),
                                  axis.title = element_text(size = 15),
                                  strip.text = element_text(size = 15)), C+
                                    theme(axis.text = element_text(size = 15),
                                          axis.title = element_text(size = 15),
                                          strip.text = element_text(size = 15),
                                          plot.margin = margin(r = 10)),
                  ncol = 1)
fig2
ggsave("/home/sofi/Insync/UNIDAD DOCTORADO/2024-RemoteSensing_Retama/6_Imagenes/2_retama_temporal_trends_classes.png",
       fig2,
       width = 8,
       height = 10,
       bg = "white",
       dpi = 300)
