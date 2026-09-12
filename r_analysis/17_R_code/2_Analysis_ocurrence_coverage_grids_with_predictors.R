# Limpieza del enviroment ------ 
ls () 
rm (list = ls()) 
ls ()

# paquetes --------

library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(DHARMa)
library(gstat) # variograma 
library(spmodel)
library(sp)
library(sf) 
library(spdep)
library(ggpubr)

# Tema para gráficos -----------------------------------------------------------
theme_mio <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.border = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      strip.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "plain", color = "black"),
      legend.key = element_blank(),
    )
}

theme_set(theme_mio())

# Project directory
project_dir <- here::here()

# Grillas disponibles -------------------------------------------------------------------
grids_dir <- file.path(
  project_dir,
  "r_analysis",
  "15_Grids_with_predictors")

grid3000 <- read.csv(file.path(grids_dir,
                               "grid_3000_predictors.csv"))
grid3000 <- as.data.frame(grid3000)

grid1000 <- read.csv(file.path(grids_dir,
                               "grid_1000_predictors.csv"))
grid1000 <- as.data.frame(grid1000)

grid7000 <- read.csv(file.path(grids_dir,
                               "grid_7000_predictors.csv"))
grid7000 <- as.data.frame(grid7000)

#Results
results_dir <- file.path(
  project_dir,
  "r_analysis",
  "19_code_interim_results",
  "19.2_Res_CODE2"
)

figure_dir <- file.path(
  project_dir,
  "r_analysis",
  "18_Images"
)
# Selección de grilla y ordenamiento de datos ---------------------------------
# formato long
grid_long <- grid1000 %>%
  pivot_longer(
    cols = matches("(retama|urbano|no_urbano)_perc_"),
    names_to = c("variable", "year"),
    names_pattern = "(.*)_perc_(\\d+)",
    values_to = "value"
  ) %>%
  mutate(year = as.numeric(year)) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  )
colnames(grid_long)

grid_long$retama_bin <- as.integer(grid_long$retama > 0)

grid_long <- grid_long %>%
  mutate(
    retama_area = (retama / 100) * available_area_m2,
    urbano_area = (urbano / 100) * available_area_m2,
    no_urbano_area = (no_urbano / 100) * available_area_m2
  )


# estandarizar variables ---

vars_cont <- c("urbano", "apn_perc", "precip_mean", 
               "coast_len","rivers_len", "roads_len", 
               "temp_mean", "elev_mean", "fires_perc")

# guardar medias y devío estandar ---- 
scale_params <- lapply(grid_long[vars_cont], function(x) {
  list(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
})

grid_scaled <- grid_long %>%
  mutate(across(all_of(vars_cont), ~ as.numeric(scale(.))))
summary(grid_scaled)

hist(grid_scaled$retama)
mean(grid_scaled$retama == 0)

# variograma 

# vg <- variogram(retama ~ 1, 
# locations = ~ lon + lat, 
# data = grid_scaled)
# plot(vg)

# variograma hay  me dice que hay correlacion positiva
# valores positivos o altos de retama se encuentran cerca de valores similars
# valores negativos o bajos se encuentran cerca de valores similares

pts <- st_as_sf(grid_scaled,
                coords = c("lon", "lat"),
                crs = 4326)

pts_utm <- st_transform(pts, 32719) # ejemplo UTM 19S

coords <- st_coordinates(pts_utm)

grid_scaled$x <- coords[,1]
grid_scaled$y <- coords[,2]

#mismo variograma de los datos, pasa lo mismo tanto para aus/pres y para la cobertura
vg <- variogram(retama_bin ~ 1,
                locations = ~ x + y,
                data = grid_scaled)
plot(vg)

vg <- variogram(retama ~ 1,
                locations = ~ x + y,
                data = grid_scaled)
plot(vg)

# Subset de coberturas positivas ---
grid_pos <- grid_scaled %>%
  group_by(id) %>%
  mutate(max_retama = max(retama, na.rm = TRUE)) %>%
  filter(max_retama > 0) %>%
  ungroup() %>%
  select(-max_retama)

vg <- variogram(retama ~ 1,
                locations = ~ x + y,
                data = grid_pos)
plot(vg)

# Delta cambio porcentaje de cobertura urbano retama ----

grid_changes <- grid_scaled %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    d_retama = retama - lag(retama),
    d_urbano = urbano - lag(urbano),
    d_retama_area = retama_area - lag(retama_area),
    d_urbano_area = urbano_area - lag(urbano_area),
  ) %>%
  ungroup()

grid_changes_pos <- grid_pos %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    d_retama = retama - lag(retama),
    d_urbano = urbano - lag(urbano)
  ) %>%
  ungroup()

# chequeo de correlacion

coords <- cbind(grid_scaled$x, grid_scaled$y)

# vecinos por distancia o contigüidad
nb <- dnearneigh(coords, 0, 2000)  # ej: 2 km
lw <- nb2listw(nb, style = "W")

moran.plot(grid_scaled$retama, lw)

# Autocorrlacion ----

# Coordenadas de las celdas
coords <- cbind(grid_scaled$x, grid_scaled$y)

# Definir vecinos: celdas dentro de 1 km
nb <- dnearneigh(coords, 0, 1000)  # dist min = 0, dist max = 1000m
lw <- nb2listw(nb, style = "W")    # pesos row-standardized

# Moran's I para retama
moran.test(grid_scaled$retama, lw)

moran.plot(grid_scaled$retama, lw,
           xlab = "% Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

grid_scaled_clean <- grid_scaled[complete.cases(
  grid_scaled[, c("urbano", "apn_perc", "fires_perc", 
                  "coast_len", "roads_len", "rivers_len", 
                  "precip_mean", "temp_mean", "elev_mean")]
), ]

# Este modelo espacial es muy pesado por eso vamos a hacer subset de controles


#  spmod0 <- spglm(retama_bin ~ urbano +  
#                 apn_perc + fires_perc + 
#                 coast_len + roads_len + rivers_len +
#                 precip_mean + temp_mean ,
#               grid_scaled_clean,
#               family = "binomial",
#               link = "logit",
#               spcov_type = "gaussian",
#               xcoord = "x",
#               ycoord = "y")


# 1: Celdas con presencia alguna vez
cells_presencia <- grid_scaled_clean %>%
  group_by(id) %>%
  summarize(max_retama = max(retama_bin)) %>%
  filter(max_retama == 1) %>%
  pull(id)

# 2: Celdas que nunca tuvieron presencia
cells_ausencia <- grid_scaled_clean %>%
  group_by(id) %>%
  summarize(max_retama = max(retama_bin)) %>%
  filter(max_retama == 0) %>%
  pull(id)

# 3: Muestreo aleatorio de celdas control # 10 veces más que celdas presencia
set.seed(123)  # reproducible
cells_ausencia_sample <- sample(cells_ausencia, length(cells_presencia)*10)

# 4: Subset del dataset completo. Uso un solo año, porque se repiten muchas celdas iguales por año
grid_subset <- grid_scaled_clean %>%
  filter(id %in% c(cells_presencia, cells_ausencia_sample),
         year == 2024)
table(grid_subset$retama_bin)

# Modelos de Ocupación ---------------------------------------------------------

spmod1_e <- spglm(retama_bin ~ urbano +  
                  apn_perc + fires_perc + 
                  coast_len + roads_len + rivers_len +
                  precip_mean + temp_mean + 
                  elev_mean + I(elev_mean^2),
                grid_subset,
                family = "binomial",
                link = "logit",
                spcov_type = "exponential",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

#saveRDS(spmod1_e, 
#        file = file.path(results_dir,
#                         "spmod1_e_retama_bin_grid1000.rds"))
        

spmod1_g <- spglm(retama_bin ~ urbano +  
                  apn_perc + fires_perc + 
                  coast_len + roads_len + rivers_len +
                  precip_mean + temp_mean + 
                  elev_mean + I(elev_mean^2),
                grid_subset,
                family = "binomial",
                link = "logit",
                spcov_type = "gaussian",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

#saveRDS(spmod1_g, 
#        file = file.path(results_dir,
#                         "spmod1_g_retama_bin_grid1000.rds"))


spmod1_n <- spglm(retama_bin ~ urbano +  
                    apn_perc + fires_perc + 
                    coast_len + roads_len + rivers_len +
                    precip_mean + temp_mean + 
                    elev_mean + I(elev_mean^2),
                  grid_subset,
                  family = "binomial",
                  link = "logit",
                  spcov_type = "none",
                  xcoord = "x",
                  ycoord = "y",
                  local = TRUE)
# saveRDS(spmod1_n, 
#        file = file.path(results_dir,
#                         "spmod1_n_retama_bin_grid1000.rds"))

spmod1_s <- spglm(retama_bin ~ urbano +  
                    apn_perc + fires_perc + 
                    coast_len + roads_len + rivers_len +
                    precip_mean + temp_mean + 
                    elev_mean + I(elev_mean^2),
                  grid_subset,
                  family = "binomial",
                  link = "logit",
                  spcov_type = "spherical",
                  xcoord = "x",
                  ycoord = "y",
                  local = TRUE)
#saveRDS(spmod1_s, 
#        file = file.path(results_dir,
#                         "spmod1_s_retama_bin_grid1000.rds"))

# Modelo de presencia
spmod1.2_s <- spglm(retama_bin ~ urbano + I(urbano^2) + 
                    apn_perc + fires_perc + 
                    coast_len + roads_len + rivers_len +
                    precip_mean + temp_mean + 
                    elev_mean + I(elev_mean^2),
                  grid_subset,
                  family = "binomial",
                  link = "logit",
                  spcov_type = "spherical",
                  xcoord = "x",
                  ycoord = "y",
                  local = TRUE)

#saveRDS(spmod1.2_s, 
#        file = file.path(results_dir,
#                         "spmod1.2_s_retama_bin_grid1000.rds"))
        
spmod1.2_e <- spglm(retama_bin ~ urbano + I(urbano^2) + 
                      apn_perc + fires_perc + 
                      coast_len + roads_len + rivers_len +
                      precip_mean + temp_mean + 
                      elev_mean + I(elev_mean^2),
                    grid_subset,
                    family = "binomial",
                    link = "logit",
                    spcov_type = "exponential",
                    xcoord = "x",
                    ycoord = "y",
                    local = TRUE)

#saveRDS(spmod1.2_e, 
#        file = file.path(results_dir,
#                         "spmod1.2_e_retama_bin_grid1000.rds"))

spmod1.2_g <- spglm(retama_bin ~ urbano + I(urbano^2) + 
                      apn_perc + fires_perc + 
                      coast_len + roads_len + rivers_len +
                      precip_mean + temp_mean + 
                      elev_mean + I(elev_mean^2),
                    grid_subset,
                    family = "binomial",
                    link = "logit",
                    spcov_type = "gaussian",
                    xcoord = "x",
                    ycoord = "y",
                    local = TRUE)

#saveRDS(spmod1.2_g, 
#        file = file.path(results_dir,
#                         "spmod1.2_g_retama_bin_grid1000.rds"))

spmod1.2_n <- spglm(retama_bin ~ urbano + I(urbano^2) + 
                      apn_perc + fires_perc + 
                      coast_len + roads_len + rivers_len +
                      precip_mean + temp_mean + 
                      elev_mean + I(elev_mean^2),
                    grid_subset,
                    family = "binomial",
                    link = "logit",
                    spcov_type = "none",
                    xcoord = "x",
                    ycoord = "y",
                    local = TRUE)

# saveRDS(spmod1.2_n, 
#        file = file.path(results_dir,
#                         "spmod1.2_n_retama_bin_grid1000.rds"))

# Lectura de modelos -----
spmod1_g <- readRDS(file.path(results_dir,
                              "spmod1_g_retama_bin_grid1000.rds"))
spmod1_s <- readRDS(file.path(results_dir,
                              "spmod1_s_retama_bin_grid1000.rds"))
spmod1_e <- readRDS(file.path(results_dir,
                              "spmod1_e_retama_bin_grid1000.rds"))
spmod1_n <- readRDS(file.path(results_dir,
                              "spmod1_n_retama_bin_grid1000.rds"))
spmod1.2_s <- readRDS(file.path(results_dir,
                                "spmod1.2_s_retama_bin_grid1000.rds"))
spmod1.2_e <- readRDS(file.path(results_dir,
                                "spmod1.2_e_retama_bin_grid1000.rds"))
spmod1.2_n <- readRDS(file.path(results_dir,
                                "spmod1.2_n_retama_bin_grid1000.rds"))
spmod1.2_g <- readRDS(file.path(results_dir,
                                "spmod1.2_g_retama_bin_grid1000.rds"))

AIC(spmod1_g,spmod1_s,spmod1_e,spmod1_n,
    spmod1.2_s, spmod1.2_g, spmod1.2_e, spmod1.2_n)

summary(spmod1.2_s)
t_spmod1.2_s <- tidy(spmod1.2_s)
g_spmod1.2_s <- glance(spmod1.2_s)
a_spmod1.2_s <- augment(spmod1.2_s)

plot(a_spmod1.2_s$.fitted, a_spmod1.2_s$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod1.2_s, type = "response")
resid_dev <- residuals(spmod1.2_s, type = "deviance")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod1.2_e$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset[, c("x", "y")])
nb <- dnearneigh(coords, 0, 2000)  # create neighbors
lw <- nb2listw(nb, style="W", zero.policy = T)

moran.test(a_spmod1.2_s$.resid, lw)

moran.plot(a_spmod1.2_s$.resid, lw,
           xlab = "Presencia Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de Presencia Retama")

coords <- grid_subset[, c("x", "y")]
df <- data.frame(x=coords$x, y=coords$y, resid=resid_dev)
coordinates(df) <- ~x+y

# Variograma experimental
vg <- variogram(resid ~ 1, data=df) 
plot(vg)

summary(spmod1.2_s)

AIC(spmod1, spmod1_g, spmod1_s, spmod1_n, spmod1.2_s)

# Predicciones modelo de presencia de retama -----

newdata <- expand.grid(
  urbano = mean(grid_subset$urbano),
  apn_perc = mean(grid_subset$apn_perc),
  fires_perc =  mean(grid_subset$fires_perc),
  coast_len = seq(min(grid_subset$coast_len), 
                  max(grid_subset$coast_len), 
                  length.out = 100),
  roads_len =mean(grid_subset$roads_len),
  rivers_len  =mean(grid_subset$rivers_len),
  precip_mean =mean(grid_subset$precip_mean),
  temp_mean = mean(grid_subset$temp_mean),
  elev_mean = mean(grid_subset$elev_mean),
  x = mean(grid_subset$x),
  y = mean(grid_subset$y)
)

pred_new <- predict(spmod1.2_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- plogis(pred_new$fit)
newdata$lower <- plogis(pred_new$lwr)
newdata$upper <- plogis(pred_new$upr)

newdata$coast_len <- newdata$coast_len * scale_params$coast_len$sd +
  scale_params$coast_len$mean

fig_coast_p <-
  ggplot(newdata, aes(x = coast_len/1000, y = fit)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Cumulative coast length (km)") +
  ylab(expression(atop("Probability of", 
                       italic("C. scoparius") ~ " occurrence"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.5)) +
  scale_x_continuous(expand = c(0,0)) 

# 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = seq(min(grid_subset$roads_len), 
                    max(grid_subset$roads_len), 
                    length.out = 100),
    rivers_len  =mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$roads_len <- newdata$roads_len * scale_params$roads_len$sd + 
    scale_params$roads_len$mean
  
fig_roads_p <-
  ggplot(newdata, aes(x = roads_len/1000, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Cumulative roads length (km)") +
  ylab(expression(atop("Probability of", 
                       italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_x_continuous(expand = c(0,0)) 

# 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = seq(min(grid_subset$rivers_len), 
                      max(grid_subset$rivers_len), 
                      length.out = 100),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$rivers_len <- newdata$rivers_len * scale_params$rivers_len$sd +
    scale_params$rivers_len$mean

fig_rivers_p <- 
  ggplot(newdata, aes(x = rivers_len/1000, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Cumulative rivers length (km)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
    scale_x_continuous(expand = c(0,0)) 
  
#   
  newdata <- expand.grid(
    urbano = seq(min(grid_subset$urbano), 
                 max(grid_subset$urbano), 
                 length.out = 100),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len =mean(grid_subset$roads_len),
    rivers_len  =mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$urbano <- newdata$urbano * scale_params$urbano$sd + 
    scale_params$urbano$mean

fig_urbano_p <-   
  ggplot(newdata, aes(x = urbano, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Urban area (%)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.1)) +
    scale_x_continuous(expand = c(0,0)) 

# 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = mean(grid_subset$rivers_len),
    precip_mean = seq(min(grid_subset$precip_mean), 
                   max(grid_subset$precip_mean), 
                   length.out = 100),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$precip_mean <- newdata$precip_mean * scale_params$precip_mean$sd +
    scale_params$precip_mean$mean

fig_pp_p <-  
  ggplot(newdata, aes(x = precip_mean, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Precipitation mean (mm)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25)) +
    scale_x_continuous(expand = c(0,0)) 
  

  # 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = seq(min(grid_subset$temp_mean), 
                    max(grid_subset$temp_mean), 
                    length.out = 100),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$temp_mean <- newdata$temp_mean * scale_params$temp_mean$sd +
    scale_params$temp_mean$mean

fig_temp_p <-   
  ggplot(newdata, aes(x = temp_mean, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Mean temperature (C°)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25)) +
    scale_x_continuous(expand = c(0,0)) 
  
  # 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = seq(min(grid_subset$elev_mean), 
                    max(grid_subset$elev_mean), 
                    length.out = 100),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$elev_mean <- newdata$elev_mean * scale_params$elev_mean$sd +
    scale_params$elev_mean$mean

fig_elev_p <-   
  ggplot(newdata, aes(x = elev_mean, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Mean elevation (m.s.n.m.)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.05)) +
    scale_x_continuous(expand = c(0,0), breaks = c(300,500,1000,1500)) 
fig_elev_p
  # 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = mean(grid_subset$apn_perc),
    fires_perc =  seq(min(grid_subset$fires_perc), 
                      max(grid_subset$fires_perc), 
                      length.out = 100),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$fires_perc <- newdata$fires_perc * scale_params$fires_perc$sd + 
    scale_params$fires_perc$mean
 
fig_fires_p <- 
  ggplot(newdata, aes(x = fires_perc, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("Fire affected areas (%)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25)) +
    scale_x_continuous(expand = c(0,0)) 
  
  # 
  newdata <- expand.grid(
    urbano = mean(grid_subset$urbano),
    apn_perc = seq(min(grid_subset$apn_perc), 
                   max(grid_subset$apn_perc), 
                   length.out = 100),
    fires_perc =  mean(grid_subset$fires_perc),
    coast_len = mean(grid_subset$coast_len),
    roads_len = mean(grid_subset$roads_len),
    rivers_len  = mean(grid_subset$rivers_len),
    precip_mean =mean(grid_subset$precip_mean),
    temp_mean = mean(grid_subset$temp_mean),
    elev_mean = mean(grid_subset$elev_mean),
    x = mean(grid_subset$x),
    y = mean(grid_subset$y)
  )
  
  pred_new <- predict(spmod1.2_s, newdata = newdata, 
                      interval = "confidence") %>% as.data.frame()
  
  newdata$fit <- plogis(pred_new$fit)
  newdata$lower <- plogis(pred_new$lwr)
  newdata$upper <- plogis(pred_new$upr)
  
  newdata$apn_perc <- newdata$apn_perc * scale_params$apn_perc$sd +
    scale_params$apn_perc$mean
  
fig_apn_p <-
   ggplot(newdata, aes(x = apn_perc, y = fit)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lower,
                    ymax = upper),
                alpha = 0.2,
                fill = "#CDAA7D") + 
    xlab("National Park area (%)") +
    ylab(expression(atop("Probability of", 
                         italic("C. scoparius") ~ " occurrence"))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,0.25)) +
    scale_x_continuous(expand = c(0,0)) 

# FIN Predicciones modelo de presencia de retama -----
#Figura OCUPACIÓN en función de todas las predictoras ---
fig6 <- ggpubr::ggarrange(
  
  fig_roads_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 15,
             y = 0.25,
             label = "***",
             size = 5),
  
  fig_rivers_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 0.75,
             label = "***",
             size = 5),
  
  fig_coast_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 0.4,
             label = "***",
             size = 5),
  
  fig_urbano_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 70,
             y = 0.075,
             label = "Linear: ***\nQuadratic: ***",
             size = 5),
  
  fig_apn_p +
    theme(axis.title.y = element_blank()),
  
  fig_fires_p +
    theme(axis.title.y = element_blank()),
  
  fig_pp_p +
    theme(axis.title.y = element_blank()),
  
  fig_temp_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 9,
             y = 0.2,
             label = "*",
             size = 5),
  
  fig_elev_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 1400,
             y = 0.04,
             label = "Linear: ns\nQuadratic: *",
             size = 5),
  
  nrow = 3,
  ncol = 3,
  align = "hv"
)

fig6 <- 
  annotate_figure(
    fig6,
    left = text_grob(expression("Probability of "* 
                                italic("C. scoparius") ~ " occurrence"),
                     rot = 90))
fig6

#Figura OCUPACIÓN en función de predictoras SIGNIFICATIVAS---
fig7 <- ggpubr::ggarrange(
  
  fig_roads_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 15,
             y = 0.25,
             label = "***",
             size = 5),
  
  fig_rivers_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 0.75,
             label = "***",
             size = 5),
  
  fig_coast_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 0.4,
             label = "***",
             size = 5),
  
  fig_urbano_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 70,
             y = 0.075,
             label = "Linear: ***\nQuadratic: ***",
             size = 5),

  fig_temp_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 9,
             y = 0.2,
             label = "*",
             size = 5),
  
  fig_elev_p +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 1400,
             y = 0.04,
             label = "Linear: ns\nQuadratic: *",
             size = 5),
  
  nrow = 2,
  ncol = 3,
  align = "hv"
)
fig7
fig7 <- 
  annotate_figure(
    fig7,
    left = text_grob(expression("Probability of " * 
                                  italic("C. scoparius") ~ " occurrence"),
                     rot = 90))

# Guardado figuras -----
ggsave(filename = file.path(figure_dir,
                 "6_retama_precense_fx_predictoras_all.png"),
       fig6,
       width = 9,
       height = 6,
       dpi = 300,
       bg = "white")

ggsave(filename = file.path(figure_dir,
                 "7_retama_precense_fx_significant_predictoras.png"),
       fig7,
       width = 9,
       height = 4,
       dpi = 300,
       bg = "white")

# Modelos de cobertura -------------------------------

# Subset solo con celdas con presencia alguna vez
grid_subset_positivos <- grid_scaled_clean %>%
  filter(id %in% cells_presencia)

spmod2 <- splm(retama_area ~ urbano +  
                  apn_perc + fires_perc + 
                  coast_len + roads_len + rivers_len +
                  precip_mean + temp_mean + elev_mean,
                grid_subset_positivos,
                spcov_type = "gaussian",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

# saveRDS(spmod2, file = file.path(results_dir,
#                 "spmod2_g_retama_bin_grid1000.rds"))
 
spmod2 <- readRDS(file = file.path(results_dir,
                                   "spmod2_g_retama_bin_grid1000.rds"))

summary(spmod2)
t_spmod2 <- tidy(spmod2)
g_spmod2 <- glance(spmod2)
a_spmod2 <- augment(spmod2)

plot(a_spmod2$.fitted, a_spmod2$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod2, type = "response")
resid_dev <- residuals(spmod2, type = "pearson")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod2$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod2$.resid, lw)

moran.plot(a_spmod2$.resid, lw,
           xlab = "Area de Retama (m2)",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")
# 
grid_subset_positivos <- grid_scaled_clean %>%
  filter(id %in% cells_presencia)

grid_subset_positivos2 <- grid_subset_positivos[grid_subset_positivos$retama_area > 0, ]
grid_subset_positivos3 <- grid_subset_positivos2
table(grid_subset_positivos2$retama == 0)
max(grid_subset_positivos2$retama)

spmod2.1 <- spglm(retama_area ~ urbano +  
                 apn_perc + fires_perc + 
                 coast_len + roads_len + rivers_len +
                 precip_mean + temp_mean + elev_mean,
               grid_subset_positivos2,
               family = "Gamma",
               link = "log",
               spcov_type = "gaussian",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)


summary(spmod2.1)
t_spmod2.1 <- tidy(spmod2.1)
g_spmod2.1 <- glance(spmod2.1)
a_spmod2.1 <- augment(spmod2.1)

plot(a_spmod2.1$.fitted, a_spmod2.1$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod2.1, type = "response")
resid_dev <- residuals(spmod2.1, type = "pearson")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod2.1$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset_positivos2[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W",zero.policy = T)

moran.test(a_spmod2.1$.resid, lw)

moran.plot(a_spmod2.1$.resid, lw,
           xlab = "Area de Retama (m2)",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

spmod2.1_n <- spglm(retama_area ~ urbano +  
                    apn_perc + fires_perc + 
                    coast_len + roads_len + rivers_len +
                    precip_mean + temp_mean + elev_mean,
                  grid_subset_positivos2,
                  family = "Gamma",
                  link = "log",
                  spcov_type = "none",
                  xcoord = "x",
                  ycoord = "y",
                  local = TRUE)

spmod2.1_e <- spglm(retama_area ~ urbano +  
                      apn_perc + fires_perc + 
                      coast_len + roads_len + rivers_len +
                      precip_mean + temp_mean + elev_mean,
                    grid_subset_positivos2,
                    family = "Gamma",
                    link = "log",
                    spcov_type = "exponential",
                    xcoord = "x",
                    ycoord = "y",
                    local = TRUE)

spmod2.1_s <- spglm(retama_area ~ urbano +  
                      apn_perc + fires_perc + 
                      coast_len + roads_len + rivers_len +
                      precip_mean + temp_mean + elev_mean,
                    grid_subset_positivos2,
                    family = "Gamma",
                    link = "log",
                    spcov_type = "spherical",
                    xcoord = "x",
                    ycoord = "y",
                    local = TRUE)

# saveRDS(spmod2.1, file = file.path(results_dir,
#                           "spmod2.1_retama_bin_grid1000.rds"))

spmod2 <- readRDS(file.path(results_dir, "spmod2_g_retama_bin_grid1000.rds"))

spmod2.1 <- readRDS(file.path(results_dir,"spmod2.1_g_retama_bin_grid1000.rds"))
spmod2.1_s <- readRDS(file.path(results_dir,"spmod2.1_s_retama_bin_grid1000.rds"))
spmod2.1_e <- readRDS(file.path(results_dir,"spmod2.1_e_retama_bin_grid1000.rds"))
spmod2.1_n <- readRDS(file.path(results_dir,"spmod2.1_n_retama_bin_grid1000.rds"))

AIC(spmod2.1, spmod2.1_s,spmod2.1_e,spmod2.1_n)

summary(spmod2.1_s)
#variograma datos
vg <- variogram(retama_area ~ 1,
                locations = ~ x + y,
                data = grid_subset_positivos2)
plot(vg)

# variograma modelo
coords <- grid_subset_positivos2[, c("x", "y")]
resid_dev <- residuals(spmod2.1_s, type = "pearson")
df <- data.frame(x=coords$x, y=coords$y, resid=resid_dev)
coordinates(df) <- ~x+y
vg <- variogram(resid ~ 1, data = df)
plot(vg)


# Predicciones modelo de área de retama -----

newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = seq(min(grid_subset_positivos2$coast_len), 
                  max(grid_subset_positivos2$coast_len), 
                      length.out = 100),
    roads_len =mean(grid_subset_positivos2$roads_len),
    rivers_len  =mean(grid_subset_positivos2$rivers_len),
    precip_mean =mean(grid_subset_positivos2$precip_mean),
    temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)
newdata$coast_len <- newdata$coast_len * scale_params$coast_len$sd + 
  scale_params$coast_len$mean
 

fig_coast <-
ggplot(newdata, aes(x = coast_len/1000, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Cumulative coast length (km)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.5)) +
  scale_x_continuous(expand = c(0,0)) 

newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = seq(min(grid_subset_positivos2$roads_len), 
                  max(grid_subset_positivos2$roads_len), 
                  length.out = 100),
  rivers_len  =mean(grid_subset_positivos2$rivers_len),
  precip_mean =mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)
newdata$roads_len <- newdata$roads_len * scale_params$roads_len$sd + 
  scale_params$roads_len$mean


fig_roads <-
ggplot(newdata, aes(x = roads_len/1000, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Cumulative roads length (km)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_x_continuous(expand = c(0,0)) 

newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = seq(min(grid_subset_positivos2$rivers_len), 
                    max(grid_subset_positivos2$rivers_len), 
                    length.out = 100),
  precip_mean =mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$rivers_len <- newdata$rivers_len * scale_params$rivers_len$sd +
  scale_params$rivers_len$mean
newdata %>%
  filter(rivers_len %in% c(min(rivers_len), max(rivers_len))) %>%
  select(rivers_len, fit, lower, upper)

fig_rivers <-
ggplot(newdata, aes(x = rivers_len/1000, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Cumulative rivers length (km)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)")))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1.5)) +
  scale_x_continuous(expand = c(0,0)) 

newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = seq(min(grid_subset_positivos2$precip_mean), 
                    max(grid_subset_positivos2$precip_mean), 
                    length.out = 100),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$precip_mean <- newdata$precip_mean * scale_params$precip_mean$sd +
  scale_params$precip_mean$mean

fig_pp <-
ggplot(newdata, aes(x = precip_mean, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Precipitation mean (mm)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) 


newdata <- expand.grid(
  urbano = seq(min(grid_subset_positivos2$urbano), 
               max(grid_subset_positivos2$urbano), 
               length.out = 100),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$urbano <- newdata$urbano * scale_params$urbano$sd + 
  scale_params$urbano$mean

fig_urbano <-
ggplot(newdata, aes(x = urbano, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Urban area (%)") +
  ylab(expression("Area occupied by " *
                  italic("C. scoparius") ~ " (ha)")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) 


newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc = seq(min(grid_subset_positivos2$fires_perc), 
                   max(grid_subset_positivos2$fires_perc), 
                   length.out = 100),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean =  mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

fig_fires <-
ggplot(newdata, aes(x = fires_perc, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Fire affected areas (%)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) #+
#  geom_point(aes(x = fires_perc, y = retama_area/10000),
#             grid_subset_positivos2,
#             alpha = 0.1)


newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = seq(min(grid_subset_positivos2$apn_perc), 
                 max(grid_subset_positivos2$apn_perc), 
                 length.out = 100),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$apn_perc <- newdata$apn_perc * scale_params$apn_perc$sd + 
  scale_params$apn_perc$mean
grid_subset_positivos3$apn_perc <- grid_subset_positivos2$apn_perc * 
  scale_params$apn_perc$sd + 
  scale_params$apn_perc$mean

fig_apn <-
ggplot(newdata, aes(x = apn_perc, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("National Park area (%)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) 


newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = mean(grid_subset_positivos2$precip_mean),
  temp_mean = seq(min(grid_subset_positivos2$temp_mean), 
                  max(grid_subset_positivos2$temp_mean), 
                  length.out = 100),
  elev_mean = mean(grid_subset_positivos2$elev_mean),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$temp_mean <- newdata$temp_mean * scale_params$temp_mean$sd + scale_params$temp_mean$mean

fig_temp <-
ggplot(newdata, aes(x = temp_mean, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Temperature mean (C°)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  ylim(0,1) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) 

newdata <- expand.grid(
  urbano = mean(grid_subset_positivos2$urbano),
  apn_perc = mean(grid_subset_positivos2$apn_perc),
  fires_perc =  mean(grid_subset_positivos2$fires_perc),
  coast_len = mean(grid_subset_positivos2$coast_len),
  roads_len = mean(grid_subset_positivos2$roads_len),
  rivers_len  = mean(grid_subset_positivos2$rivers_len),
  precip_mean = mean(grid_subset_positivos2$precip_mean),
  temp_mean = mean(grid_subset_positivos2$temp_mean),
  elev_mean = seq(min(grid_subset_positivos2$elev_mean), 
                  max(grid_subset_positivos2$elev_mean), 
                  length.out = 100),
  x = mean(grid_subset_positivos2$x),
  y = mean(grid_subset_positivos2$y)
)

pred_new <- predict(spmod2.1_s, newdata = newdata, 
                    interval = "confidence") %>% as.data.frame()

newdata$fit <- exp(pred_new$fit)
newdata$lower <- exp(pred_new$lwr)
newdata$upper <- exp(pred_new$upr)

newdata$elev_mean <- newdata$elev_mean * scale_params$elev_mean$sd + scale_params$elev_mean$mean

fig_elev <-
ggplot(newdata, aes(x = elev_mean, y = fit/10000)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower/10000,
                  ymax = upper/10000),
              alpha = 0.2,
              fill = "#CDAA7D") + 
  xlab("Elevation mean (m.s.n.m.)") +
  ylab(expression(atop("Area occupied by", 
                       italic("C. scoparius") ~ " (ha)"))) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) 

fig4 <- cowplot::plot_grid(
  
  fig_roads +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 15,
             y = 3.5,
             label = "***",
             size = 5),
  
  fig_rivers +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 1,
             label = "***",
             size = 5),
  
  fig_coast +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 3,
             y = 1,
             label = "***",
             size = 5),
  
  fig_urbano +
    theme(axis.title.y = element_blank()),
  
  fig_apn +
    theme(axis.title.y = element_blank()),
  
  fig_fires +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 75,
             y = 0.75,
             label = "•",
             size = 5),
  
  fig_pp +
    theme(axis.title.y = element_blank()) +
    annotate("text",
             x = 120,
             y = 0.75,
             label = "**",
             size = 5),
  
  fig_temp +
    theme(axis.title.y = element_blank()),
  
  fig_elev +
    theme(axis.title.y = element_blank()),
  
  nrow = 3,
  align = "hv"
)

fig4 <- annotate_figure(
  fig4,
  left = text_grob(expression("Area occupied by " *
                                italic("C. scoparius") ~ " (ha)"),
                   rot = 90))

fig4


fig5 <-
  cowplot::plot_grid(fig_roads +
              theme(axis.title.y = element_blank()
              ) +
              annotate("text",
                       x = 15,
                       y = 3.4,
                       label = "***",
                       size = 5), 
            fig_rivers +
              theme(axis.title.y = element_blank()
              ) +
              annotate("text",
                       x = 3,
                       y = 1,
                       label = "***",
                       size = 5), 
            fig_coast +
              theme(axis.title.y = element_blank()
              ) +
              annotate("text",
                       x = 3,
                       y = 1,
                       label = "***",
                       size = 5),
            fig_fires +
              theme(axis.title.y = element_blank()
              ) +
              annotate("text",
                       x = 75,
                       y = 0.75,
                       label = "•",
                       size = 5),
            fig_pp +
              theme(axis.title.y = element_blank()
              ) +
              annotate("text",
                       x = 120,
                       y = 0.75,
                       label = "**",
                       size = 5), 
            nrow = 2)
fig5

fig5 <- 
annotate_figure(
  fig5,
  left = text_grob(expression("Area occupied by " *
                                italic("C. scoparius") ~ " (ha)"),
                   rot = 90))

ggsave( filename = file.path(
        figure_dir,
        "4_Area_retama_fx_predictoras_all.png"),
       fig4,
       width = 9,
       height = 6,
       dpi = 300,
       bg = "white")

ggsave(filename = file.path(
       figure_dir,
       "5_Area_retama_fx_significant_predictoras.png"),
       fig5,
       width = 9,
       height = 4,
       dpi = 300,
       bg = "white")

#saveRDS(spmod2.1_n, file = file.path(results_dir,"spmod2.1_n_retama_bin_grid1000.rds")

# GLM con modelado de correlacion espacial para el área invadida

grid_subset_positivos$retama_beta <- grid_subset_positivos$retama / 100
grid_subset_positivos$retama_beta[grid_subset_positivos$retama_beta == 0] <- grid_subset_positivos$retama_beta[grid_subset_positivos$retama_beta == 0] + 0.00001

# este lo hice sumando 0.00001 a los ceros
spmod3 <- spglm(retama_beta ~ urbano +  
                 apn_perc + fires_perc + 
                 coast_len + roads_len + rivers_len +
                 precip_mean + temp_mean,
                grid_subset_positivos,
               family = "beta",
               link = "logit",
               spcov_type = "gaussian",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)

# saveRDS(spmod3, file = file.path(results_dir,"spmod3_retama_bin_grid1000.rds")
spmod3 <- readRDS(file.path(results_dir,"spmod3_retama_bin_grid1000.rds"))

summary(spmod3)
t_spmod3 <- tidy(spmod3)
g_spmod3 <- glance(spmod3)
a_spmod3 <- augment(spmod3)

plot(a_spmod3$.fitted, a_spmod3$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod3, type = "response")
resid_dev <- residuals(spmod3, type = "deviance")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod3$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset_positivos[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod3$.resid, lw)

moran.plot(a_spmod3$.resid, lw,
           xlab = "% Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

# este lo hice sumando 0.00001 a todos, tanto con cobertura como a los ceros -
grid_subset_positivos$retama2 <- grid_subset_positivos$retama / 100
grid_subset_positivos$retama2 <- grid_subset_positivos$retama2 + 0.00001

spmod4 <- spglm(retama2 ~ urbano +  
                  apn_perc + fires_perc + 
                  coast_len + roads_len + rivers_len +
                  precip_mean + temp_mean,
                grid_subset_positivos,
                family = "beta",
                link = "logit",
                spcov_type = "gaussian",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

# saveRDS(spmod4, file = file.path(results_dir,"spmod4_retama_bin_grid1000.rds")
spmod4 <- readRDS(file.path(results_dir,"spmod4_retama_bin_grid1000.rds"))

summary(spmod4)
t_spmod4 <- tidy(spmod4)
g_spmod4 <- glance(spmod4)
a_spmod4 <- augment(spmod4)

plot(a_spmod4$.fitted, a_spmod4$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod4, type = "response")
resid_dev <- residuals(spmod4, type = "deviance")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod4$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset_positivos[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod4$.resid, lw)

moran.plot(a_spmod4$.resid, lw,
           xlab = "% Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")
# ---
grid_subset_cob <- grid_subset_positivos[grid_subset_positivos$retama != 0,]
grid_subset_cob$retama <- grid_subset_cob$retama / 100

spmod5 <- spglm(retama ~ urbano +  
                  apn_perc + fires_perc + 
                  coast_len + roads_len + rivers_len +
                  precip_mean + temp_mean,
                grid_subset_cob,
                family = "beta",
                link = "logit",
                spcov_type = "gaussian",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

# saveRDS(spmod5, file = file.path(results_dir,"spmod5_retama_bin_grid1000.rds")
spmod5 <- readRDS(file.path(results_dir,"spmod5_retama_bin_grid1000.rds"))

summary(spmod5)

t_spmod5 <- tidy(spmod5)
g_spmod5 <- glance(spmod5)
a_spmod5 <- augment(spmod5)

plot(a_spmod5$.fitted, a_spmod5$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod5, type = "response")
resid_dev <- residuals(spmod5, type = "deviance")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod5$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_subset_cob[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W", zero.policy = TRUE)

moran.test(a_spmod5$.resid, lw)

moran.plot(a_spmod5$.resid, lw,
           xlab = "% Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

# Modelo para cambio 
grid_subset <- grid_scaled_clean %>%
  filter(id %in% c(cells_presencia, cells_ausencia_sample))

grid_changes_subset <- grid_subset %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(
    d_retama = retama - lag(retama),
    d_urbano = urbano - lag(urbano),
    d_retama_area = retama_area - lag(retama_area),
    d_urbano_area = urbano_area - lag(urbano_area),
  ) %>%
  ungroup()

spmod6 <- splm(d_retama ~ d_urbano,
                grid_changes_subset,
                spcov_type = "gaussian",
                xcoord = "x",
                ycoord = "y",
                local = TRUE)

# saveRDS(spmod6, file = file.path(results_dir,"spmod6_retama_bin_grid1000.rds")
spmod6 <- readRDS(file.path(results_dir,"spmod6_retama_bin_grid1000.rds"))



summary(spmod6)
t_spmod6 <- tidy(spmod6)
g_spmod6 <- glance(spmod6)
a_spmod6 <- augment(spmod6)

plot(a_spmod6$.fitted, a_spmod6$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod6, type = "response")
resid_dev <- residuals(spmod6, type = "response")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod6$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_changes_subset[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod6$.resid, lw)

moran.plot(a_spmod6$.resid, lw,
           xlab = "% Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

grid_total_changes <- grid_subset %>%
  group_by(id) %>%
  mutate(d_retama_total = retama[year == 2024] - retama[year == 2017],
         d_urbano_total = urbano[year == 2024] - urbano[year == 2017]) %>%
  ungroup()

spmod7 <- splm(d_retama_total ~ d_urbano_total,
               grid_total_changes,
               spcov_type = "exponential",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)

# saveRDS(spmod7, file = file.path(results_dir,"spmod7_retama_bin_grid1000.rds")
spmod7 <- readRDS(file.path(results_dir,"spmod7_retama_bin_grid1000.rds"))


summary(spmod7)
t_spmod7 <- tidy(spmod7)
g_spmod7 <- glance(spmod7)
a_spmod7 <- augment(spmod7)

plot(a_spmod7$.fitted, a_spmod7$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod7, type = "response")
resid_dev <- residuals(spmod7, type = "response")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod7$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_total_changes[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod7$.resid, lw)

moran.plot(a_spmod7$.resid, lw,
           xlab = "Delta Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")


# delta Retama principio fin 3 años promediado
grid_total_changes2 <- grid_subset %>%
  group_by(id) %>%
  mutate(d_retama_total_3y = (retama[year == 2024] + retama[year == 2023] + retama[year == 2022] )/ 3 - 
                          (retama[year == 2017] + retama[year == 2018] + retama[year == 2019] )/ 3,
         d_urbano_total_3y = (urbano[year == 2024] + urbano[year == 2023] + urbano[year == 2022] )/ 3 - 
                             (urbano[year == 2017] + urbano[year == 2018] + urbano[year == 2019] )/ 3 ,
         d_retama_area_3y = (retama_area[year == 2024] + retama_area[year == 2023] + retama_area[year == 2022] )/ 3 - 
                            (retama_area[year == 2017] + retama_area[year == 2018] + retama_area[year == 2019] )/ 3,
         d_urbano_area_3y = (urbano_area[year == 2024] + urbano_area[year == 2023] + urbano_area[year == 2022] )/ 3 - 
                            (urbano_area[year == 2017] + urbano_area[year == 2018] + urbano_area[year == 2019] )/ 3) %>%
  ungroup()

spmod8 <- splm(d_retama_total_3y ~ d_urbano_total_3y,
               grid_total_changes2,
               spcov_type = "exponential",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)

# saveRDS(spmod8, file = file.path(results_dir,"spmod8_retama_bin_grid1000.rds")
spmod8 <- readRDS(file.path(results_dir,"spmod8_retama_bin_grid1000.rds"))


summary(spmod8)
t_spmod8 <- tidy(spmod8)
g_spmod8 <- glance(spmod8)
a_spmod8 <- augment(spmod8)


plot(a_spmod8$.fitted, a_spmod8$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod8, type = "response")
resid_dev <- residuals(spmod8, type = "response")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod8$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_total_changes2[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod8$.resid, lw)

moran.plot(a_spmod8$.resid, lw,
           xlab = "Delta Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")

# 
spmod9 <- splm(d_retama_area_3y ~ d_urbano_area_3y,
               grid_total_changes2,
               spcov_type = "exponential",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)

spmod9_n <- splm(d_retama_area_3y ~ d_urbano_area_3y,
               grid_total_changes2,
               spcov_type = "none",
               xcoord = "x",
               ycoord = "y",
               local = TRUE)

spmod9_g <- splm(d_retama_area_3y ~ d_urbano_area_3y,
                 grid_total_changes2,
                 spcov_type = "gaussian",
                 xcoord = "x",
                 ycoord = "y",
                 local = TRUE)

spmod9_s <- splm(d_retama_area_3y ~ d_urbano_area_3y,
                 grid_total_changes2,
                 spcov_type = "spherical",
                 xcoord = "x",
                 ycoord = "y",
                 local = TRUE)


#saveRDS(spmod9_n, file = file.path(results_dir,"spmod9_n_retama_bin_grid1000.rds")

spmod9_e <- readRDS(file.path(results_dir,"spmod9_e_retama_bin_grid1000.rds"))
spmod9_n <- readRDS(file.path(results_dir,"spmod9_n_retama_bin_grid1000.rds"))
spmod9_g <- readRDS(file.path(results_dir,"spmod9_g_retama_bin_grid1000.rds"))
spmod9_s <- readRDS(file.path(results_dir,"spmod9_s_retama_bin_grid1000.rds"))

AIC(spmod9_e, spmod9_n, spmod9_g, spmod9_s)

spmod9 <- spmod9_e
summary(spmod9_e)
t_spmod9 <- tidy(spmod9)
g_spmod9 <- glance(spmod9)
a_spmod9 <- augment(spmod9)

plot(a_spmod9$.fitted, a_spmod9$.resid)
abline(h=0, col="red")

fitted_used <- fitted(spmod9, type = "response")
resid_dev <- residuals(spmod9, type = "response")
plot(fitted_used, resid_dev)
abline(h = 0, col = "red")

hist(a_spmod9$.resid, breaks=20, main="Residuals", xlab="Residuals")

coords <- as.matrix(grid_total_changes2[, c("x", "y")])
nb <- dnearneigh(coords, 0, 5000)  # create neighbors
lw <- nb2listw(nb, style="W")

moran.test(a_spmod9$.resid, lw)

moran.plot(a_spmod9$.resid, lw,
           xlab = "Delta Retama",
           ylab = "Promedio vecinal",
           main = "Moran scatterplot de % Retama")


#variograma datos
vg <- variogram(d_retama_area_3y ~ 1,
                locations = ~ x + y,
                data = grid_total_changes2)
plot(vg)

# variograma modelo
coords <- grid_total_changes2[, c("x", "y")]
resid_dev <- residuals(spmod9_g, type = "pearson")
df <- data.frame(x=coords$x, y=coords$y, resid=resid_dev)
coordinates(df) <- ~x+y
vg <- variogram(resid ~ 1, data = df)
plot(vg)


# Predicciones

newdata <- data.frame(
  d_urbano_area_3y = seq(
    min(grid_total_changes2$d_urbano_area_3y, na.rm = TRUE),
    max(grid_total_changes2$d_urbano_area_3y, na.rm = TRUE),
    length.out = 100
  ),
  x = mean(grid_total_changes2$x),
  y = mean(grid_total_changes2$y)
)

pred_new <- predict(spmod9_e, newdata = newdata, ,
                    interval = "prediction") %>% as.data.frame()
newdata$fit <- pred_new$fit
newdata$lower <- pred_new$lwr
newdata$upper <- pred_new$upr

set.seed(123)

muestra <- grid_total_changes2 %>%
  slice_sample(prop = 1)

#figura 3 
fig3 <-
ggplot(newdata, aes(x = d_urbano_area_3y/10000, 
                    y = (fit/1000000) * 100)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = (lower/1000000) * 100, 
        ymax = (upper/1000000) * 100),
    alpha = 0.2,
    fill = "#CDAA7D") +
  geom_point(aes(x = d_urbano_area_3y/10000, 
                 y = (d_retama_area_3y/1000000) * 100), 
             muestra, inherit.aes = F, alpha = 0.1) +
  xlab(expression(Delta * "Urbanización (ha)")) +
  ylab(expression(Delta * italic("C. scoparius") ~ " (ha)")) #+
  scale_y_continuous(breaks = seq(-20, 8, 4)) 
fig3

ggsave(filename = file.path(
       figure_dir, 
       "3_delta_retama_fx_delta_urbano.png"),
       fig3,
       width = 7,
       height = 5,
       bg = "white",
       dpi = 300)

grid_total_changes3 <- grid_subset %>%
  group_by(id) %>%
  mutate(d_retama_total_3y = (retama[year == 2024] + retama[year == 2023] + retama[year == 2022] )/ 3 - 
           (retama[year == 2017] + retama[year == 2018] + retama[year == 2019] )/ 3,
         d_retama_area_3y = (retama_area[year == 2024] + retama_area[year == 2023] + retama_area[year == 2022] )/ 3 - 
           (retama_area[year == 2017] + retama_area[year == 2018] + retama_area[year == 2019] )/ 3) %>%
  ungroup()

mapa_cambios <- grid_total_changes3 %>%
  mutate(
    cambio_retama = case_when(
      d_retama_area_3y > 0 ~ "Ganancia",
      d_retama_area_3y < 0 ~ "Pérdida",
      TRUE ~ "Sin cambio"
    )
  )

library(sf)
library(jsonlite)
library(dplyr)

geom <- geojsonsf::geojson_sfc(mapa_cambios$.geo)

mapa_cambios_sf <- st_sf(
  mapa_cambios %>% select(-.geo),
  geometry = geom
)

ggplot(mapa_cambios_sf) +
  geom_sf(aes(fill = cambio_retama), color = NA) +
  scale_fill_manual(
    values = c(
      "Ganancia" = "red",
      "Pérdida" = "blue",
      "Sin cambio" = "grey70"
    )
  ) +
  theme_void()
st_crs(mapa_cambios_sf)

spmod9_n <- readRDS(file.path(results_dir,"spmod9_n_retama_bin_grid1000.rds"))

st_write(
  mapa_cambios_sf,
  file.path(results_dir,"mapa_cambios.gpkg"),
  delete_dsn = TRUE
)

neg <- mapa_cambios$d_retama_area_3y[
  mapa_cambios$d_retama_area_3y < 0
]

pos <- mapa_cambios$d_retama_area_3y[
  mapa_cambios$d_retama_area_3y > 0
]

q_neg <- quantile(abs(neg), probs = c(1/3, 2/3), na.rm = TRUE)
q_pos <- quantile(pos, probs = c(1/3, 2/3), na.rm = TRUE)

mapa_cambios2 <- mapa_cambios %>%
  mutate(
    cambio_cat = case_when(
      d_retama_area_3y == 0 ~ "Sin cambio",
      
      d_retama_area_3y < -q_neg[2] ~ "Pérdida fuerte",
      d_retama_area_3y < -q_neg[1] ~ "Pérdida moderada",
      d_retama_area_3y < 0 ~ "Pérdida leve",
      
      d_retama_area_3y <= q_pos[1] ~ "Ganancia leve",
      d_retama_area_3y <= q_pos[2] ~ "Ganancia moderada",
      TRUE ~ "Ganancia fuerte"
    )
  )

geom <- geojsonsf::geojson_sfc(mapa_cambios$.geo)

mapa_cambios_sf2 <- st_sf(
  mapa_cambios2 %>% select(-.geo),
  geometry = geom
)

ggplot(mapa_cambios_sf2) +
  geom_sf(aes(fill = cambio_cat), color = NA) +
  scale_fill_manual(
    values = c(
      "Pérdida fuerte" = "#3300FF",
      "Pérdida moderada" = "#0099FF",
      "Pérdida leve" = "#67A9CF",
      "Sin cambio" = "#F7F7F7",
      "Ganancia leve" = "#FFBEB2",
      "Ganancia moderada" = "#F77964",
      "Ganancia fuerte" = "#AE123A"
    )
  ) +
  theme_void()

st_write(
  mapa_cambios_sf2,
  file.path(results_dir,"mapa_cambios_2_categorico.gpkg"),
  delete_dsn = TRUE
)
