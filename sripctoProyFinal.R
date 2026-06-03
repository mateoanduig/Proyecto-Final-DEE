
#  CARGAR LIBRERÍAS
library(tidyverse)
library(readxl)
library(sf)
library(mapSpain)
library(spdep)
library(rgeoda)
library(tmap)
library(leaflet)
library(knitr)
library(kableExtra)

# LECTURA Y LIMPIEZA DE DATOS (CSV)

#  Índice de Alquiler (59058.csv)
alquiler_raw <- read_delim("59058.csv", delim = ";", 
                           locale = locale(encoding = "latin1", decimal_mark = ","),
                           show_col_types = FALSE)

alquiler_clean <- alquiler_raw %>%
  select(Provincias, Total) %>%
  rename(indice_alquiler = Total) %>%
  mutate(
    Provincias = ifelse(is.na(Provincias), "00 Total Nacional", Provincias),
    cpro = str_sub(Provincias, 1, 2)
  ) %>%
  select(cpro, provincia = Provincias, indice_alquiler)

# Pernoctaciones Turísticas (67183.csv)
turismo_raw <- read_delim("67183.csv", delim = ";", 
                          locale = locale(encoding = "latin1", grouping_mark = "."),
                          show_col_types = FALSE)

turismo_clean <- turismo_raw %>%
  select(Provincias, `País de residencia`, Total) %>%
  mutate(cpro = str_sub(Provincias, 1, 2)) %>%
  pivot_wider(names_from = `País de residencia`, values_from = Total) %>%
  rename(
    turismo_total = Total,
    turismo_españoles = Españoles,
    turismo_extranjeros = Extranjeros
  ) %>%
  select(cpro, turismo_total, turismo_españoles, turismo_extranjeros)

# Unión de bases de datos
datos_finales <- alquiler_clean %>%
  left_join(turismo_clean, by = "cpro")

#  PREPARACIÓN ESPACIAL (ESPAÑA PENINSULAR)

# Códigos de territorios extrapeninsulares a excluir (Baleares, Canarias, Ceuta, Melilla)
codigos_islas_ciudades <- c("07", "35", "38", "51", "52")

# Descargamos polígonos, proyectamos a UTM (EPSG:32630) y unimos datos
mapa_peninsula <- esp_get_prov(moveCAN = FALSE) %>%
  st_transform(32630) %>%
  left_join(datos_finales, by = "cpro") %>%
  filter(!(cpro %in% codigos_islas_ciudades)) %>%
  filter(!is.na(indice_alquiler))


# EXPLORACIÓN VISUAL
# ------------------------------------------------------------------------------
library(patchwork)

mapa_alq <- ggplot(mapa_peninsula) +
  geom_sf(aes(fill = indice_alquiler), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  theme_void() +
  labs(title = "Índice de Precios de Alquiler", fill = "Índice")

mapa_tur <- ggplot(mapa_peninsula) +
  geom_sf(aes(fill = turismo_total), color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10")+ 
  theme_void() +
  labs(title = "Presión Turística Total", fill = "Pernoctaciones")

# Ensamblar y mostrar
mapas_descriptivos <- mapa_alq + mapa_tur +
  plot_annotation(title = "Distribución Espacial de las Variables de Estudio")
print(mapas_descriptivos)


#  GRÁFICO DE CORRELACIÓN SIMPLE (APARTADO 3.1.3)

library(ggrepel)

correlacion <- cor(mapa_peninsula$turismo_total, mapa_peninsula$indice_alquiler)
provincias_clave <- "Valencia|Málaga|Madrid|Alicante|Barcelona|Zamora|León|Ourense|Soria|Segovia|Granada|Cádiz|Ciudad Real|Albacete"

mapa_peninsula <- mapa_peninsula %>%
  mutate(etiqueta_destacada = ifelse(str_detect(provincia, regex(provincias_clave, ignore_case = TRUE)), provincia, NA))

grafico_correlacion <- ggplot(mapa_peninsula, aes(x = turismo_total, y = indice_alquiler)) +
  geom_point(color = "darkblue", alpha = 0.7, size = 2.5) +
  geom_text_repel(aes(label = etiqueta_destacada), size = 3.5, box.padding = 0.6, point.padding = 0.2, color = "black", fontface = "bold", max.overlaps = Inf, na.rm = TRUE) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.3) +
  scale_x_log10(labels = scales::comma_format(big.mark = ".", decimal.mark = ",")) +
  annotate("text", x = 50000, y = max(mapa_peninsula$indice_alquiler), label = paste("Correlación (r) =", round(correlacion, 3)), size = 4, color = "black", fontface = "bold") +
  theme_minimal() +
  labs(title = "Relación entre Presión Turística y Precio del Alquiler", x = "Pernoctaciones de Turistas (Escala Logarítmica)", y = "Índice de Precios de Vivienda en Alquiler")

print(grafico_correlacion)


#  MATRIZ DE CONTIGÜIDAD (CRITERIO REINA)

# Conecta provincias si comparten al menos 1 punto de frontera
vecinos_reina <- poly2nb(mapa_peninsula)

# Matriz estandarizada por filas (W)
matriz_pesos_peninsula <- nb2listw(vecinos_reina, style = "W")
plot(vecinos_reina, st_centroid(st_geometry(mapa_peninsula)), add = TRUE, col = "blue4", lwd = 1.5)

# AUTOCORRELACIÓN GLOBAL (ÍNDICE DE MORAN - MONTE CARLO)

set.seed(123) # Para reproducibilidad

# Moran - Índice de Alquiler
moran_peninsula_alq <- moran.mc(mapa_peninsula$indice_alquiler, matriz_pesos_peninsula, nsim = 999)
print("--- MORAN: ÍNDICE DE ALQUILER ---")
print(moran_peninsula_alq)

# Moran - Turismo Extranjero
moran_peninsula_tur <-moran.mc(mapa_peninsula$turismo_total, matriz_pesos_peninsula, nsim = 999)

print("--- MORAN: TURISMO EXTRANJERO ---")
print(moran_peninsula_tur)

#  AUTOCORRELACIÓN LOCAL (LISA) Y PREPARACIÓN DE MAPAS

#  Cálculo LISA
lisa_alquiler <- localmoran(mapa_peninsula$indice_alquiler, matriz_pesos_peninsula)
mapa_peninsula$lisa_i <- lisa_alquiler[, "Ii"]
mapa_peninsula$lisa_p <- lisa_alquiler[, "Pr(z != E(Ii))"] 

# Clasificación Hot-Spots y Cold-Spots
z_alquiler <- as.numeric(scale(mapa_peninsula$indice_alquiler))
w_z_alquiler <- lag.listw(matriz_pesos_peninsula, z_alquiler)
alpha <- 0.05

mapa_peninsula$cluster_lisa <- "No significativo"
mapa_peninsula$cluster_lisa[z_alquiler > 0 & w_z_alquiler > 0 & mapa_peninsula$lisa_p <= alpha] <- "Alto-Alto (Hot-spot)"
mapa_peninsula$cluster_lisa[z_alquiler < 0 & w_z_alquiler < 0 & mapa_peninsula$lisa_p <= alpha] <- "Bajo-Bajo (Cold-spot)"
mapa_peninsula$cluster_lisa[z_alquiler > 0 & w_z_alquiler < 0 & mapa_peninsula$lisa_p <= alpha] <- "Alto-Bajo (Outlier)"
mapa_peninsula$cluster_lisa[z_alquiler < 0 & w_z_alquiler > 0 & mapa_peninsula$lisa_p <= alpha] <- "Bajo-Alto (Outlier)"

mapa_peninsula$cluster_lisa <- factor(mapa_peninsula$cluster_lisa, 
                                      levels = c("Alto-Alto (Hot-spot)", "Bajo-Bajo (Cold-spot)", 
                                                 "Alto-Bajo (Outlier)", "Bajo-Alto (Outlier)", "No significativo"))

#  Variable categórica de significatividad para el Mapa 2
mapa_peninsula$nivel_sig <- cut(mapa_peninsula$lisa_p, 
                                breaks = c(0, 0.01, 0.05, 1), 
                                labels = c("Alta (p < 0.01)", "Normal (p < 0.05)", "No significativo"),
                                include.lowest = TRUE)

# 1. Para ver cuántas provincias hay en cada grupo:
table(mapa_peninsula$cluster_lisa)

# 2. Para ver el nombre exacto de las provincias significativas:
mapa_peninsula %>% 
  filter(cluster_lisa != "No significativo") %>% 
  select(provincia, cluster_lisa) %>% 
  st_drop_geometry() %>% 
  arrange(cluster_lisa)
#  CARTOGRAFÍA FINAL LISA (tmap)
tmap_mode("plot")

# Mapa de clústeres
mapa_cluster <- tm_shape(mapa_peninsula) +
  tm_polygons("cluster_lisa", palette = c("red", "blue", "pink", "lightblue", "white"), 
              title = "Clústeres (Alquiler)", border.col = "black", lwd = 0.5, colorNA = "white") +
  tm_layout(main.title = "Mapa LISA", main.title.size = 0.9, main.title.position = "center",
            legend.outside = TRUE, legend.outside.position = "right", frame = FALSE)

# Mapa de p-valores
mapa_sig <- tm_shape(mapa_peninsula) +
  tm_polygons("nivel_sig", palette = c("#1a9641", "#a6d96a", "white"), 
              title = "p-valor", border.col = "black", lwd = 0.5) +
  tm_layout(main.title = "Significatividad", main.title.size = 0.9, main.title.position = "center",
            legend.outside = TRUE, legend.outside.position = "right", frame = FALSE)

# Visualización conjunta
tmap_arrange(mapa_cluster, mapa_sig, ncol = 2)



# ------------------------------------------------------------------------------
#  GRÁFICOS DE DISPERSIÓN DE MORAN
# ------------------------------------------------------------------------------
z_alquiler <- as.numeric(scale(mapa_peninsula$indice_alquiler))
z_turismo <- as.numeric(scale(log10(mapa_peninsula$turismo_total)))

wz_alquiler <- lag.listw(matriz_pesos_peninsula, z_alquiler)
wz_turismo <- lag.listw(matriz_pesos_peninsula, z_turismo)

df_moran <- data.frame(provincia = mapa_peninsula$provincia, z_alquiler, wz_alquiler, z_turismo, wz_turismo)

graf_moran_alq <- ggplot(df_moran, aes(x = z_alquiler, y = wz_alquiler)) +
  geom_point(color = "darkblue", alpha = 0.6) + geom_smooth(method = "lm", color = "red", fill = "pink") +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() + labs(title = "Panel A: Alquiler", x = "Índice Estandarizado", y = "Retardo Espacial")

graf_moran_tur <- ggplot(df_moran, aes(x = z_turismo, y = wz_turismo)) +
  geom_point(color = "darkgreen", alpha = 0.6) + geom_smooth(method = "lm", color = "red", fill = "pink") +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() + labs(title = "Panel B: Turismo Total", x = "Log(Pernoctaciones) Est.", y = "Retardo Espacial")

print(graf_moran_alq + graf_moran_tur + plot_annotation(title = "Gráficos de Dispersión de Moran"))


#  MODELIZACIÓN ECONOMÉTRICA (OLS, DIAGNÓSTICO Y SAR)
library(spatialreg)

# 1. Creamos una nueva variable escalada (Turistas en millones)
mapa_peninsula$turismo_total_millones <- mapa_peninsula$turismo_total / 1000000

# 2. Modelo Lineal Base (OLS) con la variable escalada
modelo_ols <- lm(indice_alquiler ~ turismo_total_millones, data = mapa_peninsula)
print("--- MODELO OLS ---")
summary(modelo_ols)

# 3. Diagnóstico Espacial (LM Tests)
tests_lm <- lm.LMtests(modelo_ols, matriz_pesos_peninsula, test = "all")
print("--- TESTS DE MULTIPLICADORES DE LAGRANGE ---")
print(tests_lm)

# 4. Modelo de Retardo Espacial (SAR) con la variable escalada
modelo_sar <- lagsarlm(indice_alquiler ~ turismo_total_millones, data = mapa_peninsula, listw = matriz_pesos_peninsula)
print("--- MODELO SAR ---")
summary(modelo_sar)

#  REGRESIÓN GEOGRÁFICAMENTE PONDERADA (GWR)
# ------------------------------------------------------------------------------
library(spgwr)

# Convertimos a formato espacial clásico (sp) necesario para spgwr
mapa_sp <- as_Spatial(mapa_peninsula)

# Calculamos el ancho de banda óptimo (Bandwidth)
bw_optimo <- gwr.sel(indice_alquiler ~ turismo_total, data = mapa_sp, coords = coordinates(mapa_sp))

# Ejecutamos el modelo GWR
modelo_gwr <- gwr(indice_alquiler ~ turismo_total, data = mapa_sp, coords = coordinates(mapa_sp), 
                  bandwidth = bw_optimo, hatmatrix = TRUE)

print("--- MODELO GWR ---")
print(modelo_gwr)




##
library(tmap)
tmap_mode("plot")

# Mapa 1: Impacto Local
mapa_gwr_coef <- tm_shape(mapa_peninsula) +
  tm_polygons(fill = "coef_turismo_local", 
              fill.scale = tm_scale_intervals(style = "quantile", n = 5, values = "YlOrRd"),
              fill.legend = tm_legend(title = "Coef. Local"),
              col = "black", lwd = 0.3) +
  tm_title("GWR: Impacto Local", size = 0.9) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "bottom",
            frame = FALSE)                          

# Mapa 2: Capacidad Explicativa
mapa_gwr_r2 <- tm_shape(mapa_peninsula) +
  tm_polygons(fill = "r2_local", 
              fill.scale = tm_scale_intervals(style = "pretty", values = "Blues"),
              fill.legend = tm_legend(title = "R-cuadrado"),
              col = "black", lwd = 0.3) +
  tm_title("GWR: Capacidad Explicativa", size = 0.9) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "bottom",
            frame = FALSE)

# Juntamos los mapas
mapas_juntos <- tmap_arrange(mapa_gwr_coef, mapa_gwr_r2, ncol = 2)

# Exportar directamente al archivo grande (14x7 pulgadas)
#tmap_save(mapas_juntos, filename = "mapa_GWR_definitivo.png", width = 14, height = 7, units = "in", dpi = 300)


