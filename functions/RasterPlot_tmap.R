RasterPlot_tmap <- function(obj, showraw, title, brick_raw, Nyear, shape_list, crsID) {
  
  # 1. Raster multicapa
  r <- rast(obj)
  Nyear <- as.Date(Nyear)
  
  # 2. Extraer años
  if (inherits(Nyear, "Date")) {
    years <- as.integer(format(Nyear, "%Y"))
  } else {
    years <- as.integer(Nyear)
  }
  
  # 3. Nombres de capas
  names(r) <- paste0(years, "-", years + 1)
  
  # 4. Shapefile
  poly_sf <- if (inherits(shape_list[[1]], "sf")) {
    st_zm(shape_list[[1]])
  } else {
    st_as_sf(st_zm(shape_list[[1]]))
  }
  
  # 5. Paleta divergente centrada en 0 (global para todas las capas)
  Min_val <- floor(min(values(r), na.rm = TRUE))
  Max_val <- ceiling(max(values(r), na.rm = TRUE))
  nHalf <- 1500
  cp <- brewer.pal(11, "RdBu")
  rc1 <- colorRampPalette(c(cp[1], cp[3], cp[6]))(nHalf)
  rc2 <- colorRampPalette(c(cp[6], cp[9], cp[11]))(nHalf)
  rampcols <- c(rc1, rc2)
  
  map <- tm_shape(r) +
    tm_raster(
      style = "cont",  # Estilo continuo explícito
      palette = rampcols,  # Paleta personalizada (usa midpoint=0 por defecto en cont)
      midpoint = 0,  # Centrado en 0
      limits = c(Min_val, Max_val),  # Rango global fijo para todas las capas
      title = "Valor",  # Título de la leyenda única
      alpha = 1
    ) +
    tm_facets(
      as.layers = TRUE, 
      sync = TRUE,
      free.scales = FALSE 
    ) +
    tm_shape(poly_sf) +
    tm_borders(lwd = 2, col = "black") +
    tm_title(title) +
    tm_layout(
      legend.outside = FALSE, 
      frame = FALSE,
      legend.position = c("right", "bottom")  
    )
  
  return(map)
}