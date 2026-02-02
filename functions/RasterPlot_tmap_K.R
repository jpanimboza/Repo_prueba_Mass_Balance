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
  
  # 6. Detectar versión de tmap
  tmap_ver <- packageVersion("tmap")
  
  # 7. Crear breaks para v3
  breaks_vec <- seq(Min_val, Max_val, length.out = 100)
  
  # 8. Construir mapa según versión
  # 8. Construir mapa según versión
  map <- tryCatch({
    if (tmap_ver >= "4.0") {
      m <- tm_shape(r) +
        tm_raster(
          col.scale = tm_scale_continuous(
            values = rampcols,
            midpoint = 0,
            limits = c(Min_val, Max_val)
          ),
          col.legend = tm_legend(
            title = "Balance (m w.e./año)",
            show = c(TRUE, rep(FALSE, nlyr(r) - 1))  # ← Solo primera capa
          )
        ) +
        tm_facets(
          as.layers = TRUE,
          sync = TRUE,
          free.scales = FALSE
        ) +
        tm_shape(poly_sf) +
        tm_borders(lwd = 2, col = "black") +
        tm_layout(
          title = title,
          legend.outside = FALSE,
          frame = FALSE,
          legend.position = c("right", "bottom")
        )
      
    } else {   
      # Crear array de TRUE/FALSE para legend.show
      legend_array <- c(TRUE, rep(FALSE, nlyr(r) - 1))
      
      m <- tm_shape(r) +
        tm_raster(
          style = "cont",
          palette = rampcols,
          midpoint = 0,
          breaks = breaks_vec,
          title = "Balance (m w.e./año)",
          legend.show = legend_array  # ← Solo mostrar leyenda en primera capa
        ) +
        tm_facets(
          as.layers = TRUE,
          sync = TRUE,
          free.scales = FALSE
        ) +
        tm_shape(poly_sf) +
        tm_borders(lwd = 2, col = "black") +
        tm_layout(
          title = title,
          legend.outside = FALSE,
          frame = FALSE,
          legend.position = c("right", "bottom")
        )
    }
    
    m
  }, error = function(e) {
    stop("Error al construir mapa raster: ", e$message)
  })
  
  return(map)
}