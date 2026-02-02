ContourPlot_tmap <- function(df, method, title, brick_raw, shape_list, dem, crsID) {
  sf::sf_use_s2(FALSE)
  
  # -------------------------
  # 1) Obtener poly_sf robustamente
  # -------------------------
  build_poly_sf <- function(obj, crsID) {
    # ... (tu código existente) ...
    if (inherits(obj, "sf")) return(st_zm(obj))
    if (inherits(obj, "sfc") || inherits(obj, "sfg")) {
      return(st_zm(st_sf(geometry = obj)))
    }
    if (inherits(obj, "SpatVector")) {
      return(st_zm(st_as_sf(obj)))
    }
    if (inherits(obj, "Spatial")) {
      return(st_zm(st_as_sf(obj)))
    }
    # ... resto del código ...
    stop("No pude convertir shape_list element a sf. Clase: ", paste(class(obj), collapse = ", "))
  }
  
  # Validaciones básicas
  if (missing(shape_list) || length(shape_list) == 0) stop("shape_list no está definido o está vacío.")
  if (is.null(shape_list[[1]])) stop("shape_list[[1]] es NULL.")
  
  poly_sf <- build_poly_sf(shape_list[[1]], crsID)
  terra_poly <- terra::vect(poly_sf)
  
  # -------------------------
  # 2) Preparar DEM / grid
  # -------------------------
  if (missing(dem)) stop("Argumento 'dem' es requerido.")
  cropped_DEM <- raster(dem)
  cropped_DEM <- resample(cropped_DEM, brick_raw)
  cropped_DEM <- crop(cropped_DEM, brick_raw)
  spx <- as(cropped_DEM, "SpatialPixelsDataFrame")
  
  # -------------------------
  # 3) Puntos
  # -------------------------
  points <- na.omit(df)
  names(points) <- c("x", "y", "ai")
  coordinates(points) <- ~ x + y
  crs(points) <- crsID
  
  # -------------------------
  # 4) Interpolación
  # -------------------------
  if (method == "idw") {
    ai_predic <- idw(ai ~ 1, locations = points, newdata = spx, idp = 2)
  } else if (method == "krig") {
    Var <- variogram(ai ~ x + y, points, cutoff = 1500)
    fit_krig <- fit.variogram(
      Var,
      vgm(c("Exp", "Gau", "Sph"),
          psill = mean(Var$gamma, na.rm = TRUE),
          range = max(Var$dist, na.rm = TRUE) / 3,
          nugget = min(Var$gamma, na.rm = TRUE))
    )
    ai_predic <- krige(ai ~ x + y, points, newdata = spx, model = fit_krig)
  } else {
    stop("Método no reconocido. Use 'idw' o 'krig'.")
  }
  
  # -------------------------
  # 5) Convertir y aplicar máscara
  # -------------------------
  ai_predic_r <- tryCatch({
    terra::rast(ai_predic)
  }, error = function(e) stop("No pude convertir ai_predic a rast: ", e$message))
  
  ai_predic_r <- terra::mask(ai_predic_r["var1.pred"], terra_poly)
  
  vals <- terra::values(ai_predic_r)
  if (all(is.na(vals)) || length(vals[!is.na(vals)]) == 0) {
    stop("El raster interpolado quedó vacío (NA) tras aplicar máscara.")
  }
  rng <- range(vals, na.rm = TRUE)
  
  # -------------------------
  # 6) Contornos
  # -------------------------
  invisible({
    ctr_sf <- NULL
    if (diff(rng) > 0.05) {
      levels_vec <- pretty(range(values(ai_predic_r), na.rm = TRUE), n = 12)
      
      poly_sf_local <- tryCatch(st_transform(poly_sf, crs = crsID), error = function(e) poly_sf)
      
      clean_to_sf <- function(obj) {
        if (is.null(obj)) return(NULL)
        sfobj <- tryCatch(st_as_sf(obj), error = function(e) NULL)
        if (is.null(sfobj) || nrow(sfobj) == 0) return(NULL)
        if ("level" %in% names(sfobj)) {
          sfobj$level <- suppressWarnings(as.numeric(as.character(sfobj$level)))
          sfobj <- sfobj[!is.na(sfobj$level), ]
        }
        if (any(st_geometry_type(sfobj) == "MULTILINESTRING"))
          sfobj <- tryCatch(st_cast(sfobj, "LINESTRING"), error = function(e) sfobj)
        sfobj <- tryCatch(st_make_valid(sfobj), error = function(e) sfobj)
        return(sfobj)
      }
      invisible({
        ctr_try <- tryCatch(terra::contour(ai_predic_r, levels = levels_vec), error = function(e) NULL)
        if (!is.null(ctr_try) && nrow(ctr_try) > 0) {
          ctr_sf <- suppressWarnings(clean_to_sf(ctr_try))
        }
        
        if (is.null(ctr_sf)) {
          ctr_try2 <- tryCatch({
            rtmp <- raster::raster(ai_predic_r)
            raster::rasterToContour(rtmp, levels = levels_vec)
          }, error = function(e) NULL)
          if (!is.null(ctr_try2)) {
            ctr_sf <- suppressWarnings(clean_to_sf(ctr_try2))
          }
        }
      })
      if (!is.null(ctr_sf) && nrow(ctr_sf) > 0) {
        ctr_sf <- suppressWarnings({
          ctr_sf <- st_transform(ctr_sf, crs = crsID)
          ctr_sf <- st_intersection(ctr_sf, poly_sf_local)
          ctr_sf[!st_is_empty(ctr_sf$geometry), , drop = FALSE]
        })
      }
    }
  })
  # -------------------------
  # 7) Paleta y breaks
  # -------------------------
  Thresh <- 0
  nHalf <- 6
  Min <- floor(min(ai_predic$var1.pred, na.rm = TRUE))
  Max <- ceiling(max(ai_predic$var1.pred, na.rm = TRUE))
  
  rb1 <- seq(Min, Thresh, length.out = nHalf + 1)
  rb2 <- seq(Thresh, Max, length.out = nHalf + 1)[-1]
  rampbreaks <- c(rb1, rb2)
  
  cp <- brewer.pal(11, "RdBu")
  rc1 <- colorRampPalette(c(cp[1], cp[3], cp[6]))(nHalf)
  rc2 <- colorRampPalette(c(cp[6], cp[9], cp[11]))(nHalf)
  rampcols <- c(rc1, rc2)
  
  # -------------------------
  # 8) Construir mapa (compatible v3 y v4)
  # -------------------------
  tmap_ver <- packageVersion("tmap")
  map <- tryCatch({
    if (tmap_ver >= "4.0") {
      # Sintaxis tmap v4.2 (tu código original)
      m <- tm_shape(ai_predic_r) +
        tm_raster(
          col.scale = tm_scale_continuous(values = rampcols, midpoint = 0),
          col.legend = tm_legend(title = "Balance (m w.e./año)")
        ) +
        tm_shape(poly_sf) + 
        tm_borders(lwd = 3)
      
      if (!is.null(ctr_sf) && nrow(ctr_sf) > 0) {
        m <- m + tm_shape(ctr_sf) + tm_lines(lwd = 2.2)
      }
      
      m <- m + 
        tm_layout(
          legend.outside = TRUE,
          legend.outside.position = "right",
          frame = FALSE
        ) +
        tm_mouse_coordinates()
      
    } else {
      # Sintaxis tmap v3.3.4 (para Kaggle/notebook)
      m <- tm_shape(ai_predic_r) +
        tm_raster(
          "var1.pred",
          palette = rampcols,
          breaks = rampbreaks,
          midpoint = 0,
          title = "Balance (m w.e./año)",
          style = "fixed"
        ) +
        tm_shape(poly_sf) + 
        tm_borders(lwd = 3, col = "black")
      
      if (!is.null(ctr_sf) && nrow(ctr_sf) > 0) {
        m <- m + tm_shape(ctr_sf) + tm_lines(col = "black", lwd = 2.2)
      }
      
      m <- m + 
        tm_layout(
          title = title,
          legend.outside = TRUE,
          legend.outside.position = "right",
          frame = FALSE
        ) +
        tm_mouse_coordinates()
    }
  }, error = function(e) {
    stop("Error al construir mapa tmap: ", e$message)
  })
  dev.off()
  return(map)

}
