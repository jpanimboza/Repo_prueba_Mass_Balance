ContourPlot_tmap <- function(df, method, title, brick_raw, shape_list, dem, crsID) {
  sf::sf_use_s2(FALSE)
  
  # -------------------------
  # 1) Obtener poly_sf robustamente
  # -------------------------
  build_poly_sf <- function(obj, crsID) {
    # Si ya es sf
    if (inherits(obj, "sf")) return(st_zm(obj))
    
    # sfc / sfg
    if (inherits(obj, "sfc") || inherits(obj, "sfg")) {
      return(st_zm(st_sf(geometry = obj)))
    }
    
    # SpatVector (terra)
    if (inherits(obj, "SpatVector")) {
      return(st_zm(st_as_sf(obj)))
    }
    
    # Spatial* (sp)
    if (inherits(obj, "Spatial")) {
      return(st_zm(st_as_sf(obj)))
    }
    
    # Lista con estructura de coordenadas (varias variantes)
    if (is.list(obj)) {
      # 1) buscar una matriz numérica dentro
      found_mat <- NULL
      # revisar elementos de primer nivel y segundo nivel
      candidates <- unlist(obj, recursive = FALSE)
      for (c in candidates) {
        if (is.matrix(c) && is.numeric(c)) { found_mat <- c; break }
        if (is.numeric(c) && is.vector(c) && length(c) > 4) {
          # posible vector con x...x,y...y (mitad/mitad)
          v <- c
          if ((length(v) %% 2) == 0) {
            n <- length(v) / 2
            mat_try <- matrix(v, ncol = 2, byrow = FALSE) # columnas = first half, second half
            if (is.numeric(mat_try)) { found_mat <- mat_try; break }
          }
        }
      }
      
      # Si no, intentar acceder a la ruta común [ [2] ][[1]]
      if (is.null(found_mat) && length(obj) >= 2) {
        maybe <- obj[[2]]
        if (is.list(maybe) && length(maybe) >= 1 && is.numeric(maybe[[1]])) {
          v <- maybe[[1]]
          if ((length(v) %% 2) == 0) {
            n <- length(v) / 2
            found_mat <- matrix(v, ncol = 2, byrow = FALSE)
          }
        }
      }
      
      if (!is.null(found_mat)) {
        # Asegurar que sea matrix n x 2
        mat <- as.matrix(found_mat)
        if (ncol(mat) != 2) {
          # si quedó transpuesta, intentar corregir
          if (nrow(mat) == 2) mat <- t(mat)
          if (ncol(mat) != 2) stop("No pude convertir las coordenadas a matriz n x 2.")
        }
        poly <- st_polygon(list(mat))
        sfobj <- st_sf(geometry = st_sfc(poly, crs = crsID))
        # si existe AREA en obj[[1]] agregarlo
        if (!is.null(obj[[1]]) && is.numeric(obj[[1]])) sfobj$AREA <- obj[[1]]
        return(st_zm(sfobj))
      }
    }
    
    stop("No pude convertir shape_list element a sf. Clase: ", paste(class(obj), collapse = ", "))
  }
  
  # Validaciones básicas
  if (missing(shape_list) || length(shape_list) == 0) stop("shape_list no está definido o está vacío.")
  if (is.null(shape_list[[1]])) stop("shape_list[[1]] es NULL.")
  
  # Construir polígono (usamos el primer elemento)
  poly_sf <- build_poly_sf(shape_list[[1]], crsID)
  
  # Convertir a terra vect para operaciones de máscara
  terra_poly <- terra::vect(poly_sf)
  
  # -------------------------
  # 2) Preparar DEM / grid (igual que en tu función original)
  # -------------------------
  if (missing(dem)) stop("Argumento 'dem' es requerido.")
  cropped_DEM <- raster(dem)
  cropped_DEM <- resample(cropped_DEM, brick_raw)
  cropped_DEM <- crop(cropped_DEM, brick_raw)
  spx <- as(cropped_DEM, "SpatialPixelsDataFrame")
  
  # -------------------------
  # 3) Puntos (df -> SpatialPoints)
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
  # 5) Convertir predicción a SpatRaster y aplicar máscara
  # -------------------------
  ai_predic_r <- tryCatch({
    terra::rast(ai_predic)
  }, error = function(e) stop("No pude convertir ai_predic a rast: ", e$message))
  
  ai_predic_r <- terra::mask(ai_predic_r["var1.pred"], terra_poly)
  # Comprobación de datos válidos en el raster
  vals <- terra::values(ai_predic_r)
  if (all(is.na(vals)) || length(vals[!is.na(vals)]) == 0) {
    stop("El raster interpolado quedó vacío (NA) tras aplicar máscara. Revisa CRS / extents / puntos.")
  }
  rng <- range(vals, na.rm = TRUE)
  if (diff(rng) == 0) {
    message("El raster tiene valor constante (sin variación). No se generarán contornos.")
  }
  
  # -------------------------
  # 6) Contornos robustos (intentos múltiples)
  # -------------------------
  ctr_sf <- NULL
  if (diff(rng) > 0.05) {
    levels_vec <- pretty(range(values(ai_predic_r), na.rm = TRUE), n = 12)
    message("DEBUG contour: niveles: ", paste(levels_vec, collapse=", "))
    
    # Asegurar poly_sf en CRS local (del raster)
    poly_sf_local <- tryCatch(st_transform(poly_sf, crs = crsID), error = function(e) poly_sf)
    
    # Helper: convertir SpatVector/Spatial*-> sf y limpiar
    clean_to_sf <- function(obj) {
      if (is.null(obj)) return(NULL)
      sfobj <- tryCatch({
        st_as_sf(obj)
      }, error = function(e) {
        return(NULL)
      })
      if (is.null(sfobj) || nrow(sfobj) == 0) return(NULL)
      # forzar numeric level si existe
      if ("level" %in% names(sfobj)) {
        sfobj$level <- suppressWarnings(as.numeric(as.character(sfobj$level)))
        sfobj <- sfobj[!is.na(sfobj$level), ]
      }
      # casteos geométricos
      if (any(st_geometry_type(sfobj) == "MULTILINESTRING"))
        sfobj <- tryCatch(st_cast(sfobj, "LINESTRING"), error = function(e) sfobj)
      sfobj <- tryCatch(st_make_valid(sfobj), error = function(e) sfobj)
      return(sfobj)
    }
    
    # 1) Intento directo con terra::contour
    ctr_try <- tryCatch(terra::contour(ai_predic_r, levels = levels_vec), error = function(e) NULL)
    if (!is.null(ctr_try) && nrow(ctr_try) > 0) {
      ctr_sf <- clean_to_sf(ctr_try)
      message("terra::contour OK, geometrías: ", ifelse(is.null(ctr_sf), "NULL", paste(unique(st_geometry_type(ctr_sf)), collapse=",")))
    } else {
      message("terra::contour devolvió NULL o 0 geometrías.")
    }
    
    # 2) Intento con raster::rasterToContour (a veces más tolerante)
    if (is.null(ctr_sf)) {
      ctr_try2 <- tryCatch({
        rtmp <- raster::raster(ai_predic_r) # convertir a raster package
        raster::rasterToContour(rtmp, levels = levels_vec)
      }, error = function(e) NULL)
      if (!is.null(ctr_try2)) {
        ctr_sf <- clean_to_sf(ctr_try2)
        message("rasterToContour OK, geometrías: ", ifelse(is.null(ctr_sf), "NULL", paste(unique(st_geometry_type(ctr_sf)), collapse=",")))
      } else {
        message("rasterToContour falló o devolvió NULL.")
      }
    }
    # 3) Si se generaron contornos, transformar, recortar y limpiar final
    if (!is.null(ctr_sf) && nrow(ctr_sf) > 0) {
      # asegurarse CRS correcto
      ctr_sf <- tryCatch(st_transform(ctr_sf, crs = crsID), error = function(e) ctr_sf)
      # recortar al polígono en CRS local
      ctr_sf <- tryCatch(st_intersection(ctr_sf, poly_sf_local), error = function(e) ctr_sf)
      # eliminar rows vacías si las hay
      ctr_sf <- ctr_sf[!st_is_empty(ctr_sf$geometry), , drop = FALSE]
      # ultimo chequeo
      message("Contornos finales: filas=", nrow(ctr_sf), " tipos=", paste(unique(st_geometry_type(ctr_sf)), collapse=","))
    } else {
      message("No se pudieron generar contornos tras todos los intentos.")
    }
  }
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
  # 8) Construir mapa (tmap v4)
  # -------------------------
  
  map <- tm_shape(ai_predic_r) +
    tm_raster(
      col.scale = tm_scale_continuous(values = rampcols, midpoint = 0),
      col.legend = tm_legend(title = expression("Balance (m w.e. a"^-1*")")),
    )+
    tm_shape(poly_sf) + tm_borders(lwd = 3) +
    {      if (!is.null(ctr_sf) && nrow(ctr_sf)>0) {
        tm_shape(ctr_sf) + tm_lines(lwd = 2.2)
      } else {
        NULL # para que tmap no falle
      }
    } +
    tm_title(title) +
    tm_layout(legend.outside = TRUE,
              legend.outside.position = "right",
              frame = FALSE) +
    tm_mouse_coordinates()
  return(map)
}