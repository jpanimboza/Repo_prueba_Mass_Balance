#===================== Function that loads all the libraries ===================
# Package names (updated and modern)
packages <- c(
  # --- General purpose data manipulation and analysis ---
  "minpack.lm",
  "tidyverse",    # Colección de paquetes para manipulación, análisis y gráficos (incluye dplyr, ggplot2, tidyr, readr, stringr, etc.)
  "lubridate",    # Manejo de fechas y tiempos
  "data.table",   # Alternativa rápida a data.frames para datos grandes
  "janitor",      # Limpieza de datos (tablas, nombres de columnas, etc.)
  "stringr",      # Manipulación de strings
  "zip",          #Gestion de archivos comprimidos
  
  # --- Spatial / raster / geospatial data ---
  "devtools",
  "raster",
  "sp",
  "sf",           # Lectura y manipulación de datos espaciales vectoriales (reemplaza rgdal, sp, rgeos en gran parte)
  "terra",        # Lectura y manipulación de datos raster (reemplazo moderno de raster, rgdal)
  "gdalUtilities",# Ejecutar comandos GDAL desde R (útil para tareas geoespaciales avanzadas)
  "splines",
  "psych",
  
  # --- NetCDF / climatología / meteorología ---
  "ncdf4",        # Lectura y escritura de archivos NetCDF
  
  # --- Statistical analysis and modeling ---
  "pastecs",      # for statistics analysis "stat.desc"
  "gstat",        # Geoestadística, kriging, interpolación
  "latticeExtra", # Gráficos en lattice (complemento útil para mapas y datos espaciales)
  "performance",  # Diagnósticos y métricas de modelos estadísticos
  "Hmisc",        # Estadísticas descriptivas, tablas, utilidades diversas
  "circular",     # Estadística para datos circulares (ángulos, direcciones)
  "matrixStats",  # Estadísticas rápidas sobre filas/columnas de matrices
  
  # --- Visualization ---
  "ggpubr",       # Publicación de gráficos basados en ggplot2
  "ggprism",      # Estilo Prism (como en gráficos científicos)
  "corrplot",     # Visualización de matrices de correlación
  "RColorBrewer", # Paletas de colores
  "colorspace",   # Paletas de colores avanzadas
  "plotrix",      # Funciones gráficas adicionales (p. ej. diagramas circulares, barras 3D)
  "factoextra",   # Visualización de análisis multivariados (PCA, clustering, etc.)
  "shiny",
  "shinyjs",
  "tmap",
  "leaflet",
  "plotly",
  
  # --- Table output and reports ---
  "openxlsx",     # Leer y escribir archivos Excel sin depender de Java
  "readxl",       # Lectura de archivos Excel
  "DT",           # Tablas interactivas en HTML
  
  # --- Parallel computing ---
  "parallel",     # Procesamiento en paralelo base de R
  "doParallel",   # Backend paralelo para foreach
  "foreach",      # Bucle foreach para paralelización
  
  # --- Miscellaneous ---
  "reticulate",   # Integración con Python
  "tiff",         # Lectura y escritura de imágenes TIFF
  "here"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#install.packages("https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.6-7.tar.gz", repos = NULL, type = "source")