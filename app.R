# --- Paquetes base y manipulación de datos ---
library(minpack.lm)
library(tidyverse)       # Colección de paquetes para manipulación, análisis y gráficos (incluye dplyr, ggplot2, tidyr, readr, stringr, etc.)
library(lubridate)       # Manejo de fechas y tiempos
library(data.table)      # Alternativa rápida a data.frames para datos grandes
library(janitor)         # Limpieza de datos (tablas, nombres de columnas, etc.)
library(stringr)         # Manipulación de strings
library(zip)             # Gestion de archivos comprimidos

# --- Spatial / raster / geospatial data ---
library(devtools)
library(raster)
library(sp)
library(sf)              # Lectura y manipulación de datos espaciales vectoriales (reemplaza rgdal, sp, rgeos en gran parte)
library(terra)           # Lectura y manipulación de datos raster (reemplazo moderno de raster, rgdal)
library(gdalUtilities)   # Ejecutar comandos GDAL desde R (útil para tareas geoespaciales avanzadas)
library(splines)
library(psych)

# --- NetCDF / climatología / meteorología ---
library(ncdf4)           # Lectura y escritura de archivos NetCDF

# --- Statistical analysis and modeling ---
library(pastecs)         # for statistics analysis "stat.desc"
library(gstat)           # Geoestadística, kriging, interpolación
library(latticeExtra)    # Gráficos en lattice (complemento útil para mapas y datos espaciales)
library(performance)     # Diagnósticos y métricas de modelos estadísticos
library(Hmisc)           # Estadísticas descriptivas, tablas, utilidades diversas
library(circular)        # Estadística para datos circulares (ángulos, direcciones)
library(matrixStats)     # Estadísticas rápidas sobre filas/columnas de matrices

# --- Visualization ---
library(ggpubr)          # Publicación de gráficos basados en ggplot2
library(ggprism)         # Estilo Prism (como en gráficos científicos)
library(corrplot)        # Visualización de matrices de correlación
library(RColorBrewer)    # Paletas de colores
library(colorspace)      # Paletas de colores avanzadas
library(plotrix)         # Funciones gráficas adicionales (p. ej. diagramas circulares, barras 3D)
library(factoextra)      # Visualización de análisis multivariados (PCA, clustering, etc.)
library(shiny)
library(shinyjs)
library(tmap)
library(leaflet)
library(plotly)

# --- Table output and reports ---
library(openxlsx)        # Leer y escribir archivos Excel sin depender de Java
library(readxl)          # Lectura de archivos Excel
library(DT)              # Tablas interactivas en HTML

# --- Parallel computing ---
library(parallel)        # Procesamiento en paralelo base de R
library(doParallel)      # Backend paralelo para foreach
library(foreach)         # Bucle foreach para paralelización

# --- Miscellaneous ---
library(reticulate)      # Integración con Python
library(tiff)            # Lectura y escritura de imágenes TIFF
library(here)      

tmap::ttm()

ui <- fluidPage(
  fluidRow(
    style = "background-color: #24618E; padding: 15px 20px; margin-bottom: 20px;",  # Fondo azul oscuro + algo de padding
    
    column(
      width = 4,  
      align = "left",
      style = "display: flex; align-items: center; gap: 20px;", 
      
      # Logo 1
      img(
        src = "GEO-Mountains-Logo.png",         
        height = "60px",            
        alt = "Logo GEO"
      ),
      
      # Logo 2
      img(
        src = "./logo-yt-blanco.png",
        height = "60px",
        alt = "Logo Yachay"
      )
    ),
    
    column(
      width = 9,  # El resto para el título
      align = "left",
      style = "display: flex; align-items: center;",  # Alinea verticalmente con los logos
      
      h2(
        "BALANCE DE MASA GLACIAR",
        style = "color: white; margin: 0; font-weight: bold;"  # Texto blanco, sin margen extra
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      tags$h3("Carga de datos"),
      
      fileInput("rawData", "Archivo rawData (*.txt o *.csv)", accept = c(".txt", ".csv")),
      fileInput("shapes", "Shapefiles (sube .shp, .dbf, .shx, .prj)", multiple = TRUE, accept = c(".shp",".dbf",".shx",".prj")),
      fileInput("dem", "Archivo DEM (*.tif)", accept = ".tif"),
      fileInput("geodetic", "Balance Geodésico (*.xlsx)", accept = ".xlsx"),
      
      numericInput("range_a", "Rango de cálculo", value = 50, min = 1),
      numericInput("resolution", "Resolución del raster", value = 100),
      textInput("crsID", "CRS (código EPSG, ej. 'EPSG:32717')", value = "EPSG:32717"),
      numericInput("NObs", "N.Obs (mínimo de observaciones por celda)", value = 6),
      numericInput("summit", "Altitud de la cima (aproximada) (m)", value = 5700),
      
      actionButton("run", "Ejecutar Modelo", class = "btn-primary"),
      hr(),
      downloadButton("downloadResultados", "Descargar resultados (.zip)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Progreso",
                 div(id = "log-panel",
                     style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px;
                            padding:15px; min-height:500px; font-family: monospace;
                            white-space: pre-wrap; overflow-y:auto;",
                     uiOutput("log_messages"))
        ),
        tabPanel("Balance de masa promedio anual",
                 br(),
                 leafletOutput("contour_map", height = "500px", width = "100%")
        ),
        tabPanel("Balance anual",
                 br(),
                 leafletOutput("balance_anual", height = "500px", width = "100%")
        ),
        tabPanel("Altura vs. Balance",
                 br(),
                 plotlyOutput("balance_altura", height = "500px")
        ),
        tabPanel("Altura vs. Balance por año",
                 br(),
                 plotlyOutput("balance_altura_anual", height = "500px")
        ),
        tabPanel(
          "Soporte instituciones",
          fluidRow(
            column(2, img(src = "Andes-C2H-logo_s.jpg",               height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "Autoridad-nacional-agua-ana-logo_s.jpg", height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "GEO-Mountains-logo_s.jpg",           height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "Glacioclim-logo_s.jpg",              height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "ianigla-logo_s.jpg",                 height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "IDEAM-logo_s.jpg",                   height = "120px", style = "display: block; margin: auto;"))
          ),
          
          fluidRow(
            column(2, img(src = "Inaigem-logo_s.jpg",                 height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "Red_ecuatoriana_cambio_climatico_s.png", height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "UMSA-logo_s.jpg",                    height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "Universidad-Austral-de-Chile-logo_s.jpg", height = "120px", style = "display: block; margin: auto;")),
            column(2, img(src = "wgms-logo_s.jpg",                    height = "120px", style = "display: block; margin: auto;")),
            column(2, "") 
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  tmap_mode("view")
  mapa_reactivo <- reactiveVal(NULL)
  mapa_reactivo2 <- reactiveVal(NULL)
  grafico_reactivo <- reactiveVal(NULL)
  grafico_reactivo2 <- reactiveVal(NULL)
  geodetic_used <- reactiveVal(FALSE)
  output$contour_map <- renderLeaflet({
    req(mapa_reactivo())     
    tmap::tmap_leaflet(mapa_reactivo())        
  })
  wb_guardado <- reactiveVal(NULL)
  
  output$balance_anual <- renderLeaflet({
    req(mapa_reactivo2())     
    tmap::tmap_leaflet(mapa_reactivo2())          
  })
  
  output$balance_altura <- renderPlotly({
    req(grafico_reactivo)
    grafico_reactivo()          
  })
  output$balance_altura_anual <- renderPlotly({
    req(grafico_reactivo2)
    grafico_reactivo2()          
  })
  
  # === SISTEMA DE LOGS ===
  log_messages <- reactiveVal("<span style='color:#666'>Esperando ejecución...</span><br>")
  
  add_log <- function(text, color = "#333") {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    new_line <- paste0("<span style='color:", color, "'>[", timestamp, "] ", text, "</span><br>")
    log_messages(paste0(log_messages(), new_line))
  }
  
  output$log_messages <- renderUI({ HTML(log_messages()) })
  
  # === BOTÓN EJECUTAR ===
  observeEvent(input$run, {
    
    req(input$rawData, input$shapes, input$dem, input$range_a)
    
    log_messages("<span style='color:#007bff;font-weight:bold'>Iniciando ejecución...</span><br>")
    add_log("Creando entorno temporal...", "#17a2b8")
    
    temp_dir <- tempdir()
    dir.create(file.path(temp_dir, "shapes"), showWarnings = FALSE, recursive = TRUE)
    dir_download <- file.path(temp_dir, "Resultados")
    dir.create(dir_download, showWarnings = FALSE, recursive = TRUE)
    # Copiar archivos
    raw_path <- file.path(temp_dir, input$rawData$name)
    file.copy(input$rawData$datapath, raw_path, overwrite = TRUE)
    
    dem_path <- file.path(temp_dir, input$dem$name)
    file.copy(input$dem$datapath, dem_path, overwrite = TRUE)
    
    geodetic_path <- NULL
    if (!is.null(input$geodetic)) {
      geodetic_path <- file.path(temp_dir, input$geodetic$name)
      file.copy(input$geodetic$datapath, geodetic_path, overwrite = TRUE)
      add_log("Balance geodésico: CARGADO", "#28a745")
    } else {
      add_log("Balance geodésico: NO proporcionado (opcional)", "#ffc107")
    }
    
    shape_dir <- file.path(temp_dir, "shapes")
    for(i in 1:nrow(input$shapes)){
      file.copy(input$shapes$datapath[i], file.path(shape_dir, input$shapes$name[i]), overwrite = TRUE)
    }
    
    add_log("Archivos guardados correctamente", "#28a745")
    add_log("Cargando función de Balance de masa...", "#17a2b8")
    
    source("functions/Funcion_Geoglacier.R", local = TRUE)
    
    # Workbook global para descarga
    wb <- openxlsx::createWorkbook()
    
    tryCatch({
      add_log("Ejecutando Balance de Masa Glaciar...", "#007bff")
      args_list <- list(
        rawData = raw_path,
        shapes_path = paste0(shape_dir, "/"),
        dem = dem_path,
        Resolution = input$resolution,
        crsID = input$crsID,
        N.Obs = input$NObs,
        summit = input$summit,
        geodetic_balance = geodetic_path,
        m = input$range_a,
        wb = wb,
        temp_dir = temp_dir
      )
      
      resultado <- do.call(Geoglacier, args_list)
      
      wb_guardado(wb)
      ruta_archivo_excel <- file.path(dir_download, "Resultados_Lliboutry.xlsx")
      openxlsx::saveWorkbook(wb, ruta_archivo_excel, overwrite = TRUE)
      
      if (!is.null(resultado$mapa_contornos_krg) && 
          inherits(resultado$mapa_contornos_krg, "tmap")) {
        
        mapa_reactivo(resultado$mapa_contornos_krg) 
        
        add_log("Mapa contornos recibido y mostrado correctamente", "#28a745")
        showNotification("¡Mapa mostrado!", type = "default", duration = 5)
        
      } else {
        mapa_reactivo(NULL)
        add_log("ADVERTENCIA: mapa_contornos es NULL o no es tmap", "#ffc107")
      }
      
      if (!is.null(resultado$mapa_raster) && 
          inherits(resultado$mapa_raster, "tmap")) {
        
        mapa_reactivo2(resultado$mapa_raster)   # ← GUARDAR
        
        add_log("Mapa raster recibido y mostrado correctamente", "#28a745")
        showNotification("¡Mapa mostrado!", type = "default", duration = 5)
        
      } else {
        mapa_reactivo2(NULL)
        add_log("ADVERTENCIA: mapa_contornos es NULL o no es tmap", "#ffc107")
      }
      
      if (!is.null(resultado$balance_altura) && 
          inherits(resultado$balance_altura, "plotly")) {
        
        grafico_reactivo(resultado$balance_altura)   # ← GUARDAR
        
        add_log("Grafico Dinámico recibido correctamente", "#28a745")
        showNotification("Gráfico Generado", type = "default", duration = 5)
        
      } else {
        grafico_reactivo(NULL)
        add_log("ADVERTENCIA: Grafico Nulo", "#ffc107")
      }
      
      if (!is.null(resultado$balance_altura_anual) && 
          inherits(resultado$balance_altura_anual, "plotly")) {
        
        grafico_reactivo2(resultado$balance_altura_anual)   # ← GUARDAR
        
        add_log("Grafico Dinámico anual recibido correctamente", "#28a745")
        showNotification("Gráfico anual Generado", type = "default", duration = 5)
        
      } else {
        grafico_reactivo2(NULL)
        add_log("ADVERTENCIA: Grafico Nulo", "#ffc107")
      }
      geodetic_used(as.logical(resultado$geodetic_use))
      add_log("Modelo ejecutado correctamente", "#28a745")
      add_log("EJECUCIÓN FINALIZADA", "#28a745")
      
      showNotification("¡Modelo ejecutado con éxito!", type = "default", duration = 8)
      
    }, error = function(e) {
      add_log(paste("ERROR:", e$message), "#dc3545")
      showNotification("Error en la ejecución del modelo", type = "error", duration = NULL)
    })
  })

  # === DESCARGA EXCEL ===
  output$downloadExcel <- downloadHandler(
    filename = "Resultados_Lliboutry.xlsx",
    content = function(dir_download) {
      # Usar el workbook guardado
      wb <- wb_guardado()
      
      # Validar que existe
      req(wb)
      
      # Guardar el workbook
      openxlsx::saveWorkbook(wb, dir_download, overwrite = TRUE)
    }
  )
  output$downloadResultados <- downloadHandler(
    filename = function() {
      paste0("Resultados_", format(Sys.Date(), "%Y%m%d"), ".zip")
    },
    content = function(file) {
      temp_dir <- tempdir()
      dir_resultados <- file.path(temp_dir, "Resultados")
      
      req(dir.exists(dir_resultados))
      
      # Cambiar al directorio temporal (funciona en todos los OS)
      wd_actual <- getwd()
      setwd(temp_dir)
      
      # Crear ZIP usando rutas relativas
      zip::zip(zipfile = file, files = "Resultados")
      
      # Volver al directorio original
      setwd(wd_actual)
    }
  )
}

shinyApp(ui, server)