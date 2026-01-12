setwd(".")
source("./functions/Libraries.R", local = TRUE)

tmap::ttm()

ui <- fluidPage(
  titlePanel("BALANCE DE MASA GLACIAR"),
  
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
      downloadButton("downloadExcel", "Descargar resultados (.xlsx)")
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
        tabPanel("Balance de masa promedio",
                 br(),
                 leafletOutput("contour_map", height = "650px", width = "100%")
        ),
        tabPanel("Balance anual",
                 br(),
                 leafletOutput("balance_anual", height = "650px", width = "100%")
        ),
        tabPanel("Altura vs. Balance",
                 br(),
                 plotlyOutput("balance_altura", height = "650px")
        ),
        tabPanel("Altura vs. Balance por año",
                 br(),
                 plotlyOutput("balance_altura_anual", height = "650px")
        ),
        tabPanel(
          "Lliboutry No Ajustado",
          
          h4("Variabilidad vs Altitud"),
          plotlyOutput("sd_alt"),
          
          h4("βₜ anual"),
          plotlyOutput("beta_t"),
          
          h4("Perfil MB por año"),
          selectInput("year_mb", "Seleccionar año:", choices = NULL, width = "220px"),
          plotlyOutput("mb_profile_plot", height = "300px"),
          
          h4("Distribución de residuales"),
          plotlyOutput("residuals")
        )
        # ← Aquí puedes añadir más pestañas en el futuro si quieres (mapa, tablas, etc.)
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
  graficos_noAjustado <- reactiveVal(NULL)
  output$contour_map <- renderLeaflet({
    req(mapa_reactivo())     
    tmap::tmap_leaflet(mapa_reactivo())        
  })
  
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
  
  output$sd_alt <- renderPlotly({
    req(graficos_noAjustado())
    graficos_noAjustado()$sd_alt
  })
  
  output$beta_t <- renderPlotly({
    req(graficos_noAjustado())
    graficos_noAjustado()$beta_t
  })
  
  output$residuals <- renderPlotly({
    req(graficos_noAjustado())
    graficos_noAjustado()$residuals
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
    wb <<- openxlsx::createWorkbook()
    
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
        m = input$range_a
      )
      
      resultado <- do.call(Geoglacier, args_list)
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
      if (!is.null(resultado$Graph_NonAdjPlots) && 
          is.list(resultado$Graph_NonAdjPlots)) {
        
        graficos_noAjustado(resultado$Graph_NonAdjPlots)   # ← GUARDAR
        
        add_log("Grafico de Balance No ajustado recibido correctamente", "#28a745")
        showNotification("Gráfico de balance no ajustado Generado", type = "default", duration = 5)
        
      } else {
        graficos_noAjustado(NULL)
        add_log("ADVERTENCIA: Grafico Balance No Ajustado Nulo", "#ffc107")
      }
      add_log("Modelo ejecutado correctamente", "#28a745")
      add_log("EJECUCIÓN FINALIZADA", "#28a745")
      
      showNotification("¡Modelo ejecutado con éxito!", type = "default", duration = 8)
      
    }, error = function(e) {
      add_log(paste("ERROR:", e$message), "#dc3545")
      showNotification("Error en la ejecución del modelo", type = "error", duration = NULL)
    })
  })
  
  # --------------------------------------------------
  # Selector dinámico para perfiles MB por año
  # --------------------------------------------------
  observe({
    req(graficos_noAjustado())
    req(graficos_noAjustado()$mb_years)
    
    updateSelectInput(
      session,
      "year_mb",
      choices = graficos_noAjustado()$mb_years,
      selected = max(graficos_noAjustado()$mb_years)
    )
  })
  
  # --------------------------------------------------
  # Render del perfil MB seleccionado (Plotly)
  # --------------------------------------------------
  output$mb_profile_plot <- renderPlotly({
    req(input$year_mb, graficos_noAjustado()$mb_profiles_list)
    
    p_ggplot <- graficos_noAjustado()$mb_profiles_list[[input$year_mb]]
    
    ggplotly(p_ggplot, tooltip = c("x", "y")) %>%
      layout(
        hovermode = "closest",
        legend = list(orientation = "h", xanchor = "center", x = 0.5, y = -0.1)
      )
  })
  
  # === DESCARGA EXCEL ===
  output$downloadExcel <- downloadHandler(
    filename = "Resultados_Lliboutry.xlsx",
    content = function(file) {
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)