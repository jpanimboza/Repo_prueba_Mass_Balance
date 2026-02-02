Balance_Altura_Plots <- function(ai_df, brick_list, dem) {
  
  # =========================================
  # 1. Gráfico: Balance vs Altitud (general)
  # =========================================
  ai_df_z <- ai_df
  ai_df_z$z <- terra::extract(rast(dem), ai_df_z[c("x", "y")], method = "bilinear")[, 2]
  ai_df_z <- ai_df_z[!is.na(ai_df_z$value), ]
  
  # Calcular tendencia con loess
  ai_df_z <- ai_df_z %>%
    arrange(z) %>%
    mutate(tendencia = predict(loess(value ~ z, span = 0.4)))
  
  balance_altura <- plot_ly(
    data = ai_df_z, 
    x = ~z, 
    y = ~value,
    type = 'scatter',
    mode = 'markers',
    marker = list(size = 6, opacity = 0.7, color = ~z, colorscale = 'Viridis'),
    text = ~paste("Altitud:", round(z, 1), "m<br>Balance:", round(value, 3)),
    hoverinfo = "text",
    name = "Datos"
  ) %>%
    add_trace(
      x = ~z, 
      y = ~tendencia,
      line = list(color = 'grey', width = 2, opacity = 0.5),
      type = "scatter",
      mode = "lines",
      name = "Tendencia",
      inherit = FALSE
    ) %>%
    layout(
      title = "Relación Balance ∼ Altitud",
      xaxis = list(title = "Altitud (m)"),
      yaxis = list(title = "Balance (m w.e./año)"),
      hovermode = "closest"
    )
  
  # =========================================
  # 2. Gráfico: Balance anual vs Altitud
  # =========================================
  yearly_b_df <- as.data.frame(brick_list$LMf, xy = TRUE, na.rm = TRUE)
  yearly_b_df$z <- terra::extract(rast(dem), yearly_b_df[c("x", "y")], method = "bilinear")[, 2]
  
  # Preparar datos por separado
  yearly_data <- yearly_b_df %>%
    select(z, starts_with("X20")) %>%
    pivot_longer(
      cols = starts_with("X20"),
      names_to = "fecha",
      values_to = "balance"
    ) %>%
    mutate(año = substr(fecha, 2, 5)) %>%
    drop_na(z, balance) %>%
    group_by(año) %>%
    arrange(z) %>%
    mutate(tendencia = predict(loess(balance ~ z, span = 0.4))) %>%   
    ungroup()
  
  # Crear gráfico vacío primero
  balance_altura_anual <- plot_ly(data = yearly_data, colors = "Paired")
  
  # Agregar puntos
  balance_altura_anual <- balance_altura_anual %>%
    add_trace(
      x = ~z, 
      y = ~balance, 
      color = ~año,
      type = "scatter",
      mode = "markers",
      marker = list(size = 5, opacity = 0.7),
      text = ~paste(
        "Año:", año,
        "<br>Altitud:", round(z), "m",
        "<br>Balance:", round(balance, 3)
      ),
      hoverinfo = "text",
      name = ~año,
      legendgroup = ~año
    )
  
  # Agregar líneas SIN heredar marker
  balance_altura_anual <- balance_altura_anual %>%
    add_trace(
      x = ~z, 
      y = ~tendencia,
      color = ~año,
      type = "scatter",
      mode = "lines",
      line = list(width = 5),
      name = ~paste(año, "tendencia"),
      showlegend = TRUE,
      legendgroup = ~año,
      inherit = FALSE  # ← Clave: no heredar propiedades
    )
  
  # Layout
  balance_altura_anual <- balance_altura_anual %>%
    layout(
      title = "Balance anual vs Altitud",
      xaxis = list(title = "Altitud (m)"),
      yaxis = list(title = "Balance (m w.e./año)"),
      hovermode = "closest",
      legend = list(title = list(text = "Año"))
    )
  
  # =========================================
  # 3. Retornar ambos gráficos
  # =========================================
  return(list(
    balance_altura = balance_altura,
    balance_altura_anual = balance_altura_anual
  ))
}