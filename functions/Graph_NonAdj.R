Graph_NonAdj <- function(att, data4model, bo_t, Nyear, bit_LM, rawData, B_LMf_w, B_LMf_idw, B_LMf_krg, e_obs_mod){
  
  plots <- list()
  
  # =====================================================
  # 1) SD vs ALTITUD
  # =====================================================
  j <- .2
  k <- 50
  
  xl <- c(min(data4model$sd4model)-j-min(data4model$sd4model)%%j,
          max(data4model$sd4model)+j-max(data4model$sd4model)%%j)
  
  yl <- c(min(data4model$z)-k-min(data4model$z)%%k,
          max(data4model$z)+k-max(data4model$z)%%k)
  
  p <- ggplot(data4model, aes(x = z, y = sd4model)) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2),
                se = FALSE, colour = "red", linewidth = 1) +
    geom_point(colour = "black") +
    scale_x_continuous(limits = yl, breaks = seq(yl[1], yl[2], k)) +
    scale_y_continuous(limits = xl, breaks = seq(xl[1], xl[2], j)) +
    xlab("Altitude (m)") +
    ylab(expression("Standard deviation (m w.e.a"^-1*")")) +
    theme_test()
  
  plots$sd_altitude <- ggplotly(p)
  
  # =====================================================
  # 2) βt ANUAL
  # =====================================================
  bts <- data.frame(bo_t)
  a_comp <- data.frame(Nyear = seq(min(year(Nyear)), max(year(Nyear)), 1))
  bts <- left_join(a_comp, cbind(Nyear = year(Nyear), bts))
  
  k <- 2
  ybr <- c(min(bts$bo_t, na.rm = TRUE)-k-min(bts$bo_t, na.rm = TRUE)%%k,
           max(bts$bo_t, na.rm = TRUE)+k-max(bts$bo_t, na.rm = TRUE)%%k)
  
  p <- ggplot(bts, aes(x = Nyear, y = bo_t)) +
    geom_line(colour = "red", linewidth = 0.75) +
    geom_point(colour = "red") +
    scale_y_continuous(limits = ybr, breaks = seq(ybr[1], ybr[2], k)) +
    xlab("Year") +
    ylab(expression(paste(beta[t]," (m w.e.a"^-1*")"))) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  
  plots$beta_t <- ggplotly(p)
  
  # =====================================================
  # 3) PERFIL MB vs ALTITUD
  # =====================================================
  gg_LM_long <- bit_LM[c("z", Nyear)] %>%
    pivot_longer(all_of(Nyear), names_to = "year", values_to = "value", values_drop_na = TRUE) %>%
    mutate(year = as.character(year)) %>%
    arrange(z)
  
  gg_RW_long <- rawData[c("z", Nyear)] %>%
    pivot_longer(all_of(Nyear), names_to = "year", values_to = "value", values_drop_na = TRUE) %>%
    mutate(year = as.character(year)) %>%
    arrange(z)
  
  years_mb <- sort(unique(gg_LM_long$year))
  
  mb_plots <- list()
  for (yr in years_mb) {
    p_yr <- ggplot() +
      geom_hline(yintercept = 0, colour = "grey") +
      geom_path(data = filter(gg_LM_long, year == yr),
                aes(x = z, y = value), colour = "red", linewidth = 1) +
      geom_point(data = filter(gg_RW_long, year == yr),
                 aes(x = z, y = value), colour = "black", size = 1) +
      labs(x = "Altitud (m)",
           y = "Mass Balance (m w.e.a⁻¹)",
           title = paste("Perfil MB - Año", yr)) +
      theme_test() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5))
    
    mb_plots[[yr]] <- p_yr   # guardamos ggplot base
  }
  
  plots$mb_profiles_list <- mb_plots
  plots$mb_years <- years_mb   # para usar en selectInput
  
  # =====================================================
  # 4) ARRAYS (SIN CAMBIOS LÓGICOS)
  # =====================================================
  array <- list()
  
  array$B_a <- data.frame(
    B_LMf_w_a   = B_LMf_w$B_LMf_w_a,
    B_LMf_idw_a = B_LMf_idw$B_LMf_idw_a,
    B_LMf_krg_a = B_LMf_krg$B_LMf_krg_a
  )
  rownames(array$B_a) <- Nyear
  
  att <- data.frame(matrix(0, 1, ncol(array$B_a)))
  rownames(att) <- ymd(Nyear[1]) - years(1)
  colnames(att) <- colnames(array$B_a)
  
  array$B_a <- rbind(att, array$B_a)
  
  array$B_tot <- data.frame(
    B_LMf_w_tot   = B_LMf_w$B_LMf_w_tot,
    B_LMf_idw_tot = B_LMf_idw$B_LMf_idw_tot,
    B_LMf_krg_tot = B_LMf_krg$B_LMf_krg_tot
  )
  rownames(array$B_tot) <- Nyear
  
  att <- data.frame(matrix(0, 1, ncol(array$B_tot)))
  rownames(att) <- ymd(Nyear[1]) - years(1)
  colnames(att) <- colnames(array$B_tot)
  
  array$B_tot <- rbind(att, array$B_tot)
  
  # =====================================================
  # 5) SERIES TEMPORALES (ANUAL Y ACUMULADO)
  # =====================================================
  array_long <- lapply(array, function(x){
    tibble::rownames_to_column(x, "Nyear") |>
      pivot_longer(-Nyear, names_to = "data", values_to = "value")
  })
  
  p <- ggplot(array_long$B_a, aes(x = Nyear, y = value, color = data, group = data)) +
    geom_line(linewidth = 1) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  
  plots$annual_balance <- ggplotly(p)
  
  p <- ggplot(array_long$B_tot, aes(x = Nyear, y = value, color = data, group = data)) +
    geom_line(linewidth = 1) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  
  plots$cumulative_balance <- ggplotly(p)
  
  # =====================================================
  # 6) RESIDUALES
  # =====================================================
  array$e_obs_mod <- e_obs_mod[Nyear] |>
    cbind(z = e_obs_mod$z) |>
    pivot_longer(1:length(Nyear), names_to = "year", values_to = "value") |>
    na.omit()
  
  p <- ggplot(array$e_obs_mod, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 9,
                   fill = "grey", colour = "black") +
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = sd(array$e_obs_mod$value)),
                  colour = "red", linewidth = 0.8) +
    theme_test()
  
  plots$residuals <- ggplotly(p)
  ## =====================================================
  ## SALIDA
  ## =====================================================
  return(list(
    plots  = plots,
    arrays = array,
    bts    = bts
  ))
}
