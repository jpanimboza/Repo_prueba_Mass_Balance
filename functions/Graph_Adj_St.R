Graph_Adj <- function(att, data4model, bo_t, Nyear, bit_LM, rawData, 
                      bit_LMf_w_adj, bit_LMf_idw_adj, bit_LMf_krg_adj,
                      B_LMf_w, B_LMf_idw, B_LMf_krg,
                      e_obs_mod, temp_dir, wb){
  dir_download <- file.path(temp_dir, "Resultados")
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  ## Graph of the altitudinal gradient of the variability of point mass balance (standard deviation) 
  # Graph parameters
  j <- .2
  k <- 50
  xl <- c(min(data4model$sd4model)-j-min(data4model$sd4model)%%j, 
          max(data4model$sd4model)+j-max(data4model$sd4model)%%j)
  yl <- c(min(data4model$z)-k-min(data4model$z)%%k, 
          max(data4model$z)+k-max(data4model$z)%%k)
  
  # Graph creation
  plot <- ggplot(data4model, aes(x = z, y = sd4model)) +
    geom_smooth(colour = "red", method = "lm",
                formula = y ~ poly(x, 2), se =F, lwd = 1) +
    geom_point(colour = "black") +
    scale_x_continuous(limits = c(yl[1],yl[2]), breaks = seq(yl[1],yl[2], k)) +
    scale_y_continuous(limits = c(xl[1],xl[2]), breaks = seq(xl[1],xl[2], j)) +
    xlab("Altitude (m)") +
    ylab(expression("Standard deviation (m w.e.a"^-1*")")) +
    theme_test() + theme(axis.text.x = element_text(colour = "black",
                                                    angle = 90, vjust = 0.5),
                         axis.text.y = element_text(colour = "black"))
  # Save the plot as a JPEG
  ggsave(file.path(dir_download, "out_StdDev_Elev.jpg"), plot = plot, dpi = 2000)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Graph of the deviation of the average mass balance (beta t) 
  # Graph parameters
  bts <- data.frame(bo_t)
  a_comp <- data.frame(Nyear = seq(min(year(Nyear)), max(year(Nyear)), 1))
  bts <- left_join(a_comp, cbind(Nyear = year(Nyear), bts))
  
  k <- 2
  ybr <- c(min(bts$bo_t, na.rm = T)-k-min(bts$bo_t, na.rm = T)%%k, 
           max(bts$bo_t, na.rm = T)+k-max(bts$bo_t, na.rm = T)%%k)
  
  # Create plot
  plot <- ggplot(bts, aes(x = Nyear, y = bo_t, group = 1)) +
    geom_line(colour = "red", linewidth = 0.75) +
    geom_point(colour = "red") +
    scale_x_continuous(breaks = bts$Nyear) +
    scale_y_continuous(limits = ybr, breaks = seq(ybr[1], ybr[2], k)) +
    xlab("Year") +
    ylab(expression(paste(beta[t]," (m w.e.a"^-1*")"))) +
    theme_test() +
    theme(axis.text.x = element_text(colour = "black",
                                     angle = 90, vjust = 0.5),
          axis.text.y = element_text(colour = "black"))
  # Save the plot as a JPEG
  ggsave(file.path(dir_download, "out_Betat_anual.jpg"), plot = plot, dpi = 2000)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
  # Graph of the mass balance profile vs. elevation
  # Graph parameters
  gg_LM <- bit_LM[c("z",Nyear)] %>% 
    pivot_longer(cols = all_of(Nyear), names_to = "year", 
                 values_to = "value", values_drop_na = TRUE)
  gg_RW <- rawData[c("z",Nyear)] %>% 
    pivot_longer(cols = all_of(Nyear), names_to = "year", 
                 values_to = "value", values_drop_na = TRUE)
  gg_LM <- gg_LM[order(gg_LM$z),]
  gg_RW <- gg_RW[order(gg_RW$z),]
  
  j <- 200
  k <- 2
  xl <- c(min(gg_LM$z, na.rm = T) - min(gg_LM$z, na.rm = T)%%j, 
          max(gg_LM$z, na.rm = T)+j- max(gg_LM$z, na.rm = T)%%j)
  yl <- c(min(gg_LM$value, na.rm = T) - min(gg_LM$value, na.rm = T)%%k, 
          max(gg_LM$value, na.rm = T)+k- max(gg_LM$value, na.rm = T)+1.5%%k)
  
  # Create plot
  plot <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_path(data = gg_LM, aes(x = z, y = value, colour = "3")) +
    geom_point(data = gg_RW, aes(x = z, y = value, colour = "1"), size = 0.75) +
    scale_y_continuous(limits = yl, breaks = seq(yl[1], yl[2], k)) +
    scale_color_manual(name = "Datos", values = c("black", "red"),
                       labels = c("MB", "MB Lliboutry Model"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("blank","solid"),
                         shape = c(16, NA)))) +
    facet_wrap(~year, scales = "free_y") +
    xlab("Altitud (m)") +
    ylab(expression("Mass Balance (m w.e.a"^-1*")")) +
    theme_test() + theme(axis.text.x = element_text(colour = "black",
                                                    size = 5),
                         axis.text.y = element_text(colour = "black",
                                                    size = 5),
                         panel.spacing = unit(.65, "cm"), legend.position = "none",
                         axis.title = element_text(size = 10))
  
  ggsave(file.path(dir_download, "out_Bgrad_a.jpg"), plot = plot, dpi = 2000)
  
  
  gg_LM <- bit_LMf_w_adj[c("z",Nyear)] %>% 
    pivot_longer(cols = all_of(Nyear), names_to = "year", 
                 values_to = "value", values_drop_na = TRUE)
  gg_LM <- gg_LM[order(gg_LM$z),]
  plot <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_path(data = gg_LM, aes(x = z, y = value, colour = "3")) +
    scale_y_continuous(limits = yl, breaks = seq(yl[1], yl[2], k)) +
    scale_color_manual(name = "Datos", values = c("blue"),
                       labels = c("MB", "MB Lliboutry Model"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid"),
                         shape = c(16)))) +
    facet_wrap(~year, scales = "free_y") +
    xlab("Altitud (m)") +
    ylab(expression("Adjusted Mass Balance (m w.e.a"^-1*")")) +
    theme_test() + theme(axis.text.x = element_text(colour = "black",
                                                    size = 5),
                         axis.text.y = element_text(colour = "black",
                                                    size = 5),
                         panel.spacing = unit(.65, "cm"), legend.position = "none",
                         axis.title = element_text(size = 10))
  ggsave(file.path(dir_download, "out_Bgrad_w_adj.jpg"), plot = plot, dpi = 2000)
  
  
  gg_LM <- bit_LMf_idw_adj[c("z",Nyear)] %>% 
    pivot_longer(cols = all_of(Nyear), names_to = "year", 
                 values_to = "value", values_drop_na = TRUE)
  gg_LM <- gg_LM[order(gg_LM$z),]
  plot <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_path(data = gg_LM, aes(x = z, y = value, colour = "3")) +
    scale_y_continuous(limits = yl, breaks = seq(yl[1], yl[2], k)) +
    scale_color_manual(name = "Datos", values = c("blue"),
                       labels = c("MB", "MB Lliboutry Model"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid"),
                         shape = c(16)))) +
    facet_wrap(~year, scales = "free_y") +
    xlab("Altitud (m)") +
    ylab(expression("Adjusted Mass Balance (m w.e.a"^-1*")")) +
    theme_test() + theme(axis.text.x = element_text(colour = "black",
                                                    size = 5),
                         axis.text.y = element_text(colour = "black",
                                                    size = 5),
                         panel.spacing = unit(.65, "cm"), legend.position = "none",
                         axis.title = element_text(size = 10))
  
  ggsave(file.path(dir_download, "out_Bgrad_idw_adj.jpg"), plot = plot, dpi = 2000)
  
  
  gg_LM <- bit_LMf_krg_adj[c("z",Nyear)] %>% 
    pivot_longer(cols = all_of(Nyear), names_to = "year", 
                 values_to = "value", values_drop_na = TRUE)
  gg_LM <- gg_LM[order(gg_LM$z),]
  plot <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_path(data = gg_LM, aes(x = z, y = value, colour = "3")) +
    scale_y_continuous(limits = yl, breaks = seq(yl[1], yl[2], k)) +
    scale_color_manual(name = "Datos", values = c("blue"),
                       labels = c("MB", "MB Lliboutry Model"),
                       guide = guide_legend(override.aes = list(
                         linetype = c("solid"),
                         shape = c(16)))) +
    facet_wrap(~year, scales = "free_y") +
    xlab("Altitud (m)") +
    ylab(expression("Adjusted Mass Balance (m w.e.a"^-1*")")) +
    theme_test() + theme(axis.text.x = element_text(colour = "black",
                                                    size = 5),
                         axis.text.y = element_text(colour = "black",
                                                    size = 5),
                         panel.spacing = unit(.65, "cm"), legend.position = "none",
                         axis.title = element_text(size = 10))
  
  ggsave(file.path(dir_download, "out_Bgrad_krg_adj.jpg"), plot = plot, dpi = 2000)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  array <- list() #-> Empty list.
  # Annual Data Array
  array$B_a <- cbind(B_LMf_w_a = B_LMf_w$B_LMf_w_a,
                     B_LMf_idw_a = B_LMf_idw$B_LMf_idw_a,
                     B_LMf_krg_a = B_LMf_krg$B_LMf_krg_a,
                     B_LMfadj_an_w_a = B_LMf_w$B_LMfadj_w_a,
                     B_LMfadj_an_idw_a = B_LMf_idw$B_LMfadj_idw_a, 
                     B_LMfadj_an_krg_a = B_LMf_krg$B_LMfadj_krg_a, 
                     B_LMfadj_sp_w_a = B_LMf_w$B_LMf_w_adj_rg_a, 
                     B_LMfadj_sp_idw_a = B_LMf_idw$B_LMf_idw_adj_sp_a,
                     B_LMfadj_sp_krg_a = B_LMf_krg$B_LMf_krg_adj_sp_a)
  
  row.names(array$B_a) <- Nyear
  att <- data.frame(matrix(0, nrow=1, ncol=ncol(array$B_a)))
  rownames(att) <- c(ymd(Nyear[1]) - years(1))
  colnames(att) <- colnames(array$B_a)
  array$B_a <- rbind(att, array$B_a)
  
  # Accumulated Data Array
  array$B_tot <- cbind(B_LMf_w_tot = B_LMf_w$B_LMf_w_tot,
                       B_LMf_idw_tot = B_LMf_idw$B_LMf_idw_tot,
                       B_LMf_krg_tot = B_LMf_krg$B_LMf_krg_tot,
                       B_LMfadj_an_w_tot = B_LMf_w$B_LMfadj_w_tot, 
                       B_LMfadj_an_idw_tot = B_LMf_idw$B_LMfadj_idw_tot,
                       B_LMfadj_an_krg_tot = B_LMf_krg$B_LMfadj_krg_tot, 
                       B_LMfadj_sp_w_tot = B_LMf_w$B_LMf_w_adj_rg_tot, 
                       B_LMfadj_sp_idw_tot = B_LMf_idw$B_LMf_idw_adj_sp_tot,
                       B_LMfadj_sp_krg_tot = B_LMf_krg$B_LMf_krg_adj_sp_tot)
  
  row.names(array$B_tot) <- Nyear
  att <- data.frame(matrix(0, nrow=1, ncol=ncol(array$B_tot)))
  rownames(att) <- c(ymd(Nyear[1]) - years(1))
  colnames(att) <- colnames(array$B_tot)
  array$B_tot <- rbind(att, array$B_tot)
  
  addWorksheet(wb, "out_B_a")
  writeData(wb, sheet = "out_B_a", array$B_a)
  
  addWorksheet(wb, "out_B_tot")
  writeData(wb, sheet = "out_B_tot", array$B_tot)
  
  # Array of data to create graphs
  array <- lapply(setNames(1:2, names(array)), function(i){
    array[[i]] %>% tibble::rownames_to_column("Nyear")  %>%
      pivot_longer(2:10, names_to = "data")
  }) 
  
  # Annual Mass Balance Graph 
  k <- 0.5
  ybr <- c(min(array$B_a$value)-min(array$B_a$value)%%k, 
           max(array$B_a$value)+k-max(array$B_a$value)%%k)
  
  # Create plot
  plot <- ggplot(array$B_a[array$B_a$data %in% c("B_LMf_w_a",
                                                 "B_LMf_idw_a",
                                                 "B_LMf_krg_a",
                                                 "B_LMfadj_an_w_a",
                                                 "B_LMfadj_an_idw_a",
                                                 "B_LMfadj_an_krg_a",
                                                 "B_LMfadj_sp_w_a",
                                                 "B_LMfadj_sp_idw_a",
                                                 "B_LMfadj_sp_krg_a"),],
                 aes(x = Nyear, y = value, group = data, color = data)) +
    geom_line(linewidth=1) +
    scale_x_discrete(expand = c(0, 0.2), labels = year(c(ymd(Nyear[1]) - years(1),Nyear))) +
    scale_y_continuous(limits = ybr, breaks = seq(ybr[1],ybr[2],.5),
                       guide = waiver()) +#guide_prism_minor()
    scale_color_manual(breaks=c("B_LMf_w_a",
                                "B_LMf_idw_a",
                                "B_LMf_krg_a",
                                "B_LMfadj_an_w_a",
                                "B_LMfadj_an_idw_a",
                                "B_LMfadj_an_krg_a",
                                "B_LMfadj_sp_w_a",
                                "B_LMfadj_sp_idw_a",
                                "B_LMfadj_sp_krg_a"),
                       labels = c("B_LMf_w_a",
                                  "B_LMf_idw_a",
                                  "B_LMf_krg_a",
                                  "B_LMfadj_an_w_a",
                                  "B_LMfadj_an_idw_a",
                                  "B_LMfadj_an_krg_a",
                                  "B_LMfadj_sp_w_a",
                                  "B_LMfadj_sp_idw_a",
                                  "B_LMfadj_sp_krg_a"),
                       values = c("brown4", "red", "orange", "blue4", "dodgerblue", "cyan", "darkgreen", "green3", "green")) +
    xlab("year") + ylab(expression("Mass balance (m w.e.a"^-1*")")) +
    labs(color = "Datos") +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  # Save the plot as a JPEG
  ggsave(file.path(dir_download, "out_B_a.jpg"), plot = plot, dpi = 2000)
  
  
  # Cumulative Mass Balance Graph
  # Create plot
  plot <- ggplot(array$B_tot[array$B_tot$data %in% c("B_LMf_w_tot",
                                                     "B_LMf_idw_tot",
                                                     "B_LMf_krg_tot",
                                                     "B_LMfadj_an_w_tot",
                                                     "B_LMfadj_an_idw_tot",
                                                     "B_LMfadj_an_krg_tot",
                                                     "B_LMfadj_sp_w_tot",
                                                     "B_LMfadj_sp_idw_tot",
                                                     "B_LMfadj_sp_krg_tot"),],
                 aes(x = Nyear, y = value, group = data, color = data)) +
    geom_line(linewidth=1) +
    scale_x_discrete(expand = c(0, 0.2), labels = year(c(ymd(Nyear[1]) - years(1),Nyear))) +
    scale_y_continuous(guide = waiver()) +#guide_prism_minor()
    scale_color_manual(breaks= c("B_LMf_w_tot",
                                 "B_LMf_idw_tot",
                                 "B_LMf_krg_tot",
                                 "B_LMfadj_an_w_tot",
                                 "B_LMfadj_an_idw_tot",
                                 "B_LMfadj_an_krg_tot",
                                 "B_LMfadj_sp_w_tot",
                                 "B_LMfadj_sp_idw_tot",
                                 "B_LMfadj_sp_krg_tot"),
                       labels = c("B_LMf_w_tot",
                                  "B_LMf_idw_tot",
                                  "B_LMf_krg_tot",
                                  "B_LMfadj_an_w_tot",
                                  "B_LMfadj_an_idw_tot",
                                  "B_LMfadj_an_krg_tot",
                                  "B_LMfadj_sp_w_tot",
                                  "B_LMfadj_sp_idw_tot",
                                  "B_LMfadj_sp_krg_tot"),
                       values = c("brown4", "red", "orange", "blue4", "dodgerblue", "cyan", "darkgreen", "green3", "green")) +
    xlab("Year") + ylab(expression("Cumulative mass balance (m w.e.a"^-1*")")) +
    labs(color = "Datos") +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
  # Save the plot as a JPEG
  ggsave(file.path(dir_download, "out_B_tot.jpg"), plot = plot, dpi = 2000)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
  # Graphs of Height vs. Residuals
  # Data arrangement
  array$e_obs_mod <- e_obs_mod[Nyear] %>% cbind(z = e_obs_mod$z) %>%
    pivot_longer(1:(length(Nyear)), names_to = 'year', values_to = 'value') %>%
    na.omit()
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
  # Graph of residual distribution
  # Create plot
  k <- 1
  xbr <- c(0, max(density(array$e_obs_mod$value)$y)+k-max(density(array$e_obs_mod$value)$y)%%k)
  ybr <- c(min(array$e_obs_mod$value)-min(array$e_obs_mod$value)%%k, 
           max(array$e_obs_mod$value)+k-max(array$e_obs_mod$value)%%k)
  plot <- ggplot(array$e_obs_mod, aes(x = value)) +
    geom_histogram(aes(y = ..density..), colour = 1, fill = "gray",
                   bins = 15) +
    stat_function(fun = dnorm, args = list(mean = 0,
                                           sd = sd(array$e_obs_mod$value)),
                  col = "red", size = 0.75) +
    scale_x_continuous(expand = c(0,.5), limits = ybr) +
    scale_y_continuous(limits = xbr,
                       breaks = seq(0,25,0.25)) +
    xlab(expression("Residuals (m w.e.a"^-1*")")) + ylab("Density") +
    ggtitle("Residuals distribution: Observations vs. Lliboutry's model") +
    theme_test() + theme(axis.text.x = element_text(colour = "black"),
                         axis.text.y = element_text(colour = "black"))
  # Save the plot as a JPEG
  ggsave(file.path(dir_download, "out_e_obs_mod.jpeg"), plot = plot, dpi = 2000)
  rm(plot)
  
  return(list(
    arrays <- array,
    btss <- bts
  ))
}