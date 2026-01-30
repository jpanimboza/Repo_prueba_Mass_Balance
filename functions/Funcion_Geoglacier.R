Geoglacier <- function(rawData, shapes_path, dem, Resolution, crsID, N.Obs,summit, geodetic_balance, m, wb, temp_dir){
  source("./functions/Lliboutry_NL.r")
  NLMR<-NLM(rawData, shapes_path, dem, Resolution, crsID, N.Obs, summit,wb)
  e_obs_mod <- NLMR$e_obs_mod
  Nyear <- NLMR$Nyear
  data4model <- NLMR$data4model
  bit_LM <- NLMR$bit_LM
  bit_LMf <- NLMR$bit_LMf
  bo_t <- NLMR$bo_t
  beta_it <- NLMR$beta_it
  brick_raw <-NLMR$brick_raw
  rawData <- NLMR$rawData
  shap <- NLMR$shap
  shape_list <- NLMR$shape_list
  ai_df <- NLMR$ai_df
  XY_coords <- NLMR$XY_coords
  e_obs_mod_sd.tot <- sd(unlist(e_obs_mod[-c(1:3)]), na.rm=T)
  
  #To check the normality of the data set by means of the Shapiro-Wilk Test
  shapiro_result <- shapiro.test(unlist(e_obs_mod[Nyear]))
  # Creation of Brickfile of model outputs for Mass Balance calculation
  brick_list <- list(
    LM = bit_LM, 
    LMf = bit_LMf, 
    beta_it = beta_it, 
    e_obs_mod = e_obs_mod) # Guarda data en lista para su conversion a raster
  source("./functions/BrickBuilder.R")
  brick_list <- BrickBuilder(brick_list,brick_raw,Nyear)
  brick_list$Raw <- brick_raw
  brick_list <- brick_list[c("Raw", names(brick_list)[-5])] # Ordena la lista (reemplace -4 si se crea el beta_it despues)
  n.stk.abl <- sum(rawData$nf[which(rawData$ai<0)]/length(Nyear))
  n.stk.acc <- sum(rawData$nf[which(rawData$ai>0)]/length(Nyear))
  #######################################################################################
  #===== CALCULATION OF THE DISTRIBUTED MASS BALANCE FOR DATA WITHOUT GEODETIC ADJUSTMENT =======
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Weighted area/altitude range method for Lliboutry model output data
  source("./functions/B_LMf_W.R")
  aux = 0
  WAreaR = WArea(bit_LMf, m, shapes_path, dem, crsID, shap, Nyear,aux)
  bit_LMf_w_a <- WAreaR$bit_LMf_w_a
  B_LMf_w <- WAreaR$B_LMf_w
  stat.w.sd <- WAreaR$stat.w.sd
  # Export the DataFrame obtained through the Weighted Area Method.
  bit_LMf_w_a <-WAreaR$bit_LMf_w_a
  addWorksheet(wb, "out_bit_LMf_w_a_rg")
  writeData(wb, sheet = "out_bit_LMf_w_a_rg", bit_LMf_w_a)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2. Inverse Distance Weighted interpolation method
  aux = 0
  source("./functions/B_LMf_Idw.R")
  IdwGBR = IdwGB(bit_LMf, dem, crsID, brick_raw, Nyear, aux)
  bit_LMf_idw <- IdwGBR$bit_LMf_idw
  stat.idw.CV_LMf <- IdwGBR$stat.idw.CV_LMf
  # Creation and masking of RasterBrick objects according to the annual area (Is necessary for plotting raster cut by area)
  att <- BrickBuilder(list(LMf_idw = bit_LMf_idw), brick_raw, Nyear)
  brick_list$LMf_idw <- lapply(setNames(1:length(Nyear), Nyear), function(i){
    raster::mask(att$LMf_idw[[i]], as(st_zm(st_as_sf(shape_list[[shap[i]]])), "Spatial")) # remplaszo shape_list[[shap[i]]]
  }) %>% brick()
  
  # Calculation of the annual net mass balance for the distribution using the IDW method
  B_LMf_idw <-data.frame(area_glc = sapply(1:length(Nyear), function(i){
    # Calculate area using sf for sf objects
    sum(st_area(shape_list[[shap[i]]]))
  }), 
  LMf_idw = sapply(1:length(Nyear), function(i){
    cellStats(brick_list$LMf_idw[[i]], stat='sum', na.rm=TRUE, asSample=TRUE)
  })
  )
  
  B_LMf_idw$B_LMf_idw_a <- B_LMf_idw$LMf_idw*(Resolution^2)/B_LMf_idw$area_glc
  B_LMf_idw$B_LMf_idw_tot <- cumsum(B_LMf_idw$B_LMf_idw_a)
  B_LMf_idw <- cbind(Nyear, B_LMf_idw)
  
  # Export the DataFrame obtained through the IDW Method.
  addWorksheet(wb, "out_bit_LMf_idw")
  writeData(wb, sheet = "out_bit_LMf_idw", bit_LMf_idw)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3. Kriging method for Lliboutry model output data
  aux = 0
  source("./functions/B_LMf_Krg.R")
  KrgGBR <- KrgGB(bit_LMf, dem, crsID, brick_raw, Nyear, aux)
  varplot1_LMf <- KrgGBR$varplot1_LMf
  varplot2_LMf <- KrgGBR$varplot2_LMf
  bit_LMf_krg <- KrgGBR$bit_LMf_krg
  stat.krg.CV_LMf <- KrgGBR$stat.krg.CV_LMf
  p <- ggplot(varplot1_LMf, aes(x = dist, y = gamma)) +
    geom_line(data = varplot2_LMf, linewidth = 0.75) +
    geom_point(colour = "red") +
    facet_wrap(.~year, scales = "free") +
    ggtitle("SemiVariograms of Lliboutry fill") +
    theme_test() +
    xlab("Range distance (m)") + ylab(expression("Semivariance (m w.e.a"^-1*")"))
  
  # Creation and masking of RasterBrick objects according to the annual area
  att <- BrickBuilder(list(LMf_krg = bit_LMf_krg), brick_raw, Nyear)
  brick_list$LMf_krg <- lapply(setNames(1:length(Nyear), Nyear), function(i){
    raster::mask(att$LMf_krg[[i]], as(st_zm(st_as_sf(shape_list[[shap[i]]])), "Spatial")) # remplaszo shape_list[[shap[i]]]
  }) %>% brick()
  
  # Calculation of the distributed mass balance by Kriging
  B_LMf_krg <- data.frame(area_glc = sapply(1:length(Nyear), function(i){
    # Calculate area using sf for sf objects
    sum(st_area(shape_list[[shap[i]]]))
  }),
  LMf_krg = sapply(1:length(Nyear), function(i){
    cellStats(brick_list$LMf_krg[[i]], stat='sum', na.rm=TRUE, asSample=TRUE)
  })
  )
  
  B_LMf_krg$B_LMf_krg_a <- B_LMf_krg$LMf_krg*(Resolution^2)/B_LMf_krg$area_glc # Balance krig ponderado a la talla del pixel
  B_LMf_krg$B_LMf_krg_tot <- cumsum(B_LMf_krg$B_LMf_krg_a)
  B_LMf_krg <- cbind(Nyear, B_LMf_krg)
  
  # Export the DataFrame obtained through Lliboutry merhod interpolated by Kriging Method.
  addWorksheet(wb, "out_bit_LMf_krg")
  writeData(wb, sheet = "out_bit_LMf_krg", bit_LMf_krg)
  
  
  ##########################################################################################
  #Error Analysis - Conventional mass balance
  e.stkabl <- 0.1 #m w.e.yr-1 Gerbaux et al., 2005 e.stkabl <- 0.1
  e.densabl <- 0.00 # d_acc = 0.9
  e.pitacc <- 0.4 #m w.e.yr-1 Gerbaux et al., 2005 e.pitacc <- 0.4
  e.densacc <- 0.025 # d_acc = 0.5
  
  # error of point measurement, considering the number of stakes
  e.Abl <- e.stkabl/sqrt(n.stk.abl)
  e.Acc <- ifelse(n.stk.acc==0, 0, e.pitacc/sqrt(n.stk.acc))
  e_mesMB <- sqrt(e.Abl^2+e.Acc^2)
  # error density, considering the number of stakes
  e_dcice <- e.densabl/sqrt(n.stk.abl)
  e_dcfirm <- ifelse(n.stk.acc==0, 0, e.densacc/sqrt(n.stk.acc))
  e_dc <- sqrt(e_dcfirm^2+e_dcice^2)
  
  e_glc.point.w.tot <- sqrt(e_mesMB^2+(sum(B_LMf_w$B_LMf_w_a)^2)*(e_dc^2))*sqrt(length(Nyear))
  e_glc.point.idw.tot <- sqrt(e_mesMB^2+(sum(B_LMf_idw$B_LMf_idw_a)^2)*(e_dc^2))*sqrt(length(Nyear))
  e_glc.point.krg.tot <- sqrt(e_mesMB^2+(sum(B_LMf_krg$B_LMf_krg_a)^2)*(e_dc^2))*sqrt(length(Nyear))
  # Error spatial averaging of point measurements - local representativeness of local measurements
  e_repres.yr <- sd(unlist(bit_LMf[-c(1:3)]), na.rm = T) / sqrt(unlist(lapply(stat.idw.CV_LMf, function(sublist) sublist[[1]]), use.names=F)[1])
  e_repres.tot <- sd(unlist(bit_LMf[-c(1:3)]), na.rm = T) / sqrt(unlist(lapply(stat.idw.CV_LMf, function(sublist) sublist[[1]]), use.names=F)[1]) #e_repres.tot <- e_repres.yr * sqrt(length(Nyear))
  # lliboutry method
  #Calculado en linea 18
  #   e_obs_mod_sd.tot
  # Error interpolation method
  ## Surface-area weighthed
  std.w <- stat.w.sd^2
  n_obs <- unlist(lapply(stat.idw.CV_LMf, function(sublist) sublist[[1]]), use.names=F)
  e_w_sd.tot <- sqrt(sum(std.w)/((n_obs[1]*length(Nyear))-length(Nyear)))
  ## IDW interpolation error - CrossVal. Based on pooled variance calculation from residuals
  std.idw <- unlist(lapply(stat.idw.CV_LMf, function(sublist) sublist[[5]]), use.names=F)^2
  n_obs <- unlist(lapply(stat.idw.CV_LMf, function(sublist) sublist[[1]]), use.names=F)
  e_idw_sd.tot <- sqrt(sum(std.idw)/((n_obs[1]*length(Nyear))-length(Nyear)))
  ## Kgr interpolation error - CrossVal. Based on pooled variance calculation from residuals
  # Residual standard deviation of the kriging prediction (year) and compute the mean error according crossvalidation approach
  std.krg <- unlist(lapply(stat.krg.CV_LMf, function(sublist) sublist[[5]]), use.names=F)^2
  n_obs <- unlist(lapply(stat.krg.CV_LMf, function(sublist) sublist[[1]]), use.names=F)
  e_krg_sd.tot <- sqrt(sum(std.krg)/((n_obs[1]*length(Nyear))-length(Nyear)))
  
  e_glc.spatial.w.tot <- sqrt(e_repres.tot^2+e_w_sd.tot^2)
  e_glc.spatial.idw.tot <- sqrt(e_repres.tot^2+e_idw_sd.tot^2)
  e_glc.spatial.krg.tot <- sqrt(e_repres.tot^2+e_krg_sd.tot^2)
  
  rm(std.idw, std.krg) #std.w, n_obs
  
  # Error surface area changes 
  e_glac.area.w.tot <-  e_glc.point.w.tot * as.numeric((0.5*st_length(st_boundary(shape_list[[1]]))*3)/st_area(shape_list[[1]])) # asume 5 m de error en la determinacion de area
  e_glac.area.idw.tot <- e_glc.point.idw.tot * as.numeric((0.5*st_length(st_boundary(shape_list[[1]]))*3)/st_area(shape_list[[1]])) # asume 5 m de error en la determinacion de area
  e_glac.area.krg.tot <- e_glc.point.krg.tot * as.numeric((0.5*st_length(st_boundary(shape_list[[1]]))*3)/st_area(shape_list[[1]])) # asume 5 m de error en la determinacion de area
  
  
  # Error of Glaciological mass balance 
  E_glc.w.tot <- sqrt(e_glc.point.w.tot^2+e_glc.spatial.w.tot^2+e_glac.area.w.tot^2)
  E_glc.w.yr <-  E_glc.w.tot / sqrt(length(Nyear))
  
  E_glc.idw.tot <- sqrt(e_glc.point.idw.tot^2+e_glc.point.idw.tot^2+e_glac.area.idw.tot^2)
  E_glc.idw.yr <- E_glc.idw.tot / sqrt(length(Nyear))
  
  E_glc.krg.tot <- sqrt(e_glc.point.krg.tot^2+e_glc.point.krg.tot^2+e_glac.area.krg.tot^2)
  E_glc.krg.yr <- E_glc.krg.tot / sqrt(length(Nyear))
  message("Archivo encontrado")
  message(!is.null(geodetic_balance))
  ####################################################################################### 
  # GEODETIC MASS BALANCE (Dussaillant et al., 2019)
  if (!is.null(geodetic_balance)){
    geo_df <- data.frame(read_excel(geodetic_balance, sheet = 1))
    #######################################################################################
    #=== GLACIOLOGICAL MASS BALANCE VALIDATION - GEODETIC ADJUSTMENT (Zemp et al., 2013) ===
    #=== CALCULATION OF THE CALIBRATED DISTRIBUTED MASS BALANCE WITH GEODETIC ADJUSTMENT =====
    # Frecuentemente el periodo de monitoreo cliaciologico no coincide con el monitoreo geodesico
    # en estos casos se asume un rango de tolerancia de 3 anios (antes/despues) para suponer que la tendencia
    # observada en el periodo coincidente se mantiene
    yr_glc_o <- as.numeric(substr(Nyear[1], 1, 4))
    yr_glc_f <- as.numeric(substr(Nyear[length(Nyear)], 1, 4))
    
    ai_adj_geo <- list()
    bit_LM_w_adj <- list()
    bit_LMf_w_adj <- list()
    bit_LM_idw_adj <- list()
    bit_LMf_idw_adj <- list()
    bit_LM_krg_adj <- list()
    bit_LMf_krg_adj <- list()
    B_LMf_w$B_LMfadj_w_a <- B_LMf_w$B_LMf_w_a
    B_LMf_idw$B_LMfadj_idw_a <- B_LMf_idw$B_LMf_idw_a
    B_LMf_krg$B_LMfadj_krg_a <- B_LMf_krg$B_LMf_krg_a
    E_glc_adj_a <- list() 
    for (k in 1:(nrow(geo_df)-1)) {
      # Define the glaciological time series
      yr_geo_o <- geo_df$yr_geo_o[k] # para Antisana 15a
      yr_geo_f <- geo_df$yr_geo_f[k] # para Antisana 15a
      
      # Calculate the limits with tolerance of +/-2 years 
      # This script first expands the reference interval by a tolerance of 2 years in either direction. 
      # Then it checks whether each interval falls within this extended reference interval. 
      # If it does, it further checks if it crosses the start or end boundary of the original reference interval to determine how it could be distributed.
      # Check if the interval is within the limits and return message
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      if (yr_glc_o <= yr_geo_o & yr_glc_f >= yr_geo_f) {
        idx1 <- which(B_LMf_w$Nyear==paste(yr_geo_o, substr(Nyear[1], 5, 10), sep = ""))
        idx2 <- which(B_LMf_w$Nyear==paste(yr_geo_f, substr(Nyear[1], 5, 10), sep = ""))
        # Analisis estadistico de discrepancia para ajuste
        source("./functions/Reanal_glac_geo_Zemp.R")
        Reanal_glac_geoR = Reanal_glac_geo(B_LMf_w, idx1, idx2, e_mesMB, e_dc, geo_df, k, e_repres.tot, std.w, n_obs, shape_list)
        e_Bgl_a <- Reanal_glac_geoR$e_Bgl_a
        E_glc_adj_a[[k]] <- e_Bgl_a
        names(E_glc_adj_a[[k]]) <- c(paste0("period ", k))
        aux <- Reanal_glac_geoR$aux
        if(aux==1){
          # 1. Geodetic calibration (Annual calibration) contempla el ajuste convencional
          # Ajuste de los datos por de Bm calculado por el metodo de area ponderada
          B_LMf_w$aux <- B_LMf_w$B_LMfadj_w_a - mean(B_LMf_w$B_LMfadj_w_a[(idx1+1):idx2], na.rm=T) + (geo_df$B_geo_a[k])
          B_LMf_w$B_LMfadj_w_a[(idx1+1):idx2] <- B_LMf_w$aux[(idx1+1):idx2]
          
          B_LMf_idw$aux <- B_LMf_idw$B_LMfadj_idw_a - mean(B_LMf_idw$B_LMfadj_idw_a[(idx1+1):idx2], na.rm=T) + (geo_df$B_geo_a[k])
          B_LMf_idw$B_LMfadj_idw_a[(idx1+1):idx2] <- B_LMf_idw$aux[(idx1+1):idx2]
          
          B_LMf_krg$aux <- B_LMf_krg$B_LMfadj_krg_a - mean(B_LMf_krg$B_LMfadj_krg_a[(idx1+1):idx2], na.rm=T) + (geo_df$B_geo_a[k])
          B_LMf_krg$B_LMfadj_krg_a[(idx1+1):idx2] <- B_LMf_krg$aux[(idx1+1):idx2]
          
          # 2. Geodetic calibration (pixel by pixel) asumiendo que el bm_geo_a es el mismo para todos los anios se ajusta el ai (balance promedio)
          # hipotesis de densidad 850 kg/m3, (Huss, 2013) - periodo comun
          # LMf salida Lliboutry ajustado con geodesico LMf_adj
          ai_adj <- data.frame(rawData$ai)
          # Calculo de Delta ai ajustado con BM geodesico
          ai_adj$db_i_w_adj_geo <- (mean(B_LMf_w$B_LMf_w_a[(idx1+1):idx2])-mean(B_LMf_w$B_LMfadj_w_a[(idx1+1):idx2]))/rawData$yi # Balance area ponderada
          ai_adj$db_i_idw_adj_geo <- (mean(B_LMf_idw$B_LMf_idw_a[(idx1+1):idx2])-mean(B_LMf_idw$B_LMfadj_idw_a[(idx1+1):idx2]))/rawData$yi # Balance IDW
          ai_adj$db_i_krg_adj_geo <- (mean(B_LMf_krg$B_LMf_krg_a[(idx1+1):idx2])-mean(B_LMf_krg$B_LMfadj_krg_a[(idx1+1):idx2]))/rawData$yi # Balance kriging
          
          ai_adj$ai_w_adj_geo <- (rawData$ao_i - ai_adj$db_i_w_adj_geo)*rawData$yi
          ai_adj$ai_idw_adj_geo <- (rawData$ao_i - ai_adj$db_i_idw_adj_geo)*rawData$yi
          ai_adj$ai_krg_adj_geo <- (rawData$ao_i - ai_adj$db_i_krg_adj_geo)*rawData$yi
          
          ai_adj$ao_i_w_adj_geo <- ai_adj$ai_w_adj_geo/rawData$yi
          ai_adj$ao_i_idw_adj_geo <- ai_adj$ai_idw_adj_geo/rawData$yi 
          ai_adj$ao_i_krg_adj_geo <- ai_adj$ai_idw_adj_geo/rawData$yi 
          
          ai_adj_geo[[k]] <- ai_adj
          names(ai_adj_geo)[k] <- c(paste("geodesic", k, "period"))
          
          # Ajuste geodésico - area ponderada
          bit_SC <- rawData[Nyear]
          for (i in 1:nrow(bit_SC)){
            for (j in 1:length(Nyear)){
              bit_SC[i,j] <- ((ai_adj$ao_i_w_adj_geo[i]) + bo_t[j])*rawData$yi[i]
            }
          }
          bit_LMf_w_adj[[k]] <- bit_SC[(idx1+1):idx2]
          names(bit_LMf_w_adj)[k] <- c(paste("bi_adj_geo", k, "period"))
          
          ## Ajuste geodésico - idw
          bit_SC <- rawData[Nyear]
          for (i in 1:nrow(bit_SC)){
            for (j in 1:length(Nyear)){
              bit_SC[i,j] <- ((ai_adj$ao_i_idw_adj_geo[i]) + bo_t[j])*rawData$yi[i]
            }
          }
          bit_LMf_idw_adj[[k]] <- bit_SC[(idx1+1):idx2]
          names(bit_LMf_idw_adj)[k] <- c(paste("bi_adj_geo", k, "period"))
          
          ## Ajuste geodésico - krg
          bit_SC <- rawData[Nyear]
          for (i in 1:nrow(bit_SC)){
            for (j in 1:length(Nyear)){
              bit_SC[i,j] <- ((ai_adj$ao_i_krg_adj_geo[i]) + bo_t[j])*rawData$yi[i]
            }
          }
          bit_LMf_krg_adj[[k]] <- bit_SC[(idx1+1):idx2]
          names(bit_LMf_krg_adj)[k] <- c(paste("bi_adj_geo", k, "period"))
          
        } else {
          bit_LM_w_adj[[k]] <- bit_LM[Nyear][(idx1+1):idx2]
          bit_LMf_w_adj[[k]] <- bit_LMf[Nyear][(idx1+1):idx2]
          names(bit_LMf_w_adj)[k] <- c(paste("bi_adj_geo", k, "period"))
          
          bit_LM_idw_adj[[k]] <- bit_LM[Nyear][(idx1+1):idx2]
          bit_LMf_idw_adj[[k]] <- bit_LMf[Nyear][(idx1+1):idx2]
          names(bit_LMf_idw_adj)[k] <- c(paste("bi_adj_geo", k, "period"))
          
          bit_LM_krg_adj[[k]] <- bit_LM[Nyear][(idx1+1):idx2]
          bit_LMf_krg_adj[[k]] <- bit_LMf[Nyear][(idx1+1):idx2]
          names(bit_LMf_krg_adj)[k] <- c(paste("bi_adj_geo", k, "period"))}
      } 
    }
    if("aux" %in% colnames(B_LMf_w)){
      # Balance acumulado de los balances anuales Bm_a calculado 
      B_LMf_w$B_LMfadj_w_tot <- cumsum(B_LMf_w$B_LMfadj_w_a) # metodo area ponderada
      B_LMf_idw$B_LMfadj_idw_tot <- cumsum(B_LMf_idw$B_LMfadj_idw_a) # metodo idw
      B_LMf_krg$B_LMfadj_krg_tot <- cumsum(B_LMf_krg$B_LMfadj_krg_a) #metodo kriging
      
      # Construccion de matrices para la interpolacion
      # Area Ponderada
      bit_LMf_w_adj <- do.call(cbind, bit_LMf_w_adj)
      colnames(bit_LMf_w_adj) <- substr(colnames(bit_LMf_w_adj), 21, nchar(colnames(bit_LMf_w_adj))) # Remover los 21 primeros caracteres de los nombres de las columnas
      anos_faltantes <- setdiff(names(bit_LMf[Nyear]), names(bit_LMf_w_adj))
      for (ano in anos_faltantes) { # Rellenar los años faltantes en B usando los valores de A
        bit_LMf_w_adj[[ano]] <- bit_LMf[[ano]]
      }
      bit_LMf_w_adj <- bit_LMf_w_adj[, sort(names(bit_LMf_w_adj))]
      bit_LMf_w_adj <- cbind(rawData[c(1:3)], bit_LMf_w_adj) 
      att <- drop_na(bit_LMf_w_adj, x)
      addWorksheet(wb, "out_bit_LMf_w_adj")
      writeData(wb, sheet = "out_bit_LMf_w_adj", att)
  
      # IDW
      bit_LMf_idw_adj <- do.call(cbind, bit_LMf_idw_adj)
      colnames(bit_LMf_idw_adj) <- substr(colnames(bit_LMf_idw_adj), 21, nchar(colnames(bit_LMf_idw_adj))) # Remover los 21 primeros caracteres de los nombres de las columnas
      anos_faltantes <- setdiff(names(bit_LMf[Nyear]), names(bit_LMf_idw_adj))
      for (ano in anos_faltantes) { # Rellenar los años faltantes en B usando los valores de A
        bit_LMf_idw_adj[[ano]] <- bit_LMf[[ano]]
      }
      bit_LMf_idw_adj <- bit_LMf_idw_adj[, sort(names(bit_LMf_idw_adj))]
      bit_LMf_idw_adj <- cbind(rawData[c(1:3)], bit_LMf_idw_adj) 
      att <- drop_na(bit_LMf_idw_adj, x)
      addWorksheet(wb, "out_bit_LMf_idw_adj")
      writeData(wb, sheet = "out_bit_LMf_idw_adj", att)
      
      # Kriging
      bit_LMf_krg_adj <- do.call(cbind, bit_LMf_krg_adj)
      colnames(bit_LMf_krg_adj) <- substr(colnames(bit_LMf_krg_adj), 21, nchar(colnames(bit_LMf_krg_adj))) # Remover los 21 primeros caracteres de los nombres de las columnas
      anos_faltantes <- setdiff(names(bit_LMf[Nyear]), names(bit_LMf_krg_adj))
      for (ano in anos_faltantes) {# Rellenar los años faltantes en B usando los valores de A
        bit_LMf_krg_adj[[ano]] <- bit_LMf[[ano]]
      }
      bit_LMf_krg_adj <- bit_LMf_krg_adj[, sort(names(bit_LMf_krg_adj))]
      bit_LMf_krg_adj <- cbind(rawData[c(1:3)], bit_LMf_krg_adj) 
      att <- drop_na(bit_LMf_krg_adj, x)
      addWorksheet(wb, "out_bit_LMf_krg_adj")
      writeData(wb, sheet = "out_bit_LMf_krg_adj", att)
      
      # Mean Adjusted Mass Balance 
      ai_df_w_adj <- cbind(XY_coords, value = apply(bit_LMf_w_adj[Nyear], 1, mean, na.rm = T)) #dataframe de valores ai
      ai_df_w_adj$value[is.nan(ai_df_w_adj$value)] <- NA # Replace NaN with NA
      
      ai_df_idw_adj <- cbind(XY_coords, value = apply(bit_LMf_idw_adj[Nyear], 1, mean, na.rm = T)) #dataframe de valores ai
      ai_df_idw_adj$value[is.nan(ai_df_idw_adj$value)] <- NA # Replace NaN with NA
      
      ai_df_krg_adj <- cbind(XY_coords, value = apply(bit_LMf_krg_adj[Nyear], 1, mean, na.rm = T)) #dataframe de valores ai
      ai_df_krg_adj$value[is.nan(ai_df_krg_adj$value)] <- NA # Replace NaN with NA 
      
      aux=1
      B_LMf_w <- subset(B_LMf_w, select = -aux)
      B_LMf_idw <- subset(B_LMf_idw, select = -aux)
      B_LMf_krg <- subset(B_LMf_krg, select = -aux)
    }
    ##################################################################################
    #=== CALCULATION OF DISTRIBUTED MASS BALANCE FOR DATA WITH GEODETIC ADJUSTMENT ===
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 1. Weighted area/altitude range method for geodetically adjusted data from the Nonlinear Lliboutry Model
    # From the output files *.csv we can compute glacier variables
    if(aux==1) {
      WAreaR = WArea(bit_LMf_w_adj, m, shapes_path, dem, crsID, shap, Nyear, aux)
      bit_LMf_w_adj_rg <- WAreaR$bit_LMf_w_adj_rg
      B_LMf_w_adj_rg <- WAreaR$B_LMf_w_adj_rg
      B_LMf_w <- cbind(B_LMf_w, B_LMf_w_adj_rg)
      # Export the DataFrame obtained through the Weighted Area Method.
      addWorksheet(wb, "out_bit_LMf_w_adj_rg_a")
      writeData(wb, sheet = "out_bit_LMf_w_adj_rg_a", bit_LMf_w_adj_rg)
      
      addWorksheet(wb, "out_B_LMf_w")
      writeData(wb, sheet = "out_B_LMf_w", B_LMf_w)
      aux=1}
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 2. Inverse Distance Weighted interpolation Method for geodetically adjusted data from the Nonlinear Lliboutry Model
    if(aux==1) {
      IdwGBR = IdwGB(bit_LMf_idw_adj, dem, crsID, brick_raw, Nyear, aux)
      bit_LMf_idw_adj_sp = IdwGBR$bit_LMf_idw_adj_sp
      # Creation and masking of RasterBrick objects according to the annual area
      att <- BrickBuilder(list(LMf_idw_adj_sp = bit_LMf_idw_adj_sp), brick_raw, Nyear)
      brick_list$LMf_idw_adj_sp <- lapply(setNames(1:length(Nyear), Nyear), function(i){
        raster::mask(att$LMf_idw_adj_sp[[i]], as(st_zm(st_as_sf(shape_list[[shap[i]]])), "Spatial"))
      }) %>% brick()
      
      # Calculation of the distributed mass balance by IDW
      B_LMf_idw_adj_sp <- data.frame(area_glc = sapply(1:length(Nyear), function(i){
        # Calculate area using sf for sf objects
        sum(st_area(shape_list[[shap[i]]]))
      }),
      LMf_idw_adj_sp = sapply(1:length(Nyear), function(i){
        cellStats(brick_list$LMf_idw_adj_sp[[i]], stat='sum', na.rm=TRUE, asSample=TRUE)
      })
      )
      
      B_LMf_idw_adj_sp$B_LMf_idw_adj_sp_a <- B_LMf_idw_adj_sp$LMf_idw_adj_sp*(Resolution^2)/B_LMf_idw_adj_sp$area_glc
      B_LMf_idw_adj_sp$B_LMf_idw_adj_sp_tot <- cumsum(B_LMf_idw_adj_sp$B_LMf_idw_adj_sp_a)
      B_LMf_idw <- cbind(B_LMf_idw, B_LMf_idw_adj_sp[,2:4])
      
      # Export the DataFrame obtained
      addWorksheet(wb, "out_bit_LMf_idw_adj_sp")
      writeData(wb, sheet = "out_bit_LMf_idw_adj_sp", bit_LMf_idw_adj_sp)
      
      addWorksheet(wb, "out_B_LMf_idw")
      writeData(wb, sheet = "out_B_LMf_idw", B_LMf_idw)
      aux=1}
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Kriging method for geodetically adjusted data from the Nonlinear Lliboutry Model
    if(aux==1) {
      KrgGBR = KrgGB(bit_LMf_krg_adj, dem, crsID, brick_raw, Nyear, aux)
      varplot1_LMfadj_sp = KrgGBR$varplot1_LMfadj_sp
      varplot2_LMfadj_sp = KrgGBR$varplot2_LMfadj_sp
      bit_LMf_krg_adj_sp = KrgGBR$bit_LMf_krg_adj_sp
      # Variogram plot of adjusted annual balance
      print(
        ggplot(varplot1_LMfadj_sp, aes(x = dist, y = gamma)) +
          geom_line(data = varplot2_LMfadj_sp, size = 0.75) +
          geom_point(colour = "red") +
          facet_wrap(.~year, scales = "free") + 
          ggtitle("SemiVariograms of Lliboutry fill - adjusted") +
          xlab("Range distance (m)") + ylab(expression("Semivariance (m w.e.a"^-1*")"))
      )
      
      # Creation and masking of RasterBrick objects according to the annual area
      att <- BrickBuilder(list(LMf_krg_adj_sp = bit_LMf_krg_adj_sp), brick_raw, Nyear)
      brick_list$LMf_krg_adj_sp <- lapply(setNames(1:length(Nyear), Nyear), function(i){
        raster::mask(att$LMf_krg_adj_sp[[i]], as(st_zm(st_as_sf(shape_list[[shap[i]]])), "Spatial"))
      }) %>% brick()
      
      # Calculation of the distributed mass balance by IDW
      B_LMf_krg_adj_sp <-data.frame(
        area_glc = sapply(1:length(Nyear), function(i){
          #raster::area(shape_list[[shap[i]]])
          # Calculate area using sf for sf objects
          sum(st_area(shape_list[[shap[i]]]))
        }),
        LMf_krg_adj_sp = sapply(1:length(Nyear), function(i){
          cellStats(brick_list$LMf_krg_adj_sp[[i]], stat='sum', na.rm=TRUE, asSample=TRUE)
        })
      )
      
      B_LMf_krg_adj_sp$B_LMf_krg_adj_sp_a <- B_LMf_krg_adj_sp$LMf_krg_adj_sp*(Resolution^2)/B_LMf_krg_adj_sp$area_glc # Balance krig ponderado a la talla del pixel
      B_LMf_krg_adj_sp$B_LMf_krg_adj_sp_tot <- cumsum(B_LMf_krg_adj_sp$B_LMf_krg_adj_sp_a)
      B_LMf_krg <- cbind(B_LMf_krg, B_LMf_krg_adj_sp[2:4])
      
      # Export the DataFrame obtained
      addWorksheet(wb, "out_bit_LMf_krg_adj_sp")
      writeData(wb, sheet = "out_bit_LMf_krg_adj_sp", bit_LMf_krg_adj_sp)
      
      addWorksheet(wb, "out_B_LMf_krg")
      writeData(wb, sheet = "out_B_LMf_krg", B_LMf_krg)
      aux=1}
  }
  ################################################################################
  #====================== SET OF GRAPHS AND RASTER OUTPUTS =======================
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Contour lines of the average balance (alpha_i)
  source("./functions/Contour_tmap.R")
  if(aux==1){
    mapa_contornos_idw = ContourPlot_tmap(
      ai_df_idw_adj,
      "idw",
      "Balance de masa promedio ajustado",
      brick_raw,
      shape_list,
      dem = dem,
      crsID = crsID
    )
    
    mapa_contornos_krg = ContourPlot_tmap(
      ai_df_krg_adj,
      "krig",
      "Balance de masa promedio ajustado",
      brick_raw,
      shape_list,
      dem = dem,
      crsID = crsID
    )
  } else {
    mapa_contornos_idw = ContourPlot_tmap(
      ai_df,
      "idw",
      "Balance de masa promedio",
      brick_raw,
      shape_list,
      dem = dem,
      crsID = crsID
    )
    mapa_contornos_krg = ContourPlot_tmap(
      ai_df,
      "krig",
      "Balance de masa promedio",
      brick_raw,
      shape_list,
      dem = dem,
      crsID = crsID
    )
  }
  ai_df_z <- ai_df
  ai_df_z$z <- terra::extract(rast(dem), ai_df_z[c("x","y")], method="bilinear")[, 2]
  ai_df_z <- ai_df_z[!is.na(ai_df_z$value), ]
  balance_altura <- plot_ly(data = ai_df_z, 
                            x = ~z, 
                            y = ~value,
                            type = 'scatter',
                            mode = 'markers',
                            marker = list(size = 6, opacity = 0.7, color = ~z, colorscale = 'Viridis'),
                            text = ~paste("Altitud:", round(z,1), "m<br>Balance:", round(value,3)),
                            hoverinfo = "text") %>%
    layout(title = "Relación Balance ∼ altitud",
           xaxis = list(title = "Altitud (m)"),
           yaxis = list(title = "Balance"),
           hovermode = "closest")
  yearly_b_df <- as.data.frame(brick_list$LMf, xy = TRUE, na.rm = TRUE)
  yearly_b_df$z <- terra::extract(rast(dem), yearly_b_df[c("x","y")], method="bilinear")[, 2]
  
  balance_altura_anual <- yearly_b_df %>%
    select(z, starts_with("X20")) %>%
    pivot_longer(cols = starts_with("X20"),
                 names_to = "fecha",
                 values_to = "balance") %>%
    mutate(año = substr(fecha, 2, 5)) %>%
    drop_na(z, balance) %>%
    group_by(año) %>%
    arrange(z) %>%
    mutate(tendencia = predict(loess(balance ~ z, span = 0.4))) %>%   
    ungroup() %>%
    
    plot_ly(x = ~z, y = ~balance, color = ~año, colors = "Paired",
            type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.7),
            text = ~paste("Año:", año,
                          "<br>Altitud:", round(z), "m",
                          "<br>Balance:", round(balance, 3)),
            hoverinfo = "text",
            name = ~año) %>%
    
    add_lines(x = ~z, y = ~tendencia,
              color = ~año,
              line = list(width = 5),
              name = ~paste(año, "tendencia"),
              showlegend = TRUE) %>%
    
    layout(
      title = "Balance anual vs Altitud (2000–2009)",
      xaxis = list(title = "Altitud (m)"),
      yaxis = list(title = "Balance (m w.e./año)"),
      hovermode = "closest",
      legend = list(title = list(text = "Año"))
    )
  source("./functions/RasterPlot_tmap.R")
  mapa_raster = RasterPlot_tmap(brick_list$LMf, FALSE, "Balance de Masa por año", brick_list$Raw, Nyear, shape_list, crsID)
  source("./functions/Graph_NonAdj_St.R")
  source("./functions/Graph_Adj_St.R")
  if(aux==0){ 
    Graph_R <- Graph_NonAdj(att, data4model, bo_t, Nyear, bit_LM, 
                            rawData, B_LMf_w, B_LMf_idw, B_LMf_krg, 
                            e_obs_mod, temp_dir, wb)
  }else{ 
    Graph_Adj(att, data4model, bo_t, Nyear, bit_LM, rawData, 
              bit_LMf_w_adj, bit_LMf_idw_adj, bit_LMf_krg_adj,
              B_LMf_w, B_LMf_idw, B_LMf_krg,
              e_obs_mod, temp_dir, wb)
  }
  return(list(
    mapa_contornos_idw = mapa_contornos_idw,
    mapa_contornos_krg = mapa_contornos_krg,
    mapa_raster = mapa_raster,
    balance_altura = balance_altura,
    balance_altura_anual = balance_altura_anual,
    geodetic_used = aux
  )
  )
  
}