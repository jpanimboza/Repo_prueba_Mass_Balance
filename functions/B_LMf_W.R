# 1. en Antisana 12 con (m<200), Arteson con (m<200) se generan rangos sin datos hacia la cumbre
# debe corregirse y copiarse linealmente el valor del ultimo rango en el que existen datos
# ver bit_LMf_r. de lo contrario da una solucion singular
#================= Function applying the Weighted Area/Altitude Range Method ==================
WArea <- function(df, m, shapes_path, dem, crsID, shap, Nyear, aux){
  #aux=0
  #df=bit_LMf
  hypso=NULL
  for(i in 1:length(shap)){
    x <- st_read(dsn = shapes_path, layer = gsub('.{4}$', '', shap[i])) #ahora
    glac.dem <- raster::mask(raster(dem), st_zm(x))
    
    # Ranges over which the mass balance will be calculated using the weighting method.
    zl <- c(m*floor(minValue(glac.dem)/m), m*ceiling(maxValue(glac.dem)/m))
    ranges <- seq(zl[1], zl[2], m)
    rang <- data.frame(ranges = rep(NA, (length(ranges)-1)))
    for (j in 1:(length(ranges)-1)) {
      rang[j,] <- paste0(paste0(ranges[j],"-"),ranges[j+1])
      }
    
    # Arrangement of data by ranges of the modeled bit
    glac.dem.p <- data.frame(rasterToPoints(glac.dem))
    colnames(glac.dem.p) <- c("x","y","elev")
    coordinates(glac.dem.p) <- ~ x+y # Importing the data into the SpatialpointDataFrame
    proj4string(glac.dem.p) <- CRS(crsID)
    
    glac.dem.pmean <- NULL
    for (k in 1:(length(ranges)[1]-1)){
      glac.dem.psset <- subset(glac.dem.p, glac.dem.p$elev > ranges[k] & glac.dem.p$elev <= ranges[k+1]) # raster dem
      s_rel <- data.frame(length(glac.dem.psset)*res(glac.dem)[1]^2)
      glac.dem.pmean <- bind_rows(glac.dem.pmean, s_rel)
      }
    colnames(glac.dem.pmean) <- Nyear[i]
    glac.dem.pmean <- cbind(rang, glac.dem.pmean)
    ifelse(i==1, hypso <- bind_cols(hypso, glac.dem.pmean), hypso <- left_join(hypso, glac.dem.pmean, by = "ranges"))
  }

  # Linear hypothesis for area changes
  
  
  # Data arrangement by ranges
  # To assemble the tables the first year of available data is taken
  x <- st_read(dsn = shapes_path, layer = gsub('.{4}$', '', shap[1])) 
  glac.dem <- raster::mask(raster(dem), st_zm(x))
  zl <- c(m*floor(minValue(glac.dem)/m), m*ceiling(maxValue(glac.dem)/m)) # Range over which the Mass Balance is to be calculated using the weighting method
  
  ranges <- seq(zl[1], zl[2], m)
  rang <- data.frame(ranges = rep(NA, (length(ranges)-1)))
  for (j in 1:(length(ranges)-1)) {
    rang[j,] <- paste0(paste0(ranges[j],"-"),ranges[j+1])
    }
  
  df_r <- df %>%
    mutate(ranges = cut(z, breaks = ranges, labels = rang$ranges))
  df_r <- df_r %>% na.omit() %>% group_by(ranges) %>%
    dplyr::summarise_all("mean")
  df_r$ranges <- as.character(df_r$ranges)
  df_r <- left_join(rang, df_r, by = "ranges")
  
  # Interpolation.
  # Identifies positions to fill in ranges where no observations are available.
  pos <- data.frame(pos1 = 1:nrow(df_r))
  pos$pos2 <- pos$pos1
  for (i in order(1:nrow(pos),  decreasing = T)) {
    pos[i,2] <- ifelse(is.na(df_r[i,2]) == F, pos[i,2], pos[i+1,2])
    }
  
  pos <- pos %>% fill(names(.),.direction = "down")
  pos$pos1 <- pos$pos1-1
  
  # This loop calculates the averages of the balance in the ranges in which there are no measurements.
  # If it is at the extremes then it repeats the previous balance; if it is in an intermediate range then it averages it among those that do have data.
  for (i in 1:nrow(df_r)) {
    for (j in 2:4) {
      df_r[i,j] <- ifelse(is.na(df_r[i,j]) == F, df_r[i,j],
                          mean(c(df_r[pos$pos1[i],j],
                                 df_r[pos$pos2[i],j])))
      }
    for (k in Nyear) {
      df_r[i,k] <- ifelse(is.na(df_r[i,k]) == F, df_r[i,k],
                          predict.lm(lm(data = data.frame(
                            x = c(df_r[pos$pos1[i],4],
                                  df_r[pos$pos2[i],4]),
                            y = c(df_r[pos$pos1[i],k],
                                  df_r[pos$pos2[i],k])), y ~ x),
                            newdata = data.frame(x = df_r[i,4])))
    }
    }
  
  ## Hypsometry
  # It is outside the conditional since it is the same hypsometry for both cases, both for the annual adjustment and the pixel adjustment.
  # Conditional to fill in missing years 
  hpyear <- sapply(1:length(Nyear), function(x){
    hpyear <- grep(as.character(year(Nyear)[x]), names(hypso), value = T)
    ifelse(length(hpyear) == 0, "NA", hpyear)
    })
  hpyear[1] <- names(hypso)[-1][1]
  for(i in 2:length(hpyear)){
    ifelse(hpyear[i] =="NA", hpyear[i] <- hpyear[i-1],  hpyear[i] <- hpyear[i])
    }
  
  hpnew <- as.data.frame(sapply(setNames(1:length(Nyear), Nyear), function(i){
    hypso[[hpyear[i]]]
    }))
  hpnew <- cbind(hypso[1], hpnew)
  
  hpsum <- colSums(hpnew[-1], na.rm = T) # Total sum of the area.
  
  df_r <- left_join(hpnew[1], df_r)
  
  ## Weighted balance calculation
  mult <- df_r[Nyear]*hpnew[Nyear]
  df_w_a <-  sweep(mult,2,hpsum,`/`)
  df_w_a <- cbind.data.frame(df_r[c(2:4)], hpnew[1], df_w_a)
  B_df_w <- data.frame(B_df_w_a = colSums(df_w_a[Nyear], na.rm = T))
  B_df_w$B_df_w_tot <- cumsum(B_df_w$B_df_w_a)
  
  # This conditional assigns variables to the created dataframes, depending on the value of the auxiliary variable.
  if(aux==0) {
    colnames(B_df_w) <- c("B_LMf_w_a","B_LMf_w_tot")
    B_df_w<-tibble::rownames_to_column(B_df_w, var = "Nyear")
    bit_LMf_w_a <- df_w_a
    B_LMf_w <- B_df_w
    stat.w.sd <- apply(df_r[Nyear], 2, sd, na.rm = T)
    return(list(
      bit_LMf_w_a = bit_LMf_w_a,
      B_LMf_w = B_LMf_w,
      stat.w.sd = stat.w.sd
    ))
  } else {
    colnames(B_df_w) <- c("B_LMf_w_adj_rg_a","B_LMf_w_adj_rg_tot")
    bit_LMf_w_adj_rg <- df_w_a
    B_LMf_w_adj_rg <- B_df_w
    return(list(
      bit_LMf_w_adj_rg = bit_LMf_w_adj_rg,
      B_LMf_w_adj_rg = B_LMf_w_adj_rg
    ))
  }
  
}
