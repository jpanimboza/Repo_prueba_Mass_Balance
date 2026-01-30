#================== Function of the Nonlinear Lliboutry Model ==================
NLM <- function(rawData, shapes_path, dem, Resolution, crsID, N.Obs,summit,wb){
  # Load observed balance data
  rawData <- data.frame(read.delim(rawData, header=F, sep="\t", dec="."))
  names(rawData)<-c("date","x","y","z","value") # , "survey_start"
  rawData$date <- lubridate::ymd(rawData$date, truncated = 2) #-> Convert dates from character format to date objects and expect it to be in 'year-month-day' format.
  
  coordinates(rawData) = ~x + y #-> Assign coordinates.
  proj4string(rawData) <- CRS(crsID) #-> Assign Coordinate Reference System.
  rawData$z <- terra::extract(raster(dem), rawData, method = 'bilinear')
  rawData <-raster::as.data.frame(rawData)[,-7]
  
  addWorksheet(wb, "out_bit_observ")
  writeData(wb, sheet = "out_bit_observ", rawData)

  # Number of years of study
  Nyear <- unique(rawData$date)
  
  # Loading of the annual surfaces (shapefile)
  shape_list <- list.files(shapes_path, pattern = "\\.shp$")
  shape_list <- sapply(shape_list, function(x){
    paste0(shapes_path,x)
  })
  shape_list <- lapply(shape_list, st_read)
  
  # Assigning names to the coverages
  shap <- sapply(1:length(Nyear), function(x){
    shap <- grep(as.character(year(Nyear)[x]), names(shape_list), value = T)
    ifelse(length(shap) == 0, "NA", shap)
  })
  shap[1] <- names(shape_list)[1]
  
  for(i in 2:length(shap)){
    ifelse(shap[i] =="NA", shap[i] <- shap[i-1],  shap[i] <- shap[i])
  }
  
  # Spatial boundaries of the calculated rasters, assigned based on the glacier surface
  xylim <- as.vector(sf::st_bbox(shape_list[[1]])) #ahora
  xylim <- matrix(xylim, nrow=2, ncol=2)
  
  xylim[,1] <- Resolution*floor(xylim[,1]/Resolution)
  xylim[,2] <- Resolution*ceiling(xylim[,2]/Resolution)
  
  # Creation of empty rasters with assigned extension and resolution
  Raster <- raster(xmn=xylim[1,1], xmx=xylim[1,2], 
              ymn=xylim[2,1], ymx=xylim[2,2], res = Resolution)
  crs(Raster) <- crsID
  
  # Creation of a set of rasters (brick) for the balance observations
  Nyear <- as.character(Nyear)
  brick_raw <- list()
  for(i in 1:length(Nyear)){
    brick_raw[[Nyear[i]]] <- dplyr::filter(rawData, date == Nyear[i])
  }
  
  # Function to create rasters by averaging observations within each pixel (resolution)
  RasterBuilder <<- function(data){
    coordinates(data) = ~x + y
    data <- rasterize(data, Raster, "value", fun = mean)
  }
  
  # Brick creation
  brick_raw <- lapply(brick_raw, RasterBuilder) %>% brick
  coordinates(rawData) = ~x + y
  
  # Extraction of x, y, z coordinates for creating a set of rasters (brick) for balance calculations
  brick_coords <- list()
  brick_coords[["x"]] <- rasterize(rawData, Raster, rawData@coords[,"x"], 
                                   fun = mean)
  brick_coords[["y"]] <- rasterize(rawData, Raster, rawData@coords[,"y"], 
                                   fun = mean)
  brick_coords[["z"]] <- rasterize(rawData, Raster, "z", fun = mean)
  brick_coords <- brick(brick_coords)
  
  # Assigning pixel values to the observations table
  rawData <- data.frame(cbind(values(brick_coords),
                              values(brick_raw)))
  
  XY_coords <- cbind.data.frame(coordinates(brick_raw))
  
  rawData$x[!is.na(rawData$x)] <- XY_coords$x[!is.na(rawData$x)] #antes
  rawData$y[!is.na(rawData$y)] <- XY_coords$y[!is.na(rawData$y)] #antes
  colnames(rawData)[-c(1:3)] <- Nyear
  
  # Data counting
  rawData$n <- apply(!is.na(rawData[Nyear]),1,sum)
  rawData$n <- ifelse(rawData$n == 0, NA, rawData$n)
  
  # Calculation of standard deviation for polynomial fitting
  for(k in 1:2){
    # Final data count
    rawData$nf <- apply(!is.na(rawData[Nyear]),1,sum)
    rawData$nf <- ifelse(rawData$nf == 0, NA, rawData$nf)
    
    rawData$sd4model <- apply(rawData[Nyear],1, sd, na.rm = T)
    
    # Outliers of up to 2 standard deviations that do not fit the model are eliminated
    att <- which(rawData$sd4model>2*mean(rawData$sd4model, na.rm=T), arr.ind = T)
    ifelse(length(att)==0, NA, rawData$sd4model[att] <- NA)
    
    # Data: cells with 10 or more years of observation.
    # Calculation of average bij (raw) and standard deviation (sd) for model
    rawData$meanbij <- apply(rawData[Nyear],1, mean, na.rm = T)
    data4model <- dplyr::select(rawData, z, nf, sd4model) %>%
      dplyr::filter(nf > N.Obs)
    ifelse(summit==0, NA, data4model <- rbind(c(summit, 10, 0), data4model))

    #  Second-degree polynomial regression
    fit = lm(sd4model ~ poly(z, 2), data = data4model)
    
    # Prediction of standard deviation (sd) values as a function of the polynomial
    rawData$sdi <- predict(fit, data.frame(z = rawData$z))
    
    # Calculation of the yi factor
    rawData$yi <- rawData$sdi/max(rawData$sdi, na.rm = T)
    
    # Calculation of the initial parameter ai from bij (observations) 
    df <- drop_na(rawData,x)
    df <- df[Nyear]/df$yi # Mathematical artifice for standardizing the method
    O <- c(as.matrix(df))
    O <- matrix(O[!is.na(O)])
    # Construction of matrices for parameter estimation using the least squares method
    O <- rbind(O, 0) # Observation matrix
    
    a <- ifelse(is.na(df),0,1) # Construction of matrix A
    for(i in 1:ncol(df)){
      att <- matrix(0, nrow = nrow(a), ncol = ncol(a))
      att[,i] <- a[,i]
      b <- diag(1, nrow = nrow(a), ncol = nrow(a))
      att <- cbind(b, att)
      ifelse(i==1, c <- att, c <- rbind(c,att))
    }
    A <- c[-which(rowSums(c)==1, arr.ind = TRUE),]
    att <- cbind(matrix(0, nrow=1, ncol=nrow(df)), matrix(1, nrow=1, ncol=ncol(df)))
    A <- rbind(A, att) 
    # Solving a system of equations by the least squares method
    # Using solve() function is possible to get an error because of the absence of an inverse matrix that is the system 
    # is exactly singular i.e., determinant equals to 0. To solve this issue you have construct q non-0 matrix.
    # To solve practically check if the (e_obs_mod_sd) is within the acceptable value <1m.If that is the case just set k=1
    # you can check by computing its determinant using the det() function:
    # det(t(A) %*% A)
    x <- solve(t(A) %*% A) %*% t(A) %*% O # , tol = 1e-19
    ao_i <- x[1:nrow(df),] # Initial estimate of ai
    bo_t <- x[(nrow(df)+1):(dim(x)[1]),] # Initial estimate of bt Beta_t
    
    rawData$ao_i[!is.na(rawData$n)] <- ao_i
    
    ###### Parameter ai and bt
    # Data retrieval by applying the inverse method of mathematical artifice
    rawData$ai <- rawData$ao_i*rawData$yi
    ai_df <- cbind(XY_coords, value = rawData$ai) #dataframe de valores ai
    #borrar ai_raster <- RasterBuilder(ai_df)
    att <- drop_na(rawData,x)
    
    if(k==2) {
      addWorksheet(wb, "out_bit_observ_prom")
      writeData(wb, sheet = "out_bit_observ_prom", att)
      addWorksheet(wb, "out_meanBeta_t_llib")
      writeData(wb, sheet = "out_meanBeta_t_llib", bo_t)}
    
    # Model bit matrix(scaled prediction)
    # Data recovery by applying the inverse mathematical artifice
    bit_SC <- rawData[Nyear]
    for (i in 1:nrow(bit_SC)){
      for (j in 1:length(Nyear)){
        bit_SC[i,j] <- ((rawData$ao_i[i]) + bo_t[j])*rawData$yi[i]
      }
    }
    # Model bit matrix (prediction)
    bit_LM <- cbind(rawData[c(1:3)], bit_SC) 
    att <- drop_na(bit_LM,x)
    
    if(k==2) {
      addWorksheet(wb, "out_predic_llib")
      writeData(wb, sheet = "out_predic_llib", att)}
    
    # Residuals
    e_obs_mod <- rawData[Nyear] - bit_LM[Nyear]
    ifelse(k==1, e_obs_mod_sd <- sd(unlist(e_obs_mod), na.rm = T), NA) # Standard deviation of residuals
    e_obs_mod <- cbind(rawData[c(1:3)], e_obs_mod)
    att <- drop_na(e_obs_mod,x)
    
    if(k==2) {
      addWorksheet(wb, "out_e_obs_mod_llib")
      writeData(wb, sheet = "out_e_obs_mod_llib", att)}
    # Cleaning of observations that can be identified as measurement errors.
    # Lliboutry allows this identification, the criterion taken is to eliminate observations
    # outside 2 standard deviations. To know how many observations were eliminated, the difference
    # between the entered observations and the printed observations (n-nf) must be taken (out_observations.csv).
    att <- rbind(which(e_obs_mod[Nyear]<(-2*e_obs_mod_sd), arr.ind = T), which(e_obs_mod[Nyear]>(2*e_obs_mod_sd), arr.ind = T))
    ifelse(k==1, ifelse(length(att)!=0,rawData[Nyear][att]<-NA, NA), NA)
  }
  
  # Matrix bit fill of predicted in observed
  bit_LMf <- cbind(rawData[c(1:3)], rawData[Nyear]) # make a new copy of your target matrix rawData[-c(29:35)]
  NAs <- which(is.na(bit_LMf[Nyear]), arr.ind = T)
  bit_LMf[Nyear][NAs] <- bit_LM[Nyear][NAs]
  att <- drop_na(bit_LMf,x)
  
  addWorksheet(wb, "out_bit_LMf")
  writeData(wb, sheet = "out_bit_LMf", att)
  
  # Beta_it matrix, annual anomaly or deviation with respect to the mean annual balance (alfa_i) for each site
  # this comes from the general equation
  beta_it <- (bit_LMf[Nyear]-rawData$ai) 
  beta_it <- cbind(rawData[c(1:3)], beta_it)
  att <- drop_na(beta_it,x)
  
  addWorksheet(wb, "out_Beta_t_llib")
  writeData(wb, sheet = "out_Beta_t_llib", att)
  
  # Storage of data for later use in raster spatialization
  return(list(
    shape_list   = shape_list,
    shap         = shap,
    Nyear        = Nyear,
    data4model   = data4model,
    rawData      = rawData,
    brick_raw    = brick_raw,      
    bit_LM       = bit_LM,
    bit_LMf      = bit_LMf,
    ai_df        = ai_df,
    ao_i         = ao_i,
    bo_t         = bo_t,
    beta_it      = beta_it,
    e_obs_mod    = e_obs_mod,
    brick_coords = brick_coords,
    XY_coords    = XY_coords
  ))
}
