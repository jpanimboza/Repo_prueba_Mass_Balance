######### Distribution of Lliboutry Balance using the Inverse Weighted Distance Method ##########
#This script uses the IDW interpolation method which is a type of deterministic method 
#for multivariate interpolation with a known scattered set of points.
#for crossvalidation of IDW
#https://www.geo.fu-berlin.de/en/v/soga-r/Advances-statistics/Geostatistics/Inverse-Distance-Weighting-IDW/Model-Selection-via-Cross-Validation-for-IDW/index.html

IdwGB <- function(df, dem, crsID, brick_raw, Nyear, aux){
  #message(df)
  # Load the Digital Elevation Model (DEM)
  cropped_DEM <- raster(dem)
  # Cutting the DEM according to the study area
  cropped_DEM <- resample(cropped_DEM, brick_raw)
  cropped_DEM <- crop(cropped_DEM, brick_raw)
  # Converting the DEM into a SpatialPixelsDataFrame
  spx <- as(cropped_DEM, "SpatialPixelsDataFrame")
  
  att_pred <- df #-> Dataframe with dimensions equal to those of the df.
  idw.CV <- list()
  
  stat.idw.CV <- list()
  
  for (i in 1:length(Nyear)) {
    # Create SpatialPointsDataFrame
    # Where column 1 = x, column 2 = y, and column 3 = bm (mass balance)
    points <- df[c("x","y", Nyear[i])]
    points <- na.omit(points)
    names(points) <- c("x","y","bm")
    coordinates(points) = ~ x + y
    crs(points) <- crsID
    
    # Perform IDW interpolation
    idw.pred <- idw(formula = bm~1, locations = points, newdata = spx, idp = 2)
    neighbors <- length(points) - 1
    idw <- gstat(formula = bm~1, data = points, nmax = neighbors, set = list(idp = 2))
    idw.CV[[Nyear[i]]] <- gstat.cv(idw)
    
    
    # Calculo de estadistica para encontrar el error de ajuste idw
    stat.CV <- stat.desc(idw.CV[[Nyear[i]]])
    rmse.CV <- sqrt(sum((idw.CV[[Nyear[i]]]$residual)^2, na.rm=T)/length(idw.CV[[Nyear[i]]]$residual)) # El valor determinante sera el RMS pero se evaluaron secundariamente los otros paremetros
    Expvar <- 1-var(idw.CV[[Nyear[i]]]$residual, na.rm=T)/var(points$bm, na.rm=T) # Varianza de los datos explicada por el variograma teorico
    MPSE <- mean(subset(idw.CV[[Nyear[i]]]$residual, idw.CV[[Nyear[i]]]$residual != "NaN")^2) # Mean square prediction error MSPE, ideally small
    CORR_OP <- cor(idw.CV[[Nyear[i]]]$observed, idw.CV[[Nyear[i]]]$observed - idw.CV[[Nyear[i]]]$residual) # correlation observed and predicted, ideally 1
    CORR_OR <- cor(idw.CV[[Nyear[i]]]$observed - idw.CV[[Nyear[i]]]$residual, idw.CV[[Nyear[i]]]$residual) # correlation predicted and residual que prueba una minima diferencia entre los valores, ideally 0
    
    # Arma matriz de estadistica de ajuste kriging se guarda en comp.cv
      stat.CV <- matrix(stat.CV$residual, ncol=1)
      stat.CV <- rbind(stat.CV, rmse.CV, Expvar, MPSE, CORR_OP, CORR_OR)
      stat.CV <- stat.CV[-c(2:8)]
      #stat.CV <- stat.CV[9:19]
      stat.CV <- stat.CV[-4]
      #stat.CV <- stat.CV[-3]
      stat.CV <- stat.CV[-6]
      #stat.CV <- stat.CV[-5]
      stat.idw.CV[[Nyear[i]]] <- matrix(stat.CV, ncol=1)
      colnames(stat.idw.CV[[Nyear[i]]]) <- c("Idw.residuals")
      rownames(stat.idw.CV[[Nyear[i]]]) <- c("n_obs","Mean error ME","Mean SE", "Var","Std dev", "RMS", "GlobVar_explicada", "Mean square prediction error MSPE", "CoefCorr obs/pred", "CoefCorr obs/resid")
    
    att_pred[,i+3] <- idw.pred$var1.pred
  }
  
  # This conditional assigns variables to the created dataframes, depending on the value of the auxiliary variable.
  if(aux==0){
    bit_LMf_idw <- att_pred
    idw.CV_LMf <- idw.CV
    stat.idw.CV_LMf <- stat.idw.CV
    return(list(
      bit_LMf_idw = bit_LMf_idw,
      idw.CV_LMf = idw.CV_LMf,
      stat.idw.CV_LMf = stat.idw.CV_LMf
    ))
    } else {
      bit_LMf_idw_adj_sp <- att_pred
      return(list(bit_LMf_idw_adj_sp = bit_LMf_idw_adj_sp))
    }
  }
