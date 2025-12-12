#============== Lliboutry Balance Distribution by Kriging Interpolation ================
KrgGB <- function(df, dem, crsID, brick_raw, Nyear, aux){
  cropped_DEM <- raster(dem) #-> Load the Digital Elevation Model (DEM).
  # Cutting the DEM according to the study area.
  cropped_DEM <- resample(cropped_DEM, brick_raw)
  cropped_DEM <- crop(cropped_DEM, brick_raw)
  spx <- as(cropped_DEM, "SpatialPixelsDataFrame") #-> Converting the DEM into a SpatialPixelsDataFrame.
  
  att_pred <- df #-> Dataframe with dimensions equal to those of the df.
  varplot1 <- data.frame()
  varplot2 <- data.frame()
  krg.CV <- list()
  
  stat.krg.CV <- list() # estadisticas para estimar la incertidumbre
  
  for (i in 1:length(Nyear)) {
    # Create SpatialPointsDataFrame.
    points <- df[c("x","y", Nyear[i])]
    points <- na.omit(points)
    names(points) <- c("x","y","bm") #-> Assigns names to the columns.
    coordinates(points) = ~ x + y #-> Assign coordinates.
    crs(points) <- crsID #-> Assign Coordinate Reference System.
    
    Var <- variogram(bm ~ x+y, points, cutoff = 1500)
    emp1 <- cbind(year = matrix(Nyear[i], nrow = nrow(Var)), Var)
    varplot1 <- rbind(varplot1, emp1)
    
    fit_krig <- fit.variogram(Var,
                              vgm(c("Exp","Gau","Sph"), psill = mean(Var[1:3,"gamma"]) +
                                    mean(Var[c((nrow(Var)-4):nrow(Var)),"gamma"]),
                                  range = max(Var["dist"])/3,
                                  nugget = mean(Var[1:3,"gamma"])), fit.method = 2)
    
    var.data <- variogramLine(fit_krig, maxdist = max(Var$dist))
    emp2 <- cbind(year = matrix(Nyear[i], nrow = nrow(var.data)), var.data)
    varplot2 <- rbind(varplot2, emp2)
    
    krg.CV[[Nyear[i]]] <- krige.cv(bm ~ x+y, points, model = fit_krig)
    
    # Calculo de estadistica para encontrar el error de ajuste kriging
    stat.CV <- stat.desc(krg.CV[[Nyear[i]]])
    rmse.CV <- sqrt(mean((krg.CV[[Nyear[i]]]$var1.pred - krg.CV[[Nyear[i]]]$observed)^2, na.rm=T)) # El valor determinante sera el RMS pero se evaluaron secundariamente los otros paremetros
    Expvar <- 1-var(krg.CV[[Nyear[i]]]$residual, na.rm=T)/var(points$bm, na.rm=T) # Varianza de los datos explicada por el variograma teorico
    MPSE <- mean(subset(krg.CV[[Nyear[i]]]$residual, krg.CV[[Nyear[i]]]$residual != "NaN")^2) # Mean square prediction error MSPE, ideally small
    MSNE <- mean(subset(krg.CV[[Nyear[i]]]$zscore, krg.CV[[Nyear[i]]]$zscore != "NaN")^2) # Mean square normalized error MSNE, ideally close to 1 (No debe contener NaN)
    CORR_OP <- cor(krg.CV[[Nyear[i]]]$observed, krg.CV[[Nyear[i]]]$observed - krg.CV[[Nyear[i]]]$residual) # correlation observed and predicted, ideally 1
    CORR_OR <- cor(krg.CV[[Nyear[i]]]$observed - krg.CV[[Nyear[i]]]$residual, krg.CV[[Nyear[i]]]$residual) # correlation predicted and residual que prueba una minima diferencia entre los valores, ideally 0
    
    # Arma matriz de estadistica de ajuste kriging se guarda en comp.cv
      stat.CV <- matrix(stat.CV$residual, ncol=1)
      stat.CV <- rbind(stat.CV, rmse.CV, Expvar, MPSE, MSNE, CORR_OP, CORR_OR)
      stat.CV <- stat.CV[-c(2:8)]
      stat.CV <- stat.CV[-4]
      stat.CV <- stat.CV[-6]
      stat.krg.CV[[Nyear[i]]] <- matrix(stat.CV, ncol=1)
      colnames(stat.krg.CV[[Nyear[i]]]) <- c("Krig.residuals")
      rownames(stat.krg.CV[[Nyear[i]]]) <- c("n_obs","Mean error ME","Mean SE", "Var","Std dev", "RMS", "GlobVar_explicada", "Mean square prediction error MSPE", "Mean square normalized error MSNE", "CoefCorr obs/pred", "CoefCorr obs/resid")
    
    krg.pred <- krige(bm ~ x+y, points, newdata = spx, model = fit_krig)
    att_pred[,i+3] <- krg.pred$var1.pred
  }
  
  # This conditional assigns names to the resulting dataframe based on the given variable 'aux',
  # keeps them in the stored objects history (Environment) for use outside of this function.
  if(aux==0) {
    return(list(
      varplot1_LMf    = varplot1,
      varplot2_LMf    = varplot2,
      bit_LMf_krg     = att_pred,
      krg.CV_LMf      = krg.CV,
      stat.krg.CV_LMf = stat.krg.CV
    ))

    } else {
      return(list(
        varplot1_LMfadj_sp = varplot1,
        varplot2_LMfadj_sp = varplot2,
        bit_LMf_krg_adj_sp = att_pred
      ))
    }
  }
