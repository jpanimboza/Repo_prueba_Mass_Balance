Reanal_glac_geo <- function(B_LMf_w, idx1, idx2, e_mesMB, e_dc,geo_df, k, e_repres.tot, std.w, n_obs, shape_list){
  B_glc_tot <- sum(B_LMf_w$B_LMf_w_a[(idx1+1):idx2])
  B_glc_a <- mean(B_LMf_w$B_LMf_w_a[(idx1+1):idx2], na.rm=T)*1000 # AntisanaOCC
  
  # Uncertainty of glaciological balance at m.w.e a-1 (Gerbaux et al., 2005; Zemp et al., 2013; Basantes et al., 2016)
  e_point <- sqrt(e_mesMB^2+(sum(B_LMf_w$B_LMf_w_a)^2)*(e_dc^2))*sqrt(geo_df$diff_yrs[k])
  e_repres <- e_repres.tot
  message(e_repres)
  e_w_sd <- sqrt(sum(std.w)/((n_obs[1]*geo_df$diff_yrs[k])-geo_df$diff_yrs[k]))
  e_spatial <- sqrt(e_repres^2+e_w_sd^2)
  e_area <-  e_point * as.numeric((0.5*st_length(st_boundary(shape_list[[1]]))*3)/st_area(shape_list[[1]]))
  message(e_area)
  e_Bgl_a <- sqrt(e_point^2+e_spatial^2+e_area^2)*1000 / sqrt(geo_df$diff_yrs[k]) # 660 mm w.e para Antisana (basantes et al., 2022)
                 # 340 mm w.e (mean value Zemp et al 2013)
  # Uncertainty of Geodetic Balance in m.w.e a-1 from geo_df table
  
  Delta_PoR_a <- B_glc_a-(geo_df$B_geo_a[k]*1000) # (mm w.e y-1) B_geo_a
  Delta_PoR <- Delta_PoR_a*geo_df$diff_yrs[k] # (mm w.e y-1)
  
  Porc_s <- Delta_PoR*100/(geo_df$B_geo_a[k]*1000) # %
  e_common.PoR <- sqrt((e_Bgl_a*sqrt(geo_df$diff_yrs[k]))^2+(geo_df$e_geo.yr[k]*1000*geo_df$diff_yrs[k])^2) # (m w.e)
  S_disc_reduc <- Delta_PoR/e_common.PoR
  
  Ho_0.05 <- ifelse(S_disc_reduc > -1.96, ifelse(S_disc_reduc < 1.96, "yes", "no"), "no")
  Beta_0.05 <- 100 * (pnorm(1.96 - S_disc_reduc, mean = 0, sd = 1, lower.tail = TRUE) +
                        - pnorm(-1.96 - S_disc_reduc, mean = 0, sd = 1, lower.tail = TRUE))
  epsilon_limit_0.05 <- 3.606*e_common.PoR/geo_df$diff_yrs[k]  # (m w.e)
  
  Ho_0.10 <- ifelse(S_disc_reduc > -1.64, ifelse(S_disc_reduc < 1.64, "yes", "no"), "no")
  Beta_0.10 <- 100 * (pnorm(1.64 - S_disc_reduc, mean = 0, sd = 1, lower.tail = TRUE) +
                        - pnorm(-1.64 - S_disc_reduc, mean = 0, sd = 1, lower.tail = TRUE))
  epsilon_limit_0.10 <- 2.928*e_common.PoR/geo_df$diff_yrs[k] # (m w.e)
  
  calib <- if(Ho_0.05 == "no" & Ho_0.10 == "no"){
    message('REQUIRES GEODETIC CALIBRATION')
    aux <- 1
    }else{
      if(Ho_0.05 == "yes" & Ho_0.10 == "no"){
        message('REQUIRES GEODETIC CALIBRATION')
        aux <- 1
        }else{
          if(Ho_0.05 == "no" & Ho_0.10 == "yes"){
            message('NO GEODETIC CALIBRATION REQUIRED')
            aux <- 0
            }else{
              if(Ho_0.05 == "yes" & Ho_0.10 == "yes"){
                message('NO GEODETIC CALIBRATION REQUIRED')
                aux <- 0
              }
            }
        }
      }
  aux <- aux
  message(paste0("Error glac anual para periodo de reanalisis: ", e_Bgl_a))
  message(aux)
  return(list(e_Bgl_a = e_Bgl_a,aux=aux))
  #B_glc_a <<- B_glc_a
  #B_glc_tot <<- B_glc_tot
  }
