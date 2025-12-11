#======================= Function to create bricks =======================
BrickBuilder <- function(list,brick_raw,Nyear){
  # Data array
  arraybrick <- function(x,y){
    temp <- data.frame(cbind(coordinates(brick_raw), y[,x]))
    names(temp) <- c("x","y","value")
    temp <- na.omit(temp)
  }
  
  lab <- names(list)
  
  # Creation of bricks and storage in a list
  brick_list <- lapply(setNames(lab, lab), function(z){
    lapply(setNames(Nyear, Nyear), arraybrick, y = list[[z]]) %>% 
      lapply(RasterBuilder) %>% brick
  })
}
