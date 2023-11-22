# Author: Yick Hin Ling, Carl Wu lab, Johns Hopkins University
# Email: yhinling@gmail.com

###########################################################################################################################################################################
#Get the common name of files with filename extension
common_name <- function(x) {
  # Remove the filename extension
  x <- tools::file_path_sans_ext(x)
  
  # Function to split the file name into parts based on "_", " ", "-"
  splitName <- function(n) strsplit(n, "(?<=_| |\\-)", perl = TRUE)[[1]]
  
  # Split each filename into parts
  x_split <- lapply(x, splitName)
  
  # The length of the longest split filename
  max_len <- max(sapply(x_split, length))
  
  # Ensure all split filenames are the same length (for easy comparison)
  x_split <- lapply(x_split, function(x) {
    c(x, rep("", max_len - length(x)))
  })
  
  # Transpose the list of split filenames
  x_split_t <- do.call(rbind, x_split)
  
  # A function to check if an element is common to a column of the transposed list
  isCommon <- function(x) {
    # Remove leading and trailing white spaces and delimiters
    x <- trimws(gsub("^[_ -]+|[_ -]+$", "", x))
    # If the column contains only digits or is empty, return FALSE
    if(all(grepl("^\\d+$", x)) | all(x == "")) {
      return(FALSE)
    }
    # Check if all elements in the column are the same, ignoring leading/trailing spaces and delimiters
    all(trimws(gsub("^[_ -]+|[_ -]+$", "", x)) == trimws(gsub("^[_ -]+|[_ -]+$", "", x[1])))
  }
  
  # Find the common parts
  common_cols <- apply(x_split_t, 2, isCommon)
  
  # Concatenate the common parts
  common_name <- paste(x_split_t[1, common_cols], collapse = "")
  
  # If the first or last characters are " ", "_", or "-", remove them
  common_name <- sub("^[_ -]+", "", common_name)
  common_name <- sub("[_ -]+$", "", common_name)
  
  return(common_name)
}
###########################################################################################################################################################################
# Define the function to keep the part of the string after the last space or underscore
extract_final_substring <- function(labels) {
  sapply(labels, function(label) {
    stringr::str_extract(label, "(?<=_|\\s)[^_\\s]+$")
  })
}
###########################################################################################################################################################################
# Merge track list insdie a trackll, input is a list of list
track.merge <- function(list_of_list, id = "Trajectory"){
  lapply(list_of_list, function (x){
    names(x) <- seq_along(x)
    y <- dplyr::bind_rows(x, .id = id)
    if(!is.null(id)){
      y[[id]] <- as.numeric(y[[id]])
    }
    return(y)
  })
}
###########################################################################################################################################################################
# This function calculate the displacement and angle from the mask_centroid
# Input should contains list with variables "x", "y", "mask_centroid_x" and "mask_centroid_y"
# Angle is in -180 to 180. From x-axis, positive: anti-clockwise, negative: clockwise
centroid.disp.angle <- function(input){
  lapply(input, function(x){
    # Calculate the displacement from the centroid to the detection
    to_mask_centroid_disp <- sqrt( ( x[["x"]] - x[["mask_centroid_x"]] )^2 + 
                                     ( x[["y"]] - x[["mask_centroid_y"]] )^2 
    )
    
    # Calculate the angle from the centroid to the detection
    # https://stackoverflow.com/questions/2676719/calculating-the-angle-between-the-line-defined-by-two-points
    xy <- cbind(x[["x"]], x[["y"]])
    mask_centroid_xy <- cbind(x[["mask_centroid_x"]], x[["mask_centroid_y"]])
    delta_x <- (xy-mask_centroid_xy)[,1]
    delta_y <- (xy-mask_centroid_xy)[,2]
    to_mask_centroid_theta <- atan2(delta_y, delta_x)*180/pi #From x-axis, positive: anti-clockwise, negative: clockwise
    
    cbind(to_mask_centroid_disp, to_mask_centroid_theta)
  })
}
###########################################################################################################################################################################
# Function to calculate KDE
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
###########################################################################################################################################################################
# Calculate the distance from one point to a line
# The following function can be used to calculate the distance d of the point a from the line defined by the two points b and c
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 
###########################################################################################################################################################################
# This function calculates the angle between two vector, each vector contains c(x, y)
# Full = 0 to 360. if false, -180 to 180
# If Degree is false, return radian
# From A to B
angle2 <- function(A, B, Full = TRUE, Degree = TRUE){
  A <- unlist(A)
  B <- unlist(B)
  
  theta <- (atan2(B[2], B[1]) - atan2(A[2], A[1]))
  
  if(Full){
    if(theta < 0 & !is.na(theta)){
      theta <-  2*pi + theta
    }
  }
  
  if(Degree){
    return(theta*180/pi)
  }else{
    return(theta)
  }    
}
###########################################################################################################################################################################
# Counting the number of neighbors around each point (https://github.com/LKremer/ggpointdensity)
count_neighbors <- function(x, y, r2, xy) {
  .Call("count_neighbors_", x, y, r2, xy, "ggpointdensity")
}    #This refers to a "C" function inside ggpointdensity

# Implementation of count_neighbors in R. Not actually used, just for clarity
count_neighbors_r <- function(x, y, r2, xy) {
  yx <- 1 / xy
  sapply(1:length(x), function(i) {
    sum((yx * (x[i] - x) ^ 2) + (xy * (y[i] - y) ^ 2) < r2)
  })
}
###########################################################################################################################################################################
# cbind without equal row number
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
###########################################################################################################################################################################
# Function to draw a circle with radius = 1
circleFun <- function(center = c(0,0),diameter = 2, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
###########################################################################################################################################################################
