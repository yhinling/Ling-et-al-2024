# Author: Yick Hin Ling, Carl Wu lab, Johns Hopkins University
# Email: yhinling@gmail.com

##################################  Step 1: Read DiaTrack .mat files ##################################################################################################################################################################
readDiaTrack <- function(folder){
    
    mat_files <- list.files(file.path(folder), pattern = ".mat$")
    trackll <- list()
    
    #Foreach
    cl <- parallel::makeCluster(detectCores()-1)
    doSNOW::registerDoSNOW(cl)
    iterations <- length(mat_files)
    
    cat(paste0("\nCommon Name:  ", common_name(mat_files), "\n"))
    cat(paste0("Reading number of Diatrack files: ", iterations, "\n"))
    
    pb <- pbmcapply::progressBar(max = iterations, style = "ETA")
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    
    trackll <- foreach(DiaTrack_file = 1:iterations,
                       .packages=c("R.matlab", "tidyverse", "pbmcapply", "doSNOW"), 
                       .options.snow = opts) %dopar% {
                           
                           data <- R.matlab::readMat(mat_files[DiaTrack_file])
                           track <- c()
                           
                           for (i in 1:length(data$tracks)){
                               
                               # Skip the loop if there is no detection in that particular frame
                               if (length(data$tracks[[i]][[1]][[1]]) == 0){
                                   next
                               }
                               
                               # Get the index in .mat file
                               x_index <- which(rownames(data$tracks[[i]][[1]]) == "RefinedCooX")
                               y_index <- which(rownames(data$tracks[[i]][[1]]) == "RefinedCooY")
                               z_index <- which(rownames(data$tracks[[i]][[1]]) == "RefinedCooZ")
                               Intensity_index <- which(rownames(data$tracks[[i]][[1]]) == "Intensity")
                               Successor_index <- which(rownames(data$tracks[[i]][[1]]) == "Successor")
                               Predecessor_index <- which(rownames(data$tracks[[i]][[1]]) == "Predecessor")
                               
                               x <- c(unlist(data$tracks[[i]][[1]][[x_index]]))
                               y <- c(unlist(data$tracks[[i]][[1]][[y_index]]))
                               z <- c(unlist(data$tracks[[i]][[1]][[z_index]]))
                               Intensity <- c(unlist(data$tracks[[i]][[1]][[Intensity_index]]))
                               Successor <- c(unlist(data$tracks[[i]][[1]][[Successor_index]]))
                               Predecessor <- c(unlist(data$tracks[[i]][[1]][[Predecessor_index]]))
                               Frame <- c(rep(i, length(y)))
                               
                               track <- rbind(track, cbind(Frame, x, y, z, Intensity, Successor, Predecessor))
                           }
                           
                           track <- as_tibble(track)
                           track <- cbind(Trajectory = integer(nrow(track)), track)
                           track$Trajectory[1] <- 1
                           
                           for (i in 2:nrow(track)){     
                               # If Predecessor = 0 >> New Track
                               if (track$Predecessor[i] == 0){
                                   track$Trajectory[i] <- max(track$Trajectory) + 1
                               }else{
                                   # Label track based on Predecessor number
                                   track$Trajectory[i] <- track$Trajectory[track$Frame == track$Frame[i]-1][track$Predecessor[i]]
                               }
                           }
                           
                           track <- track %>% 
                               arrange(Trajectory, Frame)
                           
                           # Remove Trajectory, Successor and Predecessor columns
                           track_list <- split(track[, -c(1, 7, 8)], track$Trajectory)
                           
                           # Concatenate the trackll list
                           trackll[[DiaTrack_file]] <- track_list
                           
                       }
    close(pb)
    parallel::stopCluster(cl) 
    # End of Foreach
    
    # Name the track list with filename
    names(trackll) <- tools::file_path_sans_ext(mat_files)
    # Reorder the names
    trackll <- trackll[gtools::mixedsort(names(trackll))]
    
    return(trackll)
}
trackll_raw <- readDiaTrack(folder)

#################################### Step 2: Read mask and landmark .tif files ########################################################################################################################################################
# MinPts: DBSCAN parameters
# eps: DBSCAN parameters

readMask <- function(folder, resolution = c(128, 128), pixel_size_um = 16/150, MinPts = 10, eps = 3, mask.plot = TRUE, mask_outline.plot = TRUE, maskName = "MASK", landmark = FALSE){
  
  # Read _MASK.tif files and save the matrix in mask_list
  mask_files <- list.files(file.path(folder), pattern = paste0(maskName, ".tif$"))
  mask_list <- list()
  for (i in seq_along(mask_files)){
    img <- EBImage::readImage(file.path(folder, mask_files[i]))
    mat <- which(img != 0, arr.ind = TRUE, useNames = FALSE)  # Find nucleus with pixel = 1
    mask_list[mask_files[i]] <- list(mat)
  }
  
  # DBSCAN to find mask
  mask <- lapply(mask_list, function(x){
    model <- fpc::dbscan(x, MinPts = MinPts, eps = eps)
    classification <- data.frame(x = x[,1], y = x[,2], mask = as.numeric(model$cluster))
    classification <- classification %>% 
      dplyr::group_by(mask) %>% 
      dplyr::mutate(mask_centroid_x = mean(x), mask_centroid_y = mean(y))
    classification <- classification[order(classification$mask),]
    return(classification)
  })
  
  # Reorder the name of the mask
  names(mask) <- gsub(paste0("_", maskName, "$"), "", tools::file_path_sans_ext(names(mask)))
  mask <- mask[gtools::mixedsort(names(mask))]
  
  cat(paste0("\nCommon Name:  ", common_name(names(mask)), "\n"))
  cat(paste0("Number of ", maskName," files: ", length(mask), "\n"))
  
  # Plot mask
  if(mask.plot){
    # Merge mask data
    mask_merged_list <- bind_rows(mask, .id = "mask_files")
    mask_merged_list$mask_files <- factor(mask_merged_list$mask_files, levels = gtools::mixedsort(names(mask)))
    
    # Covert pixel to um
    mask_merged_list$x <- mask_merged_list$x * pixel_size_um
    mask_merged_list$y <- mask_merged_list$y * pixel_size_um
    mask_merged_list$mask_centroid_x <- mask_merged_list$mask_centroid_x * pixel_size_um
    mask_merged_list$mask_centroid_y <- mask_merged_list$mask_centroid_y * pixel_size_um
    
    p <- ggplot2::ggplot(mask_merged_list, aes(x,y, color = factor(mask))) +
      geom_point() +
      geom_point(aes(mask_centroid_x, mask_centroid_y), size = 2, color = "black") +
      xlim(0,resolution[1] * pixel_size_um) + ylim(0,resolution[2] * pixel_size_um) +
      xlab("x (\u03BCm)") + ylab("y (\u03BCm)") +
      guides(color = guide_legend(title = maskName)) +
      coord_fixed(ratio = 1) +
      facet_wrap(~ mask_files)
    
    cat("Plotting...\n")
    
    print(p)
  }
  
  # Get the outline of the mask
  mask_outline <- lapply(mask, function(x) split(x, x[["mask"]]))
  
  mask_outline <- lapply(mask_outline, function(x){
    lapply(x, function(k){
      y_outline_min <- min(k[["y"]])
      y_outline_max <- max(k[["y"]])
      outline_min <- c()
      outline_max <- c()
      
      # Find min and max of x_outline for every value of y
      for(i in y_outline_min:y_outline_max){
        x_outline_min <- min(k[["x"]][ k[["y"]] == i ])
        x_outline_max <- max(k[["x"]][ k[["y"]] == i ])
        outline_min <- c(outline_min, x_outline_min, i) # i refer to y_outline
        outline_max <- c(outline_max, x_outline_max, i) # i refer to y_outline
      }
      
      outline_min <- matrix(outline_min, ncol = 2, byrow = T)
      outline_max <- matrix(outline_max, ncol = 2, byrow = T)
      
      mask_outline <- rbind(outline_max, outline_min[nrow(outline_min):1, ], outline_max[1,])
      colnames(mask_outline) <- c("x", "y")  
      as.data.frame(mask_outline)
    })
    
  })
  
  mask_outline <- lapply(mask_outline, function(x) dplyr::bind_rows(x, .id = "mask"))
  mask_outline <- lapply(mask_outline, function(x) {
    x$mask <- as.numeric(x$mask)
    return(x)
  })            
  
  # Add mask_centroid from mask data to mask_outline data, and calculate nuclear mask radius
  mask_merged_list <- dplyr::bind_rows(mask, .id = "mask_files")
  mask_outline_merged_list <- dplyr::bind_rows(mask_outline, .id = "mask_files")
  mask_outline_merged_list <- mask_outline_merged_list %>%
    dplyr::inner_join(mask_merged_list, by = c("mask_files", "mask", "x", "y")) %>%
    dplyr::group_by(mask_files, mask) %>% 
    dplyr::mutate(mask_radius_um = mean(sqrt((x - mask_centroid_x)^2 + (y - mask_centroid_y)^2)) * pixel_size_um)
  
  mask_outline <- split(mask_outline_merged_list, mask_outline_merged_list$mask_files)
  
  # Add mask_radius_um from mask_outline to mask (There maybe a cleaner way to do it...)
  mask_merged_list <- mask_merged_list %>% dplyr::inner_join(mask_outline_merged_list[, c("mask_files", "mask","mask_radius_um")], by = c("mask_files", "mask"))
  mask_merged_list <- mask_merged_list[!duplicated(mask_merged_list), ]  #inner_join() creats a lot of duplicated rows
  mask <- split(mask_merged_list, mask_merged_list$mask_files)
  mask <- lapply(mask, function(x) x[, -1])
  
  
  if(mask_outline.plot){
    # Covert pixel to um
    mask_outline_merged_list$x <- mask_outline_merged_list$x * pixel_size_um
    mask_outline_merged_list$y <- mask_outline_merged_list$y * pixel_size_um
    mask_outline_merged_list$mask_centroid_x <- mask_outline_merged_list$mask_centroid_x * pixel_size_um
    mask_outline_merged_list$mask_centroid_y <- mask_outline_merged_list$mask_centroid_y * pixel_size_um
    
    p <- ggplot2::ggplot(mask_outline_merged_list, aes(x,y, color = factor(mask))) +
      geom_path() + # Unlike geom_point(), geom_path() connects consecutive data
      geom_point(aes(mask_centroid_x, mask_centroid_y), size = 2, color = "black") +
      xlim(0,resolution[1] * pixel_size_um) + ylim(0,resolution[2] * pixel_size_um) +
      xlab("x (\u03BCm)") + ylab("y (\u03BCm)") +
      guides(color = guide_legend(title = paste0(maskName, "\noutline"))) +
      coord_fixed(ratio = 1) +
      facet_wrap(~ mask_files)
    
    cat(paste0("\nCommon Name:  ", common_name(names(mask_outline)), "\n"))
    cat(paste0("Plotting...Number of ", maskName, " outline files: ", length(mask_outline), "\n"))
    print(p)
  }
  
  # Rename the name of the mask according to the maskName parameter
  
  if(landmark){
    final_maskName <- "landmark"
  }else{
    final_maskName <- "mask"
  }
  
  mask <- lapply(mask, function(x){ 
    names(x) <- stringr::str_replace(names(x), "mask", final_maskName) 
    return(x)
  })
  
  mask_outline <- lapply(mask_outline, function(x){ 
    names(x) <- stringr::str_replace(names(x), "mask", final_maskName) 
    return(x)
  })
  
  
  # Order trackll and mask (Just in case)
  mask_filename <- gtools::mixedsort(names(mask))
  mask_outline_filename <- gtools::mixedsort(names(mask_outline))
  mask <- mask[mask_filename]
  mask_outline <- mask_outline[mask_outline_filename]
  
  
  # mask_filename <- gtools::mixedsort(gsub("_MASK.tif", "", names(mask)))
  
  if(landmark){
    return(list(landmark = mask, landmark_outline = mask_outline))
  }else{
    return(list(mask = mask, mask_outline = mask_outline))
  }
  
  # do.call("<<-",list(final_maskName, mask))
  # do.call("<<-",list(paste0(final_maskName, "_outline"), mask_outline))
  
}
mask_list <- readMask(folder, resolution = c(128, 128), pixel_size_um = 16/150, maskName = "MASK")
landmark_list <- readMask(folder, resolution = c(128, 128), pixel_size_um = 16/150, maskName = "NUCLEOLUS", landmark = TRUE)

# Match landmark with mask (The function finds mask and landmark in the environment)
matchMasknLandmark <- function(mask_list, landmark_list){

cat("Matching mask and landmark...\n")
  
landmark.area.split <- lapply(landmark_list[[1]], function(x) split(x, x$landmark))
landmark.outline.split <- lapply(landmark_list[[2]], function(x) split(x, x$landmark))

landmark.area.match <- list()
for(i in 1:length(mask_list[[1]])){
      
  landmark.area.match. <- lapply(landmark.area.split[[i]], function(x){
    
    match <- plyr::match_df(mask_list[[1]][[i]], x, on = c("x", "y"))
    x$landmark <- unique(match$mask)
    return(x)
  })
  landmark.area.match[[names(mask_list[[1]])[i]]] <- landmark.area.match.
}    

landmark <- track.merge(landmark.area.match, id = NULL)

landmark.outline.match <- list()
for(i in 1:length(mask_list[[1]])){
  
  landmark.outline.match. <- lapply(landmark.outline.split[[i]], function(x){
    
    match <- plyr::match_df(mask_list[[1]][[i]], x, on = c("x", "y"))
    x$landmark <- unique(match$mask)
    return(x)
  })
  landmark.outline.match[[names(mask_list[[1]])[i]]] <- landmark.outline.match.
}    


landmark_outline <- track.merge(landmark.outline.match, id = NULL)

landmark_list <- list(landmark = landmark, landmark_outline = landmark_outline)
return(landmark_list)
}
landmark_list <- matchMasknLandmark(mask_list = mask_list, landmark_list = landmark_list)

####################################### Step 3: Mask tracks based on nuclear mask #####################################################################################################################################################
track.masking <- function(trackll, mask_list, resolution = c(128, 128), pixel_size_um = 16/150, mask_track.plot = TRUE, plot_index = 1:length(trackll), maskName = "MASK", preprocess = FALSE){
    mask <- mask_list$mask
    mask_outline <- mask_list$mask_outline
  
    # Sort the movie filename in trackll and mask data
    trackll_filename <- gtools::mixedsort(names(trackll))
    mask_filename <- gtools::mixedsort(names(mask))
    
    # Order trackll and mask (Just in case)
    trackll <- trackll[trackll_filename]
    mask <- mask[mask_filename]
    
    # Check whether all the movies have the corresponding _MASK.tif
    if(!identical(trackll_filename, mask_filename)){
        cat("\nMovies do not match with _MASK.tif files!\n")
        stop()
    }
    
    if(mask_track.plot){
        # Merge mask data
        mask_outline_merged_list <- dplyr::bind_rows(mask_outline[plot_index], .id = "mask_files")
        mask_outline_merged_list$mask_files <- factor(mask_outline_merged_list$mask_files, levels = gtools::mixedsort(names(mask_outline)))

        # Merge trackll data
        trackll_merged_list <- track.merge(trackll[plot_index], id = "Trajectory")
        trackll_merged_merged_list <- bind_rows(trackll_merged_list, .id = "mask_files")
        trackll_merged_merged_list$mask_files <- factor(trackll_merged_merged_list$mask_files, levels = gtools::mixedsort(names(trackll)))

        # Covert pixel to um
        mask_outline_merged_list$x <-  mask_outline_merged_list$x * pixel_size_um
        mask_outline_merged_list$y <- mask_outline_merged_list$y * pixel_size_um
        mask_outline_merged_list$mask_centroid_x <- mask_outline_merged_list$mask_centroid_x * pixel_size_um
        mask_outline_merged_list$mask_centroid_y <- mask_outline_merged_list$mask_centroid_y * pixel_size_um
        trackll_merged_merged_list$x <- trackll_merged_merged_list$x * pixel_size_um
        trackll_merged_merged_list$y <- trackll_merged_merged_list$y * pixel_size_um
        
        # Plot mask and track on the same plot (y has a negative for flipping the plot vertically)
        p <- ggplot2::ggplot() +
            geom_point(data = trackll_merged_merged_list, aes(x,-y), color ="grey", alpha = 0.3, size = 0.1) +
            geom_path(data = mask_outline_merged_list, aes(x,-y, color = factor(mask))) +
            geom_point(data = mask_outline_merged_list, aes(mask_centroid_x, -mask_centroid_y), size = 2, color = "black") +
            xlim(0,resolution[1] * pixel_size_um) + ylim(-resolution[2] * pixel_size_um, 0) +
            xlab("x (\u03BCm)") + ylab("y (\u03BCm)") +
            guides(color = guide_legend(title = maskName)) +
            coord_fixed(ratio = 1) +
            facet_wrap(~ mask_files)

        cat(paste0("\nCommon Name:  ", common_name(trackll_filename), "\n"))
        cat(paste0("Number of trackll-mask pairs: ", length(trackll), "\n"))
        cat("Ploting mask and track overlay...\n")
        
        print(p)
    }    

if(preprocess){
  stop()
}        

cl <- parallel::makeCluster(detectCores()-1)
doSNOW::registerDoSNOW(cl)
iterations <- length(trackll)

cat("Masking tracks...\n")

pb <- pbmcapply::progressBar(max = iterations, style = "ETA")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Master for-loop
trackll_masked <- list()
trackll_masked <- foreach(i = 1:iterations,
                   .packages=c("tidyverse","plyr", "pbmcapply", "doSNOW"), 
                   .options.snow = opts) %dopar% {

# Match mask with track
trackll_masked_list <- lapply(trackll[[i]], function(x){
    match <- plyr::match_df(round(x), mask[[i]], on = c("x", "y"))
    
    # Remove track if it spans two masks
    match_mask <- plyr::match_df(mask[[i]], round(x), on = c("x", "y"))
    if(length(unique(match_mask$mask)) != 1){
        x <-list()
    }else{
    # Retain matched
    x <- x[x$Frame %in% match$Frame, ]
    
    # Label the trackll with mask, mask_centroid_x and mask_centroid_y
    x$mask <- rep(match_mask$mask[1], length(match$Frame))
    x$mask_centroid_x<- rep(match_mask$mask_centroid_x[1], length(match$Frame))
    x$mask_centroid_y<- rep(match_mask$mask_centroid_y[1], length(match$Frame)) 
    x$mask_radius_um<- rep(match_mask$mask_radius_um[1], length(match$Frame)) 
    return(x)
    }
 
})    
trackll_masked[[trackll_filename[i]]] <- trackll_masked_list
                   }

close(pb)
parallel::stopCluster(cl) 

names(trackll_masked) <- names(trackll)

return(trackll_masked)

}
trackll_masked <- track.masking(trackll_raw, mask_list = mask_list, resolution = c(128, 128), pixel_size_um = 16/150)

######################################## Step 4: Remove gaps generated by masking #####################################################################################################################################################
track.gapRemove <- function(trackll){

# Remove list with no element due to masking
trackll_masked <- lapply(trackll, function(x){
    lapply(x, function(y){
        if(purrr::is_empty(y)){
            y <- y
        }else if(nrow(y) == 0){
            y <- list()
        }else{
            y <- y
        } 
    }) 
})
trackll_masked <- lapply(trackll_masked, function(x) purrr::compact(x))

# Merge masked trackll, and add Trajectory column, and then split according to Trajectory
# track.merge is defined outside of this function
trackll_masked_merged <- track.merge(trackll_masked)
trackll_masked_temp <- lapply(trackll_masked_merged, function(x) split(x, x$Trajectory))

# Check which Trajectory has gap
gap <- lapply(trackll_masked_temp, function(x){
    lapply(x, function(y){
        seq(y[["Frame"]][1], y[["Frame"]][length(y[["Frame"]])]) != y[["Frame"]]
    })
})

# Set Trajectory with no gap to NULL
gap <- lapply(gap, function(x){
    lapply(x, function(y){
        if(sum(y) != 0){
            y <- y
        }else{
            y <- NULL
        }
    })
})
gap <- lapply(gap, function(x) purrr::compact(x))

# Trajectory with gaps
gap <- lapply(gap, function(x){names(x)})

# Remove gap
for(i in seq_along(trackll_masked_temp)){     # i: movie
    # In case no gaps in all the trackll
    if(purrr::is_empty(gap[[i]])){
        next
    }
    for (j in as.numeric(unlist(gap[i]))){        # j: track
        if (length(trackll_masked_temp[[i]][[j]]$Frame) == 1){
            next
        }else{
            for (k in 1: (length(trackll_masked_temp[[i]][[j]]$Frame)-1) ){
                if (trackll_masked_temp[[i]][[j]]$Frame[k]+1 != trackll_masked_temp[[i]][[j]]$Frame[k+1]){
                    trackll_masked_temp[[i]][[j]]$Trajectory[(k+1):length(trackll_masked_temp[[i]][[j]]$Frame)] <- 
                        trackll_masked_temp[[i]][[j]]$Trajectory[(k+1):length(trackll_masked_temp[[i]][[j]]$Frame)] + 0.0001
                }
            }
        }
    }
}

# Tidy the trackll_masked_gapRemoved data
trackll_masked_gapRemoved <- lapply(trackll_masked_temp, function(x) dplyr::bind_rows(x))
trackll_masked_gapRemoved <- lapply(trackll_masked_gapRemoved, function(x) split(x, x$Trajectory))
trackll_masked_gapRemoved <- lapply(trackll_masked_gapRemoved, function(x){
    lapply(x, function(y){ subset( y, select = -c(Trajectory) )
    })
})
trackll_masked_gapRemoved <- lapply(trackll_masked_gapRemoved, function(x){ 
    names(x) <- seq_along(x)
    return(x)
    })

return(trackll_masked_gapRemoved)

}
trackll_masked_gapRemoved <- track.gapRemove(trackll_masked)

#####################################  Step 5: Correlate tracks to individual cell  ###################################################################################################################################################
track.singleCell <- function(trackll_masked_gapRemoved, mask_list){
  
mask_outline <- mask_list$mask_outline
  
# Sort the movie filename in trackll and mask data
trackll_masked_gapRemoved_filename <- gtools::mixedsort(names(trackll_masked_gapRemoved))
mask_outline_filename <- gtools::mixedsort(names(mask_outline))

# Order trackll and mask (Just in case)
trackll_masked_gapRemoved <- trackll_masked_gapRemoved[trackll_masked_gapRemoved_filename]
mask_outline <- mask_outline[mask_outline_filename]

# Check whether all the movies have the corresponding _MASK.tif
if(!identical(trackll_masked_gapRemoved_filename, mask_outline_filename)){
    cat("\nTrackll do not match with mask_outline files!")
    stop()
}    
    
# Parallel proccessing by foreach  
cl <- parallel::makeCluster(detectCores()-1)
# centroid.disp.angle() is defined outside this function
parallel::clusterExport(cl, c("centroid.disp.angle")) # Allow foreach find the function centroid.disp.angle
doSNOW::registerDoSNOW(cl)
iterations <- length(trackll_masked_gapRemoved)

cat(paste0("\nCommon Name:  ", common_name(trackll_masked_gapRemoved_filename), "\n"))
cat(paste0("Number of trackll_masked_gapRemoved-mask_outline pairs: ", length(trackll_masked_gapRemoved_filename), "\n"))
cat("Working on single cell analysis...\n")

pb <- pbmcapply::progressBar(max = iterations, style = "ETA")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Master for-loop
trackll_singleCell <- list()
trackll_singleCell <- foreach(iterFile = 1:iterations,
                          .packages=c("tidyverse","dplyr", "pbmcapply", "doSNOW", "sf"), 
                          .options.snow = opts) %dopar% {
    
    # Add displacement and angle from the mask_centroid to trackll and mask_outline
    trackll_singleCell_list <- mapply(cbind, 
                                trackll_masked_gapRemoved[[iterFile]], 
                                centroid.disp.angle(input = trackll_masked_gapRemoved[[iterFile]]),
                                SIMPLIFY = F)
    
    mask_outline_singeCell <- mapply(cbind, 
                                     split(mask_outline[[iterFile]], mask_outline[[iterFile]]$mask), 
                                     centroid.disp.angle(input = split(mask_outline[[iterFile]], mask_outline[[iterFile]]$mask)),
                                     SIMPLIFY = F)
    mask_outline_singeCell <- dplyr::bind_rows(mask_outline_singeCell)
    
    # Find the intersection between the detection vector and the mask outline polygon
    max_disp_xy <- c()
    max_disp_xy <- lapply(trackll_singleCell_list, function(x){
        
        for(k in 1:nrow(x)){ # k: number of row in each track
            
        # Match mask_outline with the current track (for-loop: k)
        mask_matched <- mask_outline_singeCell$mask == x$mask[1]
        
        # Extend the detection vector to make a line for intersection
        x2 = x[k,]["mask_centroid_x"] + ((x[k,]["to_mask_centroid_disp"]+99999) * cos((x[k,]["to_mask_centroid_theta"])*pi/180))
        y2 = x[k,]["mask_centroid_y"] + ((x[k,]["to_mask_centroid_disp"]+99999) * sin((x[k,]["to_mask_centroid_theta"])*pi/180))
        
        # Calculate the intersection point >> maximum point
        line <- sf::st_linestring(as.matrix(rbind(unlist(x[k, ][, c("x", "y")]), unlist(c(x2, y2)))))
        polygon <- sf::st_polygon(list(as.matrix(mask_outline_singeCell[mask_matched, ][, c("x", "y")])))
        
        # If the detection is on the outline, then sf::st_intersection(polygon, line) return a vector of two components, instead of four
        # If the detection is "just" outside the outline, then sf::st_intersection(polygon, line) return a NULL vector; GEOMETRYCOLLECTION EMPTY
        max_x <- sf::st_intersection(polygon, line)[2] # Max. x touching the outline
        max_y <- sf::st_intersection(polygon, line)[4] # Max. y touching the outline

        if(length(sf::st_intersection(polygon, line)) != 4){
            max_x <- x[k,]$x
            max_y <- x[k,]$y
            to_mask_outline_max_disp <- x[k,]$to_mask_centroid_disp
        }else{
            to_mask_outline_max_disp <- sqrt( (max_x - x[k,]["mask_centroid_x"])^2 + (max_y - x[k,]["mask_centroid_y"])^2)
        }   
 
        # Save the data
        max_disp_xy. <- data.frame(to_mask_outline_max_disp, max_x, max_y)
        names(max_disp_xy.) <- c("to_mask_outline_max_disp", "max_x", "max_y")
        max_disp_xy <- rbind(max_disp_xy, max_disp_xy.)
        }
    
        return(max_disp_xy)
    })
        
    # Update trackll_singleCell_list with max_disp_xy information
    trackll_singleCell_list <- mapply(cbind, trackll_singleCell_list, max_disp_xy, SIMPLIFY = F)
    
    # Normalize x and y (assume the mask outline is a circle with radius = 1 AU)
    normalized_disp_xy <- lapply(trackll_singleCell_list, function(x){

        # Caclulate the normalized displacement: to_mask_centroid_disp.normalized
        to_mask_centroid_disp.normalized <- x$to_mask_centroid_disp/x$to_mask_outline_max_disp
        to_mask_centroid_disp.normalized[to_mask_centroid_disp.normalized > 1] <- 1
        
        # Calculate normalized x,y coordinate: x.normalized, y.normalized
        x.normalized <- x$mask_centroid_x + (to_mask_centroid_disp.normalized * cos((x$to_mask_centroid_theta)*pi/180))   
        y.normalized <- x$mask_centroid_y + (to_mask_centroid_disp.normalized * sin((x$to_mask_centroid_theta)*pi/180))  

        # Calculate normalized and centered x,y coordinate: x.normalized.centered, y.normalized.centered
        x.normalized.centered <- x.normalized - x$mask_centroid_x
        y.normalized.centered <- y.normalized - x$mask_centroid_y
        
        normalized_disp_xy <- cbind(to_mask_centroid_disp.normalized, x.normalized, y.normalized, x.normalized.centered, y.normalized.centered)
        colnames(normalized_disp_xy) <- c("to_mask_centroid_disp.normalized", "x.normalized", "y.normalized", "x.normalized.centered", "y.normalized.centered")
        return(normalized_disp_xy)
    })
    
    # Update trackll_singleCell_list with normalized_disp_xy information
    trackll_singleCell_list <- mapply(cbind, trackll_singleCell_list, normalized_disp_xy, SIMPLIFY = F)
    
    trackll_singleCell[[trackll_masked_gapRemoved_filename[iterFile]]] <- trackll_singleCell_list
} # End of master for-loop
                              
close(pb)
parallel::stopCluster(cl) 

# Need to rename the trackll after foreach
names(trackll_singleCell) <- names(trackll_masked_gapRemoved)

return(trackll_singleCell)
}
trackll_singleCell <- track.singleCell(trackll_masked_gapRemoved, mask_list = mask_list)

####################################   Step 6: Roatate and align all the tracks based on the landmark (nucleolus)   ###################################################################################################################
track.singleCell.landmark.rotate <- function(trackll_singleCell, mask_list, landmark_list, landmark_mask_track.plot = TRUE, plot_index = 1:length(trackll_singleCell), resolution = c(128, 128), pixel_size_um = 16/150){

  mask = mask_list$mask
  mask_outline = mask_list$mask_outline
  landmark_outline = landmark_list$landmark_outline
  
# Sort the movie filename in trackll_singleCell, mask_outline and landmark_outline
trackll_singleCell_filename <- gtools::mixedsort(names(trackll_singleCell))
mask_outline_filename <- gtools::mixedsort(names(mask_outline))
landmark_outline_filename <- gtools::mixedsort(names(landmark_outline))

# Order trackll_singleCell, mask_outline and landmark_outline (Just in case)
trackll_singleCell <- trackll_singleCell[trackll_singleCell_filename]
mask_outline <- mask_outline[mask_outline_filename]
landmark_outline <- landmark_outline[mask_outline_filename]

# Check whether all the movies have the corresponding _MASK.tif
if(!identical(trackll_singleCell_filename, mask_outline_filename) |
   !identical(trackll_singleCell_filename, landmark_outline_filename) |
   !identical(mask_outline_filename, landmark_outline_filename)
   ){
    cat("\trackll_singleCell, mask_outline, landmark_outline do not match!")
    stop()
}

if(landmark_mask_track.plot){
    # Merge landmark_outline data
    landmark_outline_merged_list <- dplyr::bind_rows(landmark_outline[plot_index], .id = "files")
    landmark_outline_merged_list$files <- factor(landmark_outline_merged_list$files, levels = gtools::mixedsort(names(landmark_outline)))
    landmark_outline_merged_list$landmark <- factor(as.numeric(landmark_outline_merged_list$landmark))
    
    # Merge mask_outline data
    mask_outline_merged_list <- dplyr::bind_rows(mask_outline[plot_index], .id = "files")
    mask_outline_merged_list$files <- factor(mask_outline_merged_list$files, levels = gtools::mixedsort(names(mask_outline)))
    mask_outline_merged_list$mask <- factor(as.numeric(mask_outline_merged_list$mask))
    
    # Merge trackll_singleCell data
    # track.merge() is defined outside the function
    trackll_singleCell_merged_list <- track.merge(trackll_singleCell[plot_index], id = "Trajectory")
    trackll_singleCell_merged_merged_list <- dplyr::bind_rows(trackll_singleCell_merged_list, .id = "files")
    trackll_singleCell_merged_merged_list$files <- factor(trackll_singleCell_merged_merged_list$files, levels = gtools::mixedsort(names(trackll_singleCell)))
    
    # Covert pixel to um
    mask_outline_merged_list$x <- mask_outline_merged_list$x * pixel_size_um
    mask_outline_merged_list$y <- mask_outline_merged_list$y * pixel_size_um
    mask_outline_merged_list$mask_centroid_x <- mask_outline_merged_list$mask_centroid_x * pixel_size_um
    mask_outline_merged_list$mask_centroid_y <- mask_outline_merged_list$mask_centroid_y * pixel_size_um
    
    landmark_outline_merged_list$x <- landmark_outline_merged_list$x * pixel_size_um
    landmark_outline_merged_list$y <- landmark_outline_merged_list$y * pixel_size_um
    landmark_outline_merged_list$landmark_centroid_x <- landmark_outline_merged_list$landmark_centroid_x * pixel_size_um
    landmark_outline_merged_list$landmark_centroid_y <- landmark_outline_merged_list$landmark_centroid_y * pixel_size_um
    
    trackll_singleCell_merged_merged_list$x <- trackll_singleCell_merged_merged_list$x * pixel_size_um
    trackll_singleCell_merged_merged_list$y <- trackll_singleCell_merged_merged_list$y * pixel_size_um

    # Plot landmark, mask and track on the same plot
    p <- ggplot2::ggplot() +
        geom_point(data = trackll_singleCell_merged_merged_list, aes(x,y), color ="grey", alpha = 0.3, size = 0.1) +
        
        geom_path(data = landmark_outline_merged_list, aes(x, y,  color = landmark)) +
        geom_point(data = landmark_outline_merged_list, aes(landmark_centroid_x, landmark_centroid_y), size = 1, color = "black", shape = 4) +
        
        geom_path(data = mask_outline_merged_list, aes(x,y, color = mask)) +
        geom_point(data = mask_outline_merged_list, aes(mask_centroid_x, mask_centroid_y), size = 2, color = "black") +
        
        xlim(0,resolution[1] * pixel_size_um) + ylim(0,resolution[2] * pixel_size_um) +
        xlab("x (\u03BCm)") + ylab("y (\u03BCm)") +
        guides(color = guide_legend(title = "Mask + Landmark")) +
        coord_fixed(ratio = 1) +
        facet_wrap(~ files)
    
    # common_name() is defined outside this function
    cat(paste0("\nCommon Name:  ", common_name(trackll_singleCell_filename), "\n"))
    cat(paste0("Number of trackll_singleCell-mask-landmark pairs: ", length(trackll_singleCell_filename), "\n"))
    cat("Ploting mask, landmark and track overlay...\n")
    
    print(p)
}    

# Construct landmark_mask file which match landmark with the corresponding mask
landmark_mask <- list()
for (i in seq_along(landmark_outline)){

    # Prepare a summary for match the landmark with the mask
    landmark_summary <- landmark_outline[[i]] %>% 
        dplyr::group_by(landmark) %>% 
        dplyr::summarize(mean(landmark_centroid_x), mean(landmark_centroid_y), .groups = 'drop')

    names(landmark_summary) <- c("landmark","x","y")
    landmark_summary <- split(landmark_summary, landmark_summary$landmark)
    
    
    # Match landmark with mask
    landmark_mask_pairs <- lapply(landmark_summary, function(x){
        match_mask <- plyr::match_df(mask[[i]], round(x), on = c("x", "y"))

            # Label the landmark_summary with mask, mask_centroid_x and mask_centroid_y
            x$mask <- match_mask$mask[1]
            x$mask_centroid_x<- match_mask$mask_centroid_x[1]
            x$mask_centroid_y<- match_mask$mask_centroid_y[1]
            return(x)

    })  
    # # Calculate landmark_centroid_to_mask_centroid_disp and landmark_centroid_to_mask_centroid_theta
    # landmark_mask_list <- mapply(cbind, landmark_mask_pairs, centroid.disp.angle(landmark_mask_pairs), SIMPLIFY = F)
    
    landmark_mask[[names(landmark_outline)[i]]] <- landmark_mask_pairs
}    

# This is a slight variation of trackll.singleCell function, but use specially for landmark_mask file
landmark.singleCell <- function(landmark_mask, mask_outline){
    
    # Sort the movie filename in trackll and mask data
    landmark_mask_filename <- gtools::mixedsort(names(landmark_mask))
    mask_outline_filename <- gtools::mixedsort(names(mask_outline))
    
    # Order trackll and mask (Just in case)
    landmark_mask <- landmark_mask[landmark_mask_filename]
    mask_outline <- mask_outline[mask_outline_filename]
    
    # Check whether all the movies have the corresponding _MASK.tif
    if(!identical(landmark_mask_filename, mask_outline_filename)){
        cat("\nTrackll do not match with mask_outline files!")
        stop()
    }    
    
    #####    
    cl <- parallel::makeCluster(detectCores()-1)
    parallel::clusterExport(cl, c("centroid.disp.angle")) # Allow foreach find the function centroid.disp.angle
    doSNOW::registerDoSNOW(cl)
    iterations <- length(landmark_mask)
    
    # cat(paste0("\nCommon Name:  ", common_name(landmark_mask_filename), "\n"))
    # cat(paste0("Number of landmark_mask-mask_outline pairs: ", length(landmark_mask_filename), "\n"))
    # cat("Working on single cell analysis...\n")
    # 
    # pb <- pbmcapply::progressBar(max = iterations, style = "ETA")
    # progress <- function(n) setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)
    
    # Master for-loop
    trackll_singleCell <- list()
    trackll_singleCell <- foreach(iterFile = 1:iterations,
                                  .packages=c("tidyverse","dplyr", "pbmcapply", "doSNOW", "sf", "data.table") 
                                  ) %dopar% {
                                      
                                      # Add displacement and angle from the mask_centroid to trackll and mask_outline
                                      trackll_singleCell_list <- mapply(cbind, 
                                                                        landmark_mask[[iterFile]], 
                                                                        centroid.disp.angle(input = landmark_mask[[iterFile]]),
                                                                        SIMPLIFY = F)
                                      
                                      mask_outline_singeCell <- mapply(cbind, 
                                                                       split(mask_outline[[iterFile]], mask_outline[[iterFile]]$mask), 
                                                                       centroid.disp.angle(input = split(mask_outline[[iterFile]], mask_outline[[iterFile]]$mask)),
                                                                       SIMPLIFY = F)
                                      mask_outline_singeCell <- dplyr::bind_rows(mask_outline_singeCell)
                                      
                                      # Find the intersection between the detection vector and the mask outline polygon
                                      max_disp_xy <- c()
                                      max_disp_xy <- lapply(trackll_singleCell_list, function(x){
                                          
                                          for(k in 1:nrow(x)){ # k: number of row in each track
                                              
                                              # Match mask_outline with the current track (for-loop: k)
                                              mask_matched <- mask_outline_singeCell$mask == x$mask[1]
                                              
                                              # Extend the detection vector to make a line for intersection
                                              x2 = x[k,]["mask_centroid_x"] + ((x[k,]["to_mask_centroid_disp"]+99999) * cos((x[k,]["to_mask_centroid_theta"])*pi/180))
                                              y2 = x[k,]["mask_centroid_y"] + ((x[k,]["to_mask_centroid_disp"]+99999) * sin((x[k,]["to_mask_centroid_theta"])*pi/180))
                                              
                                              # Calculate the intersection point >> maximum point
                                              line <- sf::st_linestring(as.matrix(rbind(unlist(x[k, ][, c("x", "y")]), unlist(c(x2, y2)))))
                                              polygon <- sf::st_polygon(list(as.matrix(mask_outline_singeCell[mask_matched, ][, c("x", "y")])))
                                              
                                              max_x <- sf::st_intersection(polygon, line)[2] # Max. x touching the outline
                                              max_y <- sf::st_intersection(polygon, line)[4] # Max. y touching the outline
                                              
                                              # If the detection is already on the outline
                                              if(is.null(max_x[[1]]) | is.null(max_y[[1]])){
                                                  max_x <- x[k,]$x
                                                  max_y <- x[k,]$y
                                                  to_mask_outline_max_disp <- x[k,]$to_mask_centroid_disp
                                              }else{
                                                  to_mask_outline_max_disp <- sqrt( (max_x - x[k,]["mask_centroid_x"])^2 + (max_y - x[k,]["mask_centroid_y"])^2)
                                              }   
                                              
                                              # Save the data
                                              max_disp_xy. <- data.frame(to_mask_outline_max_disp, max_x, max_y)
                                              names(max_disp_xy.) <- c("to_mask_outline_max_disp", "max_x", "max_y")
                                              max_disp_xy <- rbind(max_disp_xy, max_disp_xy.)
                                          }
                                          
                                          return(max_disp_xy)
                                      })
                                      
                                      # Update trackll_singleCell_list with max_disp_xy information
                                      trackll_singleCell_list <- mapply(cbind, trackll_singleCell_list, max_disp_xy, SIMPLIFY = F)
                                      
                                      
                                      
                                      # Normalize x and y (assume the mask outline is a circle with radius = 1 AU)
                                      normalized_disp_xy <- lapply(trackll_singleCell_list, function(x){
                                          
                                          # Caclulate the normalized displacement: to_mask_centroid_disp.normalized
                                          to_mask_centroid_disp.normalized <- x$to_mask_centroid_disp/x$to_mask_outline_max_disp
                                          to_mask_centroid_disp.normalized[to_mask_centroid_disp.normalized > 1] <- 1
                                          
                                          # Calculate normalized x,y coordinate: x.normalized, y.normalized
                                          x.normalized <- x$mask_centroid_x + (to_mask_centroid_disp.normalized * cos((x$to_mask_centroid_theta)*pi/180))   
                                          y.normalized <- x$mask_centroid_y + (to_mask_centroid_disp.normalized * sin((x$to_mask_centroid_theta)*pi/180))  
                                          
                                          # Calculate normalized and centered x,y coordinate: x.normalized.centered, y.normalized.centered
                                          x.normalized.centered <- x.normalized - x$mask_centroid_x
                                          y.normalized.centered <- y.normalized - x$mask_centroid_y
                                          
                                          normalized_disp_xy <- cbind(to_mask_centroid_disp.normalized, x.normalized, y.normalized, x.normalized.centered, y.normalized.centered)
                                          colnames(normalized_disp_xy) <- c("to_mask_centroid_disp.normalized", "x.normalized", "y.normalized", "x.normalized.centered", "y.normalized.centered")
                                          return(normalized_disp_xy)
                                      })
                                      
                                      
                                      # Update trackll_singleCell_list with normalized_disp_xy information
                                      trackll_singleCell_list <- mapply(cbind, trackll_singleCell_list, normalized_disp_xy, SIMPLIFY = F)
                                      
                                      # Rename for landmark file
                                      trackll_singleCell_list <- lapply(trackll_singleCell_list, function(x){
                                      data.table::setnames(x, 
                                                           old = c("x", "y", "to_mask_centroid_disp", "to_mask_centroid_theta", "to_mask_outline_max_disp", "max_x", "max_y", 
                                                                   "to_mask_centroid_disp.normalized", "x.normalized", "y.normalized", "x.normalized.centered", "y.normalized.centered"), 
                                                           new = c("landmark_centroid_x", "landmark_centroid_y", "landmark_centroid_to_mask_centroid_disp", "landmark_centroid_to_mask_centroid_theta", "landmark_centroid_to_mask_outline_max_disp", "landmark_centroid_max_x", "landmark_centroid_max_y", 
                                                                   "landmark_centroid_to_mask_centroid_disp.normalized", "landmark_centroid_x.normalized", "landmark_centroid_y.normalized", "landmark_centroid_x.normalized.centered", "landmark_centroid_y.normalized.centered")
                                      )
                                          return(x)
                                      })
                                      
                                      trackll_singleCell[[landmark_mask_filename[iterFile]]] <- trackll_singleCell_list
                                  } # End of master for-loop
    
    # close(pb)
    parallel::stopCluster(cl) 
    
    names(trackll_singleCell) <- names(landmark_mask)
    return(trackll_singleCell)
    
}
landmark_singleCell <- landmark.singleCell(landmark_mask, mask_outline)
landmark_singleCell_merged <- track.merge(landmark_singleCell, id = NULL)

# Master for loop
cl <- parallel::makeCluster(detectCores()-1)
doSNOW::registerDoSNOW(cl)
iterations <- length(trackll_singleCell)

cat("Rotating all detection based on the centroid of the landmark...\n")

pb <- pbmcapply::progressBar(max = iterations, style = "ETA")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

trackll_singleCell_landmarked_rotated <- list()
trackll_singleCell_landmarked_rotated <- foreach(iterFile = 1:iterations,
                              .packages=c("tidyverse","dplyr", "pbmcapply", "doSNOW", "sf", "SpatialGraph", "data.table"), 
                              .options.snow = opts) %dopar% {

    
    # Match landmark with trackll_singleCell
    trackll_singleCell_landmarked <- lapply(trackll_singleCell[[iterFile]], function(x){
        y <- c()
        for(i in 1:nrow(x)){
            y. <- dplyr::left_join(x[i, ],
                                       landmark_singleCell_merged[[iterFile]][x[i, ]$mask == landmark_singleCell_merged[[iterFile]]$mask, ],
                                       by = c("mask", "mask_centroid_x", "mask_centroid_y"))
         y <- rbind(y, y.)   
        }
        return(y)
    }) 
        
    # Rotate the detection. Rotation angle is 180 - theta of the landmark
    trackll_singleCell_landmarked_rotated_list <- lapply(trackll_singleCell_landmarked, function(x){
            x$x.normalized.centered.rotated <-  SpatialGraph::rotation(as.matrix(x[, c("x.normalized.centered", "y.normalized.centered")]), (180 - x$landmark_centroid_to_mask_centroid_theta[1])*pi/180)[, 1]
            x$y.normalized.centered.rotated <- SpatialGraph::rotation(as.matrix(x[, c("x.normalized.centered", "y.normalized.centered")]), (180 - x$landmark_centroid_to_mask_centroid_theta[1])*pi/180)[, 2]
            x$landmark_centroid_x.normalized.centered.rotated <- SpatialGraph::rotation(as.matrix(x[, c("landmark_centroid_x.normalized.centered", "landmark_centroid_y.normalized.centered")]), (180 - x$landmark_centroid_to_mask_centroid_theta[1])*pi/180)[, 1]
            x$landmark_centroid_y.normalized.centered.rotated <- SpatialGraph::rotation(as.matrix(x[, c("landmark_centroid_x.normalized.centered", "landmark_centroid_y.normalized.centered")]), (180 - x$landmark_centroid_to_mask_centroid_theta[1])*pi/180)[, 2]
            x$x.normalized.centered.rotated.mid <- mean(x$x.normalized.centered.rotated)
            x$y.normalized.centered.rotated.mid <- mean(x$y.normalized.centered.rotated)
            
            return(x)
    })

    trackll_singleCell_landmarked_rotated[[trackll_singleCell_filename[iterFile]]] <- trackll_singleCell_landmarked_rotated_list
    
}# End of master for loop
close(pb)
parallel::stopCluster(cl) 

names(trackll_singleCell_landmarked_rotated) <- names(trackll_singleCell)

return(trackll_singleCell_landmarked_rotated)
}
trackll_singleCell_landmarked_rotated <- track.singleCell.landmark.rotate(trackll_singleCell, mask_list = mask_list, landmark_list = landmark_list, resolution = c(128, 128), pixel_size_um = 16/150)

#########################################   Step 7: Merge data from individual movies into a single dataset  ##########################################################################################################################
files.merge <- function(trackll){
    
# Sort the movie filename in trackll
trackll_filename <- gtools::mixedsort(names(trackll))

# Order trackll (Just in case)
trackll <- trackll[trackll_filename]

# Merge track
    a <- track.merge(trackll)
    b <- bind_rows(a, .id = "filename")
    b <- b %>% 
        dplyr::select(-filename, filename)
    b$fileNum <- match(b$filename, trackll_filename)
    b$trackID <- paste0(b$fileNum , "_", b$mask, "_", b$Trajectory)
    c <- split(b, b$trackID)
    c <- c[gtools::mixedsort(names(c))]
    names(c) <- seq_along(names(c))

d <- lapply(c, function(x){
    x <- x %>% dplyr::select(-Trajectory)
    return(x)
})

trackll <- list(d)
names(trackll) <- common_name(trackll_filename)

return(trackll)
}
trackll_spatial <- files.merge(trackll_singleCell_landmarked_rotated)

###########################################   Step 8: 2-state HMM  ####################################################################################################################################################################
# Filter track with only 1 frame
track.filter <- function(trackll, minFrame = 2, maxFrame = Inf){
  
  Frame.filter <- lapply(trackll, function(x){
    lapply(x, function(y){
      if (nrow(y) < minFrame | nrow(y) > maxFrame){
        y <- NULL
      }else{
        y <- y
      }
    })
  })
  
  track.filtered <- lapply(Frame.filter, function(x){
    x[lengths(x) != 0]   
  })
  
  return(track.filtered)
}
trackll_spatial.filtered <- track.filter(trackll_spatial)

# Export trackll to vbSPT folder
trackll_export4vbSPT <- track.merge(trackll_spatial.filtered)[[1]]
vbSPT_filename <- "vbSPT_trackll"
write.csv(trackll_export4vbSPT[, c("Trajectory","Frame", "x","y","z")],"vbSPT path/Input/vbSPT_trackll.csv", row.names = TRUE)

# Run vbSPT using Matlab
matlabr::get_matlab(try_defaults = TRUE, desktop = FALSE, splash = FALSE, display = FALSE, wait = TRUE)
system("matlab -wait -nodesktop -nosplash  -r \"run('vbSPT path\\Step_1_Process_trackll.m'); exit\"")
system("matlab -wait -nodesktop -nosplash  -r \"run('vbSPT path\\Step_2_vbSPT_classify.m'); exit\"")

# Read vbSPT files
vbSPT_file <- R.matlab::readMat(file.path("vbSPT path", "vbSPT_classified", paste0(vbSPT_filename, "_classified.mat")))
states <- lapply(vbSPT_file[["CellTrackViterbiClass"]], function(x) as.numeric(c(unlist(x)[1], unlist(x))))

transit <- lapply(states, function(x){
  if (sum(x) == x[1]*length(x)){
    x <- rep(0, length(x))
  }else{
    x <- rep(1, length(x))
  }
})

HMM. <- mapply(cbind, split(trackll_export4vbSPT, trackll_export4vbSPT$Trajectory), states, transit, SIMPLIFY=FALSE)
HMM <- lapply(HMM., setNames, c(colnames(trackll_export4vbSPT), "states", "transit"))
trackll_spatial_HMM <- list(HMM)
names(trackll_spatial_HMM) <- names(trackll_spatial.filtered)

############################################# Step 9 (Optional): Split transitioning track after HMM states analysis  ################################################################################################################
trackll.HMM.splitTrack <- function(trackll_HMM){
  
  track_HMM <- track.merge(trackll_HMM)
  track_HMM <- bind_rows(track_HMM)
  track_HMM$Trajectory <- as.numeric(track_HMM$Trajectory)
  
  # This is the input data, but add a "Trajectory" column
  trackll_HMM. <- split(track_HMM, track_HMM$Trajectory)
  
  # Add a tag on transitioning track
  track_HMM_splitTrack <- lapply(trackll_HMM., function(x){
    if( length(unique(x$states)) == 1 ){
      return(x)
    }else{
      for(i in 1:(nrow(x)-1) ){
        if (x$states[i] != x$states[i+1]){
          x$Trajectory[(i+1) :nrow(x)] <-  x$Trajectory[(i+1) :nrow(x)] + 0.000001  # 0.000001 is just a marker for splitting track
          return(x)
        }
        
        
      }
    }
    
  })
  
  # Merge
  track_HMM_splitTrack <- dplyr::bind_rows(track_HMM_splitTrack)
  
  # Split while drop the Trajectory column, and then convert track to trackll
  track_HMM_splitTrack <- split( track_HMM_splitTrack[-which(names(track_HMM_splitTrack)=="Trajectory")], track_HMM_splitTrack$Trajectory )
  names(track_HMM_splitTrack) <- seq_along(names(track_HMM_splitTrack))
  trackll_HMM_splitTrack <- list(track_HMM_splitTrack)
  names(trackll_HMM_splitTrack) <- names(trackll_HMM)
  
  return(trackll_HMM_splitTrack)
}
trackll_spatial_HMM_splitTrack <- trackll.HMM.splitTrack(trackll_spatial_HMM[1])

################################################## Step 10: Temporal analysis ########################################################################################################################################################
# Standard parameters:
# pixel_size_um = 16/150   # Camera pixel size (um) / magnification
# dt = 30   # Total number of time interval analyzed (just set a large number)
# frameRate_sec = 10/1000
# DispFeature_dt_atMost = 1  # Min, Max and Mean of displacement of more than 1 dt actually does not make sense

# ----------Dcoef iterFit start from low to high dt----------
# DcoefFit_dt_atLeast = 6
# DcoefFit_dt_atMost = 10   # End point of iterative fitting
# DcoefFitRsqLimit = 0.8
# locError_um = 30/1000    # One of the anchor point for Dcoef iterative fitting

# Short time diffusion coefficient is used to calculated Boundedness and Trappedness.
# Current setting fit Dshort from dt2 to dt6 with iterative fitting

# ----------Alpha, Rc and Drift iterFit start from high to low dt----------
# AlphaFit_dt_atLeast = 6  # End point of iterative fitting
# AlphaFit_dt_atMost = 20  # Start point of fitting (if data has less dt than AlphaFit_dt_atMost, then follow the data)
# AlphaFitRsqLimit = 0.8

# RcFit_dt_atLeast = 6     # End point of iterative fitting
# RcFit_dt_atMost = 20     # Start point of fitting (if data has less dt than Rc_dt_atMost, then follow the data)
# Rc_max = 1               # Set it as the estimated radius of the nucleus (um))

# DriftFit_dt_atLeast = 6   # End point of iterative fitting
# DriftFit_dt_atMost = 20   # Start point of fitting (if data has less dt than Rc_dt_atMost, then follow the data)

# AngleBiasWeight = 25      # Change weight for extreme angle
# AngleBias_dt_atLeast = 1  # If set to 10, if the data does not have 10 dt, return NA. Note: 2 displacements make 1 angle

# pixel_size_um = 16/150; dt = 30; frameRate_sec = 10/1000;
# DcoefFit_dt_atLeast = 4; DcoefFit_dt_atMost = 10; DcoefFitRsqLimit = 0.8; locError_um = 30/1000;
# AlphaFit_dt_atLeast = 4; AlphaFit_dt_atMost = 20; AlphaFitRsqLimit = 0.8;
# RcFit_dt_atLeast = 4; RcFit_dt_atMost = 20; Rc_max = 1;
# DriftFit_dt_atLeast = 4; DriftFit_dt_atMost = 20;
# AngleBias_dt_atLeast = 1; AngleBiasWeight = 25

track.analyze <- function(trackll,
                          pixel_size_um = 16/150, dt = 30, frameRate_sec = 10/1000,
                          DcoefFit_dt_atLeast = 4, DcoefFit_dt_atMost = 10, DcoefFitRsqLimit = 0.8, locError_um = 30/1000, 
                          AlphaFit_dt_atLeast = 4, AlphaFit_dt_atMost = 20, AlphaFitRsqLimit = 0.8, 
                          RcFit_dt_atLeast = 4, RcFit_dt_atMost = 20, Rc_max = 1,
                          DriftFit_dt_atLeast = 4, DriftFit_dt_atMost = 20,
                          AngleBias_dt_atLeast = 1, AngleBiasWeight = 25){
  
  
  # Sort the movie filename in trackll
  trackll_filename <- gtools::mixedsort(names(trackll))
  
  # Order trackll (Just in case)
  trackll <- trackll[trackll_filename]
  
  # Master for loop (In order for using this function to caluclate trackll of mulitple factors, cannot use foreach here)
  for (iterFile in seq_along(trackll)){
    
    cat( paste0("\nWorking on ", trackll_filename[iterFile], "...Total of ", iterFile, "/", length(trackll_filename), " files\n"))    
    trackll_analyzed <- list()
    
    # Check if the trackll have the required columns
    colCheck <- bind_rows(trackll[[iterFile]])
    # Check if there are missing columns due to people using the very basic trackll data format with just Frame, x, y and z
    if(!any(names(colCheck) ==  "trackID") & !any(names(colCheck) ==  "states" & !any(names(colCheck) ==  "mask_radius_um"))){
      trackll_basic <- lapply(trackll[[iterFile]], function(x){
        x <- data.frame(Frame = x$Frame, x = x$x, y = x$y, z = x$z)
        return(x)
      })
    }
    # Check if the trackll is processed by trackll.HMM()
    if(!any(names(colCheck) ==  "states") & any(names(colCheck) ==  "trackID") & any(names(colCheck) ==  "mask_radius_um")){
      trackll_basic <- lapply(trackll[[iterFile]], function(x){
        x <- data.frame(Frame = x$Frame, x = x$x, y = x$y, z = x$z, trackID = x$trackID, mask_radius_um = x$mask_radius_um)
        return(x)
      })
    }
    # Only retain Frame, x, y, z, states, transit, trackID and mask_radius_um column for faster manipulation
    if(any(names(colCheck) ==  "states") & any(names(colCheck) ==  "trackID") & any(names(colCheck) ==  "mask_radius_um")){
      trackll_basic <- lapply(trackll[[iterFile]], function(x){
        x <- data.frame(Frame = x$Frame, x = x$x, y = x$y, z = x$z, trackID = x$trackID, states = x$states, transit = x$transit, mask_radius_um = x$mask_radius_um)
        return(x)
      })
    }
    
    # Detect if gaps in the track, if yes, the whole function will stop
    track.gapDetect <- function(trackll){    
      
      # Detect if there are gaps in the track  
      gap <- lapply(trackll, function(x){
        seq(x[["Frame"]][1], x[["Frame"]][length(x[["Frame"]])]) != x[["Frame"]]
      })
      
      if (sum(unlist(gap)) > 0){    
        stop("There are gaps in the following tracks, please remove gaps before proceeding.
            \nTracks with gaps (May not show all if the dataset is too large): ",  
             paste0(as.character(which(lapply(gap, function(x) sum(x) > 0) == T)), sep = ", "))
      }    
    }
    track.gapDetect(trackll_basic)
    
    # Calculate Track length
    cat("  Calculating Track length...\n")
    getTrackLength <- function(trackll, frameRate_sec = 10/1000){
      TrackLength <- lapply(trackll, function(x){
        TrackLength <- nrow(x) - 1
        DwellTime <- TrackLength * frameRate_sec
        
        data.frame(TrackLength, DwellTime)
        
      })
      
      return(TrackLength)
    }
    TrackLength <- getTrackLength(trackll_basic, frameRate_sec = frameRate_sec)
  
    
    
    # Calculate displacement for each lagtime (dt)
    cat("  Calculating Displacement...\n")
    getDisplacement <- function(trackll, pixel_size_um = 16/150, dt = 30){
      
      # Calculate displacement for different lagtimes
      displacement <-pbapply::pblapply(trackll, function(x) {
        
        displacement <- c()
        for (i in 1:dt){
          v <- x[c("x", "y")] - dplyr::lag(x[c("x", "y")], i)
          d <- sqrt(v$x^2 + v$y^2)*pixel_size_um
          
          displacement <- cbind(displacement, d)
        }    
        colnames(displacement) <- sapply(seq(1:dt), function(x) paste0("disp_dt", x))  
        return(displacement)
        
      })   
      trackll_basic_disp <- mapply(cbind, trackll, displacement, SIMPLIFY = FALSE)
      return(trackll_basic_disp)
    }
    trackll_basic_disp <- getDisplacement(trackll_basic, pixel_size_um = pixel_size_um, dt = dt)
    
    # Calculate single-step Dcoef
    # first: Fill the first point with SSDcoef of the second point (In theory, the first point has no SSDcoef)
    cat("  Calculating single-step Dcoef...\n")
    getSSDcoef <- function(trackll_basic_disp, frameRate_sec = 10/1000, first = TRUE){
      
      # Find column with name containing "disp_dt
      disp_dt_column <- which(!is.na(stringr::str_match(names(trackll_basic_disp[[1]]), "disp_dt1$"))) 
      
      # Calculate single-step Dcoef from dt1, no fitting
      SSDcoef <- pbapply::pblapply(trackll_basic_disp, function(x){
        
        SSDcoef <- cbind(x[, disp_dt_column]^2 / (4 * frameRate_sec), log10(x[, disp_dt_column]^2 / (4 * frameRate_sec)))
        
        # Fill the first point with SSDcoef of the second point (In theory, the first point has no SSDcoef)
        if(nrow(SSDcoef) > 1 & first){
          SSDcoef[1, ] <- SSDcoef[2, ]
        }
        
        colnames(SSDcoef) <- c("SSDcoef", "logSSDcoef") 
        return(SSDcoef)
        
        
      })   
      trackll_basic_disp <- mapply(cbind, trackll_basic_disp, SSDcoef, SIMPLIFY = FALSE)
      return(trackll_basic_disp)
      
    }
    trackll_basic_disp <- getSSDcoef(trackll_basic_disp, frameRate_sec = frameRate_sec, first = TRUE)
    
    # Calculate Displacement features
    cat("  Calculating Displacement features...\n")
    getDispFeature <- function(trackll_basic_disp, pixel_size_um = 16/150){
      
      # Find column with name containing "disp_dt
      # Always use 1 for DispFeature_dt_atMost !!! (Unless you have special purpose)
      disp_dt_column <- which(!is.na(stringr::str_match(names(trackll_basic_disp[[1]]), "disp_dt1$")))   #[1: DispFeature_dt_atMost]    
      
      # Calculate Displacement features
      DispFeature <- pbapply::pblapply(trackll_basic_disp, function(x){
        
        # max, min, mean and straightness are calculated for the first dt
        maxDisplacement <- max(unlist(x[, disp_dt_column]), na.rm = T)
        minDisplacement <- min(unlist(x[, disp_dt_column]), na.rm = T)
        meanDisplacement <- mean(unlist(x[, disp_dt_column]), na.rm = T)
        medianDisplacement <- median(unlist(x[, disp_dt_column]), na.rm = T)
        
        # Use the following lines if you want to calculate Displacement statistic for each dt                             
        # maxDisplacement <- lapply(trackll_basic_disp, function(x) apply(x[, disp_dt_column], 2, max, na.rm = TRUE))
        # minDisplacement <- lapply(trackll_basic_disp, function(x) apply(x[, disp_dt_column], 2, min, na.rm = TRUE))
        # avgDisplacement <- lapply(trackll_basic_disp, function(x) apply(x[, disp_dt_column], 2, mean, na.rm = TRUE))
        
        
        # Span = maxDisplacement for all dt
        Span <- max(dist(x[c("x", "y")])* pixel_size_um)
        
        # Straightness is also called effectiveness in the Diatrack manual
        Straightness <- 
          sqrt( (x$x[1] - x$x[nrow(x)])^2 + (x$y[1] - x$y[nrow(x)])^2 )*pixel_size_um / 
          sum(unlist(x[, disp_dt_column]), na.rm = T)
        
        # Ratio of squared end-to-end distance and the sum of squared distances
        # Based on Helmuth et al. (2007)
        Efficiency <- ( sqrt( (x$x[1] - x$x[nrow(x)])^2 + (x$y[1] - x$y[nrow(x)])^2 )*pixel_size_um )^2 / 
          sum( unlist(x[, disp_dt_column])^2, na.rm = T)
        
        # The fractal path dimension was implemented according Katz et al. (1985).
        # It takes between around 1 for straight paths, values around 2 for random paths, and values around 3 for constrained paths
        FractalDim <- log(nrow(x) - 1) / 
          (log(nrow(x) - 1) + log(Span / sum(unlist(x[, disp_dt_column]), na.rm = T)))
        
        
        
        #eigenvalues of radius of gyration tensor T(covariance matrix of x and y)
        if(nrow(x) == 1){
          Asymmetry1 <- NA
          Asymmetry2 <- NA
          Asymmetry3 <- NA
          Elongation <- NA
          AspectRatio <- NA
        }else{
          
          lambda1 <-  max(eigen(cov(cbind(x$x, x$y)))$values) #larger principal radius of gyration (eigenvalue)
          lambda2 <- min(eigen(cov(cbind(x$x, x$y)))$value) #smaller principal radius of gyration (eigenvalue)
          Asymmetry1 <- (lambda1 - lambda2)^2 / (lambda1 + lambda2)^2
          Asymmetry2 <- lambda2 / lambda1
          Asymmetry3 <- -log10( 1 - (lambda1 - lambda2)^2 / (2*(lambda1 + lambda2)^2))
          
          # minBox: Minimum bounding rectangle
          minBox <- shotGroups::getMinBBox(matrix(c(x$x,x$y), ncol = 2))
          Elongation <- 1 - min(c(minBox$width, minBox$height)) / max(c(minBox$width, minBox$height))
          AspectRatio <- max(c(minBox$width, minBox$height)) / min(c(minBox$width, minBox$height))
          
        }
        
        data.frame(maxDisplacement, minDisplacement, meanDisplacement, medianDisplacement, 
                   Span, Straightness, Efficiency, FractalDim,
                   Asymmetry1, Asymmetry2, Asymmetry3,
                   Elongation, AspectRatio)
      })
      
      return(DispFeature)
    }
    DispFeature <- getDispFeature(trackll_basic_disp, pixel_size_um = pixel_size_um)
    
    # Calculate MSD for each lagtime (dt)
    cat("  Calculating Mean sqaure displacement (MSD)...\n")
    getMSD <- function(trackll_basic_disp){
      
      # Find column with name containing "disp_dt
      disp_dt_column <- which(!is.na(stringr::str_match(names(trackll_basic_disp[[1]]), "disp_dt")))    
      
      # Calculate MSD
      MSD <- pbapply::pblapply(trackll_basic_disp, function(x) apply(x[, disp_dt_column]^2, 2, mean, na.rm = TRUE))
      
      # Name the MSD column
      MSD_dt_name <- sapply(seq(1:dt), function(x) paste0("MSD_dt", x))        
      MSD <- lapply(MSD, setNames, MSD_dt_name)
      MSD <- lapply(MSD, t)   # This is important for maintaining proper data structure!
      
      return(MSD)
    }
    MSD <- getMSD(trackll_basic_disp)
    
    # Gaussianity. Whereas normal diffusion shows values of 0, other motion types show deviations of zero.
    # It return values with multiple lag times, so MEDIAN Gaussianity is calculated instead
    cat("  Calculating Gaussianity...\n")
    getGaussianity <- function(trackll_basic){
      
      #Trajectory's quartic moment
      getQuarticMoment <- function(trackll, pixel_size_um = 16/150, dt = 30){
        
        # Calculate displacement for different lagtimes
        QuarticMoment <-pbapply::pblapply(trackll, function(x) {
          
          QuarticMoment <- c()
          for (i in 1:dt){
            v <- x[c("x", "y")] - dplyr::lag(x[c("x", "y")], i)
            q <- (v$x * pixel_size_um)^4 + (v$y * pixel_size_um)^4
            
            QuarticMoment <- cbind(QuarticMoment, q)
          }    
          colnames(QuarticMoment) <- sapply(seq(1:dt), function(x) paste0("quad_dt", x))  
          return(QuarticMoment)
          
        })   
        trackll_basic_QuarticMoment <- mapply(cbind, trackll, QuarticMoment, SIMPLIFY = FALSE)
        return(trackll_basic_QuarticMoment)
      }
      trackll_basic_QuarticMoment <- getQuarticMoment(trackll_basic, pixel_size_um = pixel_size_um, dt = dt)
      #Mean Trajectory's quartic moment
      getMQM <- function(trackll_basic_QuarticMoment){
        
        # Find column with name containing "quad_dt
        quad_dt_column <- which(!is.na(stringr::str_match(names(trackll_basic_QuarticMoment[[1]]), "quad_dt")))    
        
        # Calculate MQD
        MQM <- lapply(trackll_basic_QuarticMoment, function(x) apply(x[, quad_dt_column], 2, mean, na.rm = TRUE))
        
        # Name the MQD column
        MQM_dt_name <- sapply(seq(1:dt), function(x) paste0("MQM_dt", x))        
        MQM <- lapply(MQM, setNames, MQM_dt_name)
        MQM <- lapply(MQM, t)   # This is important for maintaining proper data structure!
        
        return(MQM)
      }
      MQM <- getMQM(trackll_basic_QuarticMoment)
      
      Gaussianity <- mapply(x = MQM, y = MSD, function(x, y) data.frame(median(((2*x) / (3*y*y)) - 1, na.rm = TRUE)) ,SIMPLIFY = FALSE)
      Gaussianity <- lapply(Gaussianity, function(x) {
        names(x) <- "Gaussianity"
        return(x)
      })
      
      return(Gaussianity)
    }
    Gaussianity <- getGaussianity(trackll_basic)
    
    # Calculate Diffusion coefficient, MSD = 4 * Dcoef * dt
    # iterFit, extend 1 dt per iteration, and then iterFit again with three anchor point (dt2, locError_um, 0)
    # According to my observation, achor point does not improve the number of fit
    cat("  Calculating Diffusion coefficient (Dcoef)...\n")
    get2D_Dcoef <- function(MSD, frameRate_sec = 10/1000, DcoefFit_dt_atMost = 20, DcoefFit_dt_atLeast = 6, DcoefFitRsqLimit = 0.8, locError_um = 30/1000){
      
      #MSD linear regression to fit Diffusion coefficient (Dcoef)
      Dcoef <- pbapply::pblapply(MSD, function(x){   #pblapply?
        
        if(sum(is.na(x[1:DcoefFit_dt_atLeast])) > 0){
          Dcoef <- Inf
          DcoefFit_intercept <- Inf
          DcoefFitRsq <- Inf
          DcoefFitSummary <- "None"
        }else{
          DcoefFit <- lm( x[1:DcoefFit_dt_atLeast] ~ c((1:DcoefFit_dt_atLeast) * frameRate_sec) )
          if(length(DcoefFit$coef) == 2){
            Dcoef <- DcoefFit$coef[2]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- DcoefFit$coef[1]
          }else{
            Dcoef <- DcoefFit$coef[1]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- 0
          }
          DcoefFitRsq <- summary(DcoefFit)$r.squared
          DcoefFitSummary <- "Good"
        }
        
        # Set poor fit value to -1 to avoid error in the while-loop
        if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
          Dcoef <- -1
          DcoefFit_intercept <- -1
          DcoefFitRsq <- -1
        }
        
        # Basic iterFit, without anchor
        j <- 1
        while( (Dcoef <= 0 | DcoefFitRsq < DcoefFitRsqLimit) &
               sum(!is.na(x)) > DcoefFit_dt_atLeast & j + DcoefFit_dt_atLeast <= DcoefFit_dt_atMost ){
          
          a <- x[1:(j+DcoefFit_dt_atLeast)] # Extend 1 dt for MSD
          b <- c(1:(j+DcoefFit_dt_atLeast)) * frameRate_sec  # Extend 1 dt for lagtime
          
          DcoefFit <- lm( a ~ b )
          if(length(DcoefFit$coef) == 2){
            Dcoef <- DcoefFit$coef[2]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- DcoefFit$coef[1]
          }else{
            Dcoef <- DcoefFit$coef[1]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- 0
          }
          
          DcoefFitRsq <- summary(DcoefFit)$r.squared
          DcoefFitSummary <- paste("Extended", j, "dt")
          
          if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
            Dcoef <- -1
            DcoefFit_intercept <- -1
            DcoefFitRsq <- -1
          }
          j <- j + 1
        }
        
        # iterFit, anchor to second dt
        j <- 0
        while( (Dcoef <= 0 | DcoefFitRsq < DcoefFitRsqLimit) &
               sum(!is.na(x)) > DcoefFit_dt_atLeast & j + DcoefFit_dt_atLeast <= DcoefFit_dt_atMost ){
          
          a <- x[1:(j+DcoefFit_dt_atLeast)] # Extend 1 dt for MSD
          b <- c(1:(j+DcoefFit_dt_atLeast)) * frameRate_sec  # Extend 1 dt for lagtime
          
          DcoefFit <- lm(I(a - a[2]) ~ I(b - b[2]) + 0) # Anchor to second dt
          if(length(DcoefFit$coef) == 2){
            Dcoef <- DcoefFit$coef[2]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- DcoefFit$coef[1]
          }else{
            Dcoef <- DcoefFit$coef[1]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- 0
          }
          DcoefFitRsq <- summary(DcoefFit)$r.squared
          DcoefFitSummary <- paste("Extended", j, "dt, and by anchoring to dt2")
          
          if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
            Dcoef <- -1
            DcoefFit_intercept <- -1
            DcoefFitRsq <- -1
          }
          j <- j + 1
        }
        
        # iterFit, anchor to 4 * localization error^2
        j <- 0
        while( (Dcoef <= 0 | DcoefFitRsq < DcoefFitRsqLimit) &
               sum(!is.na(x)) > DcoefFit_dt_atLeast & j + DcoefFit_dt_atLeast <= DcoefFit_dt_atMost ){
          
          a <- x[1:(j+DcoefFit_dt_atLeast)] # Extend 1 dt for MSD
          b <- c(1:(j+DcoefFit_dt_atLeast)) * frameRate_sec  # Extend 1 dt for lagtime
          
          DcoefFit <- lm(a ~ b + 0, offset = rep((locError_um)^2, j+DcoefFit_dt_atLeast))
          if(length(DcoefFit$coef) == 2){
            Dcoef <- DcoefFit$coef[2]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- DcoefFit$coef[1]
          }else{
            Dcoef <- DcoefFit$coef[1]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- 0
          }
          DcoefFitRsq <- summary(DcoefFit)$r.squared
          DcoefFitSummary <- paste("Extended", j, "dt, and by anchoring to 4 * localization error^2")
          
          if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
            Dcoef <- -1
            DcoefFit_intercept <- -1
            DcoefFitRsq <- -1
          }
          j <- j + 1
        }
        
        # iterFit, anchor to 0
        j <- 0
        while( (Dcoef <= 0 | DcoefFitRsq < DcoefFitRsqLimit) &
               sum(!is.na(x)) > DcoefFit_dt_atLeast & j + DcoefFit_dt_atLeast <= DcoefFit_dt_atMost ){
          
          a <- x[1:(j+DcoefFit_dt_atLeast)] # Extend 1 dt for MSD
          b <- c(1:(j+DcoefFit_dt_atLeast)) * frameRate_sec  # Extend 1 dt for lagtime
          
          DcoefFit <- lm(a ~ b + 0, offset = rep(0, j+DcoefFit_dt_atLeast))
          if(length(DcoefFit$coef) == 2){
            Dcoef <- DcoefFit$coef[2]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- DcoefFit$coef[1]
          }else{
            Dcoef <- DcoefFit$coef[1]/4 # MSD = 4 * Dcoef * dt
            DcoefFit_intercept <- 0
          }
          DcoefFitRsq <- summary(DcoefFit)$r.squared
          DcoefFitSummary <- paste("Extended", j, "dt, and by anchoring to 0")
          
          if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
            Dcoef <- -1
            DcoefFit_intercept <- -1
            DcoefFitRsq <- -1
          }
          j <- j + 1
        }
        
        # Set poor fit value to -1 to avoid error in the while-loop
        if(sum(is.na(c(Dcoef, DcoefFit_intercept, DcoefFitRsq))) > 0){
          Dcoef <- -1
          DcoefFit_intercept <- -1
          DcoefFitRsq <- -1
        }
        
        # If poor fit    
        if(Dcoef <= 0 | DcoefFitRsq < DcoefFitRsqLimit){
          Dcoef <- Inf
          DcoefFit_intercept <- Inf
          DcoefFitRsq <- Inf
          DcoefFitSummary <- "Fail to obatin good fit"
        }
        
        # Calculate logDcoef
        logDcoef <- log10(Dcoef)
        
        #Summarize
        data.frame(Dcoef, logDcoef, DcoefFit_intercept, DcoefFitRsq, DcoefFitSummary)
        
      })    
      return(Dcoef)
    }
    Dcoef <- get2D_Dcoef(MSD, frameRate_sec = frameRate_sec, DcoefFit_dt_atMost = DcoefFit_dt_atMost, DcoefFit_dt_atLeast = DcoefFit_dt_atLeast, DcoefFitRsqLimit = DcoefFitRsqLimit, locError_um = locError_um)
    
    cat("  Calculating Short time diffusion coefficient...\n")
    Dshort <- get2D_Dcoef(MSD, frameRate_sec = frameRate_sec, DcoefFit_dt_atMost = 6, DcoefFit_dt_atLeast = 2, DcoefFitRsqLimit = DcoefFitRsqLimit, locError_um = locError_um)
    
    cat("  Calculating Boundedness...\n")
    Boundedness <- pbapply::pbmapply(x = Dshort, y = DispFeature, z = TrackLength, function(x,y,z) { 
      Boundedness <- x$Dcoef * z$DwellTime / (y$Span / 2)^2
      if(is.na(Boundedness) | Boundedness == Inf){
        Boundedness <- NA
      }else{
        Boundedness <- Boundedness
      }
      names(Boundedness) <- "Boundedness"
      return(data.frame(Boundedness))
      
    }, SIMPLIFY = FALSE)
    cat("  Calculating Trappedness...\n")
    Trappedness <- pbapply::pblapply(Boundedness, function(x){ 
      
      Trappedness <- 1- exp(0.2048-2.5117*x)
      if(is.na(Trappedness)){
        Trappedness <- NA
      }else if (Trappedness < 0){
        Trappedness <- 0
      }else{
        Trappedness <- Trappedness
      }
      names(Trappedness) <- "Trappedness"
      return(data.frame(Trappedness))
      
    })
    cat("  Calculating MSDratio...\n")
    MSDratio <- pbapply::pblapply(MSD, function(x){
      
      MSDratio <- c()
      for(i in 1 : (sum(!is.na(x))-1) ){
        MSDratio. <- (x[i] / x[i + 1]) - ((i) / (i + 1))
        MSDratio <- c(MSDratio, MSDratio.)
      }
      MSDratio <-median(MSDratio)
      names(MSDratio) <- "MSDratio"
      return(data.frame(MSDratio))
      
      
    })
    
    
    # Fit for alpha, log(MSD) = Alpha * log(dt) + offset (# Add 1e-100 to MSD to avoid error after log?)
    # iterFit, decrease 1 dt per iteration, til AlphaFit_dt_atLeast, and then iterFit again with anchor point: dt2
    # According to my observation, achor point does not improve the number of fit
    cat("  Calculating Anomalous diffusion exponent (Alpha)...\n")
    getAlpha <- function(MSD, frameRate_sec = 10/1000, AlphaFit_dt_atMost = 20, AlphaFit_dt_atLeast = 6, AlphaFitRsqLimit = 0.8){
      
      #log(MSD) linear regression to fit Alpha, anomalous diffusion exponent
      Alpha <- pbapply::pblapply(MSD, function(x){   #pblapply?
        
        if(sum(!is.na(x)) < AlphaFit_dt_atLeast){
          Alpha <- Inf
          Alphafit_intercept <- Inf
          AlphaFitRsq <- Inf
          AlphaFitSummary <- "None"
        }else{
          Alphafit <- lm(log10(x[1: min(sum(!is.na(x)), AlphaFit_dt_atMost) ]) ~ log10(c(1: min(sum(!is.na(x)), AlphaFit_dt_atMost) ) * frameRate_sec))
          Alpha <- Alphafit$coef[2]
          Alphafit_intercept <- Alphafit$coef[1]
          AlphaFitRsq <- summary(Alphafit)$r.squared
          AlphaFitSummary <- "Good"
        }
        
        # Set poor fit value to -1 to avoid error in the while-loop
        if(sum(is.na(c(Alpha, Alphafit_intercept, AlphaFitRsq))) > 0){
          Alpha <- -1
          Alphafit_intercept <- -1
          AlphaFitRsq <- -1
        }
        
        # Basic iterFit, without anchor
        # 0.02 mean extreme value of alpha, which can be caused by the algo. fitting the last part of the curve
        j <- 1
        while( (Alpha <= 0.02 | AlphaFitRsq < AlphaFitRsqLimit) &
               (min(sum(!is.na(x)), AlphaFit_dt_atMost) - j) >= AlphaFit_dt_atLeast & (min(sum(!is.na(x)), AlphaFit_dt_atMost) - j) > 2 ){
          
          a <- log10( x[1:( min(sum(!is.na(x)), AlphaFit_dt_atMost) - j)] )
          b <- log10( c(1:( min(sum(!is.na(x)), AlphaFit_dt_atMost) - j)) * frameRate_sec )
          
          Alphafit <- lm(a ~ b)
          Alpha <- Alphafit$coef[2]
          Alphafit_intercept <- Alphafit$coef[1]
          AlphaFitRsq <- summary(Alphafit)$r.squared
          AlphaFitSummary <- paste("Decreased", j, "dt")
          
          if(sum(is.na(c(Alpha, Alphafit_intercept, AlphaFitRsq))) > 0){
            Alpha <- -1
            Alphafit_intercept <- -1
            AlphaFitRsq <- -1
          }
          j <- j + 1
        }
        
        # iterFit with anchor: dt2
        j <- 0
        while( (Alpha <= 0.02 | AlphaFitRsq < AlphaFitRsqLimit) &
               (min(sum(!is.na(x)), AlphaFit_dt_atMost) - j) >= AlphaFit_dt_atLeast & (min(sum(!is.na(x)), AlphaFit_dt_atMost) - j) > 2 ){
          
          a <- log10( x[1:( min(sum(!is.na(x)), AlphaFit_dt_atMost) - j)] )
          b <- log10( c(1:( min(sum(!is.na(x)), AlphaFit_dt_atMost) - j)) * frameRate_sec )
          
          Alphafit <- lm(I(a - a[2]) ~ I(a - a[2]) + 0)
          Alpha <- Alphafit$coef[2]
          Alphafit_intercept <- Alphafit$coef[1]
          AlphaFitRsq <- summary(Alphafit)$r.squared
          AlphaFitSummary <- paste("Decreased", j, "dt, and by anchoring to dt2")
          
          if(sum(is.na(c(Alpha, Alphafit_intercept, AlphaFitRsq))) > 0){
            Alpha <- -1
            Alphafit_intercept <- -1
            AlphaFitRsq <- -1
          }
          j <- j + 1
        }
        
        # Set poor fit value to -1 to avoid error
        if(sum(is.na(c(Alpha, Alphafit_intercept, AlphaFitRsq))) > 0){
          Alpha <- -1
          Alphafit_intercept <- -1
          AlphaFitRsq <- -1
        }
        
        # If poor fit    
        if(Alpha <= 0.02 | AlphaFitRsq < AlphaFitRsqLimit){
          Alpha <- Inf
          Alphafit_intercept <- Inf
          AlphaFitRsq <- Inf
          AlphaFitSummary <- "Fail to obatin good fit"
        }
        
        #Summarize
        data.frame(Alpha, Alphafit_intercept, AlphaFitRsq, AlphaFitSummary)
      })    
      
      return(Alpha)
      
    }
    Alpha <- getAlpha(MSD, frameRate_sec = frameRate_sec, AlphaFit_dt_atMost = AlphaFit_dt_atMost, AlphaFit_dt_atLeast = AlphaFit_dt_atLeast, AlphaFitRsqLimit = AlphaFitRsqLimit)
    
    # Combine MSD and Alpha (and radius of the nucleus) data for subdiffusion analysis
    if(any(names(colCheck) ==  "mask_radius_um")){
      
      mask_radius_um <- lapply(trackll_basic, function(x) {
        mask_radius_um <- x[1,"mask_radius_um"]
        names(mask_radius_um) <- "mask_radius_um"
        return(data.frame(mask_radius_um))
      })
      MSD_Alpha <- mapply(cbind, MSD, Alpha, mask_radius_um, SIMPLIFY = FALSE)
      
      
    }else{
      MSD_Alpha <- mapply(cbind, MSD, Alpha, SIMPLIFY = FALSE)
    }
    
    
    
    
    
    # Radius of confinement (This function take a while to run, maybe > 10 min)
    # iterFit, decrease 1 dt per iteration, til RcFit_dt_atLeast
    # mle2 works better than nls for Rc fitting
    cat("  Calculating Radius of confinement (Rc)...\n")
    getRc <- function(MSD_Alpha, frameRate_sec = 10/1000, RcFit_dt_atMost = 20, RcFit_dt_atLeast = 6, Rc_max = 1, locError_um = 30/1000){
      
      
      cl <- parallel::makeCluster(detectCores()-1)
      parallel::clusterExport(cl, c("MSD_Alpha", "frameRate_sec", "RcFit_dt_atMost", "RcFit_dt_atLeast", "Rc_max", "locError_um"), envir=environment())
      
      
      Rc <- pbapply::pblapply(MSD_Alpha, function(k){   #Use k as the index here, because x is used to define the nll function
        
        m <- k[which(!is.na(stringr::str_match(names(k), "MSD")))] # Subset the MSD_dt column
        
        if(any(colnames(k) == "mask_radius_um")){
          Rc_max <- k$mask_radius_um
        }else{
          Rc_max <- Rc_max
        }
        
        if((sum(!is.na(m)) < RcFit_dt_atLeast) | k$Alpha >= 1){
          Rc <- Inf
          RcDmicro <- Inf
          Rc_offset <- Inf
          Rc_Nuc_ratio <- Inf
          RcFitSummary <- "None"
        }else{
          nll <- function(Rc, RcDmicro, Rc_offset) { #nll: negative log likelihood function 
            x <- c(1:( min(sum(!is.na(m)), RcFit_dt_atMost)  )) * frameRate_sec
            y <- m[1:( min(sum(!is.na(m)), RcFit_dt_atMost) )]
            MSD = I(Rc^2 * (1 - exp((-4 * RcDmicro* x)/Rc^2)) + Rc_offset) # Confinement diffusion equation
            -sum(y*(log(MSD)) - MSD) #This is the nll function
          }
          RcFit. <- try(RcFit <- bbmle::mle2(nll, start = list(Rc = 0.05, RcDmicro = 0.03, Rc_offset = 0.001), control=list(maxit=5000))
                        ,silent = TRUE)
          
          if(class(RcFit.) == "try-error"){
            Rc <- -1
            RcDmicro <- -1
            Rc_offset <- -1
          }else{
            Rc <- RcFit@coef[1]
            RcDmicro <- RcFit@coef[2]
            Rc_offset <- RcFit@coef[3]
            Rc_Nuc_ratio <- Rc / Rc_max
            RcFitSummary <- "Good"
          }
        }
        
        # Set poor fit value to -1 to avoid error in the while-loop
        if(sum(is.na(c(Rc, RcDmicro, Rc_offset))) > 0){
          Rc <- -1
          RcDmicro <- -1
          Rc_offset <- -1
        }
        
        # Basic iterFit, without anchor. Set Rc to be larger than 0.03 um (localization error); For HHT1, bound H3 has a Rc of 0.095 um.
        j <- 1
        while( (Rc <= locError_um | Rc > Rc_max & Rc < Inf | RcDmicro <=0) & (min(sum(!is.na(m)), RcFit_dt_atMost) - j) >= RcFit_dt_atLeast & (min(sum(!is.na(m)), RcFit_dt_atMost) - j) > 2 ){
          
          nll <- function(Rc, RcDmicro, Rc_offset) { #nll: negative log likelihood function
            x <- c(1:( min(sum(!is.na(m)), RcFit_dt_atMost) - j)) * frameRate_sec
            y <- m[1:( min(sum(!is.na(m)), RcFit_dt_atMost) - j)]
            MSD = I(Rc^2 * (1 - exp((-4 * RcDmicro* x)/Rc^2)) + Rc_offset) # Confinement diffusion equation
            -sum(y*(log(MSD)) - MSD) #This is the nll function
          }
          RcFit. <- try(RcFit <- bbmle::mle2(nll, start = list(Rc = 0.05, RcDmicro = 0.03, Rc_offset = 0.001), control=list(maxit=5000))
                        ,silent = TRUE)
          
          if(class(RcFit.) == "try-error"){
            Rc <- -1
            RcDmicro <- -1
            Rc_offset <- -1
          }else{
            Rc <- RcFit@coef[1]
            RcDmicro <- RcFit@coef[2]
            Rc_offset <- RcFit@coef[3]
            Rc_Nuc_ratio <- Rc / Rc_max
            RcFitSummary <- paste("Decreased", j, "dt")
          }
          
          if(sum(is.na(c(Rc, RcDmicro, Rc_offset))) > 0){
            Rc <- -1
            RcDmicro <- -1
            Rc_offset <- -1
          }
          j <- j + 1
        }
        
        # Set poor fit value to -1 to avoid error
        if(sum(is.na(c(Rc, RcDmicro, Rc_offset))) > 0){
          Rc <- -1
          RcDmicro <- -1
          Rc_offset <- -1
        }
        
        # If poor fit
        if(Rc > Rc_max & Rc < Inf){
          Rc <- Inf
          RcDmicro <- Inf
          Rc_offset <- Inf
          Rc_Nuc_ratio <- Inf
          RcFitSummary <- "Fail to obatin good fit (Rc too large)"
        }
        
        if(Rc <= locError_um | RcDmicro <= 0){
          Rc <- NA
          RcDmicro <- NA
          Rc_offset <- NA
          Rc_Nuc_ratio <- NA
          RcFitSummary <- "Fail to obatin good fit (Rc too small or RcDmicro is negative)"
        }
        
        #Summarize
        return(data.frame(Rc, RcDmicro, Rc_offset, Rc_Nuc_ratio, RcFitSummary))
      }, cl = cl)    
      
      parallel::stopCluster(cl)
      
      return(Rc)
    }
    Rc <- getRc(MSD_Alpha, frameRate_sec = frameRate_sec, RcFit_dt_atMost = RcFit_dt_atMost, RcFit_dt_atLeast = RcFit_dt_atLeast, Rc_max = Rc_max)
    
    # Drift velocity, MSD = 4 * D *dt + (v * dt)^2
    # There is chance that D becomes zero, and the equation becomes MSD = (v * dt)^2
    # iterFit, decrease 1 dt per iteration, til DriftFit_dt_atLeast
    # For unknown reason, mle2 does not work for Drift fitting, need to use nls, although the fitting is not very accurate
    cat("  Calculating Drift velocity (Drift)...\n")
    getDrift <- function(MSD_Alpha, frameRate_sec = 10/1000, DriftFit_dt_atMost = 20, DriftFit_dt_atLeast = 6){
      
      
      Drift <- pbapply::pblapply(MSD_Alpha, function(k){   #Use k as the index here, becuase x is used to define the nll function
        
        # Subset the MSD_dt column
        m <- k[which(!is.na(stringr::str_match(names(k), "MSD")))]
        
        if((sum(!is.na(m)) < DriftFit_dt_atLeast) | k$Alpha <= 1){
          Drift <- Inf
          DriftDcoef <- Inf
          Drift_offset <- Inf
          DriftFitSummary <- "None"
        }else{
          
          
          x <- c(1:( min(sum(!is.na(m)), DriftFit_dt_atMost) )) * frameRate_sec
          y <- unlist(m[1:( min(sum(!is.na(m)), DriftFit_dt_atMost) )])
          
          DriftFit. <- try(
            DriftFit <- nls(y ~ 4 * DriftDcoef * x + (Drift * x)^2 + Drift_offset, trace = F,
                            start=c(Drift = 5, DriftDcoef = 1.5, Drift_offset = 0),
                            lower=c(0,0,0),
                            upper=c(50,50,0.5), nls.control(maxiter = 2000), algorithm = "port"), 
            silent = TRUE)
          
          if(class(DriftFit.) == "try-error"){
            Drift <- -1
            DriftDcoef <- -1
            Drift_offset <- -1
          }else{
            Drift <- coefficients(DriftFit)[1]
            DriftDcoef <- coefficients(DriftFit)[2]
            Drift_offset <- coefficients(DriftFit)[3]
            DriftFitSummary <- "Good"
          }
        }
        
        # Set poor fit value to -1 to avoid error in the while-loop
        if(sum(is.na(c(Drift, DriftDcoef, Drift_offset))) > 0){
          Drift <- -1
          DriftDcoef <- -1
          Drift_offset <- -1
        }
        
        # Basic iterFit, without anchor
        j <- 1
        while( (Drift <= 0) & (min(sum(!is.na(m)), DriftFit_dt_atMost) - j) >= DriftFit_dt_atLeast & 
               (min(sum(!is.na(m)), DriftFit_dt_atMost) - j) > 2 ){
          
          x <- c(1:( min(sum(!is.na(m)), DriftFit_dt_atMost) - j)) * frameRate_sec
          y <- unlist(m[1:( min(sum(!is.na(m)), DriftFit_dt_atMost) - j)])
          
          DriftFit. <- try(
            DriftFit <- nls(y ~ 4 * DriftDcoef * x + (Drift * x)^2 + Drift_offset, trace = F,
                            start=c(Drift = 5, DriftDcoef = 1.5, Drift_offset = 0),
                            lower=c(0,0,0),
                            upper=c(50,50,0.5), nls.control(maxiter = 2000), algorithm = "port"), 
            silent = TRUE)
          
          if(class(DriftFit.) == "try-error"){
            Drift <- -1
            DriftDcoef <- -1
            Drift_offset <- -1
            
          }else{
            Drift <- coefficients(DriftFit)[1]
            DriftDcoef <- coefficients(DriftFit)[2]
            Drift_offset <- coefficients(DriftFit)[3]
            DriftFitSummary <- paste("Decreased", j, "dt")
          }
          
          if(sum(is.na(c(Drift, DriftDcoef, Drift_offset))) > 0){
            Drift <- -1
            DriftDcoef <- -1
            Drift_offset <- -1
          }
          j <- j + 1
        }
        
        # Set poor fit value to -1 to avoid error
        if(sum(is.na(c(Drift, DriftDcoef, Drift_offset))) > 0){
          Drift <- -1
          DriftDcoef <- -1
          Drift_offset <- -1
        }
        
        
        # If poor fit    
        if(Drift <= 0 |DriftDcoef < 0){
          Drift <- Inf
          DriftDcoef <- Inf
          Drift_offset <- Inf
          DriftFitSummary <- "Fail to obatin good fit"
        }
        
        # If it is a simple flow model without DriftDcoef 
        if(DriftDcoef == 0){
          DriftFitSummary <- "Drift with no diffusion"
        }
        
        #Summarize
        data.frame(Drift, DriftDcoef, Drift_offset, DriftFitSummary)
      })    
      
      return(Drift)
    }
    Drift <- getDrift(MSD_Alpha, frameRate_sec = frameRate_sec, DriftFit_dt_atMost = DriftFit_dt_atMost, DriftFit_dt_atLeast = DriftFit_dt_atLeast)
    
    # Calculate the angle between two displacements for each dt
    cat("  Calculating Angle...\n")
    getAngle <- function(trackll, dt = 30){    
      
      
      Angle <-pbapply::pblapply(trackll, function(x) {
        
        Angle <- c()
        for (i in 1:dt){
          v1 <- x[c("x", "y")] - dplyr::lag(x[c("x", "y")], i)
          v2 <- dplyr::lag(v1, i)
          
          theta <- (atan2(v1[,2], v1[,1]) - atan2(v2[,2], v2[,1]))
          # From -180:180 to 0:360
          theta[theta < 0 & !is.na(theta)] <- 2*pi + theta[theta < 0 & !is.na(theta)]
          theta <- theta*180/pi
          
          # Note the relationship between d1, d2, v1, v2
          d1 <- sqrt(v2$x^2 + v2$y^2)*pixel_size_um
          d2 <- sqrt(v1$x^2 + v1$y^2)*pixel_size_um
          d2[which(is.na(d1))] <- NA # Just to make sure that there are values for both displacements
          dAvg <- (d1 + d2) / 2
          
          Angle <- cbind(Angle, theta, d1, d2, dAvg)
          
        }
        
        colnames(Angle) <- matrix(sapply(seq(1:dt), function(x) c(paste0("angle_dt", x),
                                                                  paste0("angleDisp1_dt", x),
                                                                  paste0("angleDisp2_dt", x),
                                                                  paste0("angleDispAvg_dt", x))), ncol = 1)
        
        return(Angle)
        
        
      })
      
      trackll_angle <- mapply(cbind, trackll, Angle, SIMPLIFY = FALSE)
      
      
      return(trackll_angle)
    }
    trackll_basic_angle <- getAngle(trackll_basic, dt = dt)
    
    # Calculate AngleBias

    cat("  Calculating Angle bias...\n")
    getAngleBias <- function(trackll_basic_angle, AngleBiasWeight = 25, AngleBias_dt_atLeast = 1){
      
      AngleBias <- pbapply::pblapply(trackll_basic_angle, function(x){
        
        angle_dt <- x[, !is.na(stringr::str_match(names(x), "angle_dt"))]
        
        # To remove column with All NAs
        all_na <- function(x) any(!is.na(x))
        angle_dt_complete <- angle_dt %>% dplyr::select_if(all_na)
        
        if(ncol(angle_dt_complete) < AngleBias_dt_atLeast){
          AngleBias <- NA
          AngleBiasWeight <- NA
          cosAngle <- NA
          sinAngle <- NA
          f180.0 <- NA
        }else{
          theta <- angle_dt_complete[!is.na(angle_dt_complete)]
          # Change angles from 0:360 to 0:180 [angleRemap() changes angles from 0:360 to -180:180]
          theta <- abs(oce::angleRemap(theta))
          AngleBias <-  mean( (1 - 2/(1 + exp((-theta + 90) / AngleBiasWeight)))/(1 - 2/(1 + exp(90 / AngleBiasWeight))) )
          cosAngle <- mean(cos(theta * pi / 180))
          sinAngle <- mean(sin(theta * pi / 180))
          f180.0 <- sum(theta >= 150) / sum(theta <= 30)
        }    
        
        # Just get the angles from dt1
        if(ncol(angle_dt_complete) == 0){
          AngleBias_dt1 <- NA
          cosAngle_dt1 <- NA
          sinAngle_dt1 <- NA
          f180.0_dt1 <- NA
        }else{
          theta_dt1 <- angle_dt_complete$angle_dt1[!is.na(angle_dt_complete$angle_dt1)]
          theta_dt1 <- abs(oce::angleRemap(theta_dt1))
          AngleBias_dt1 <-  mean( (1 - 2/(1 + exp((-theta_dt1 + 90) / AngleBiasWeight)))/(1 - 2/(1 + exp(90 / AngleBiasWeight))) )
          cosAngle_dt1 <- mean(cos(theta_dt1 * pi / 180))
          sinAngle_dt1 <- mean(sin(theta_dt1 * pi / 180))
          f180.0_dt1 <- sum(theta_dt1 >= 150) / sum(theta_dt1 <= 30)
        }
        
        
        data.frame(AngleBias, AngleBias_dt1, AngleBiasWeight, cosAngle, cosAngle_dt1, sinAngle, sinAngle_dt1, f180.0, f180.0_dt1)
      })
      return(AngleBias)
    }
    AngleBias <- getAngleBias(trackll_basic_angle, AngleBiasWeight = AngleBiasWeight, AngleBias_dt_atLeast = AngleBias_dt_atLeast)
    
    
    
    
    
    
    
    # Combine all the parameters
    cat("  Combining results...\n")
    
    # Get the trackID, states and transit
    # For states, if the data is not processed by trackll.HMM.splitTrack(), it only shows the states for the first displacement
    trackID <- lapply(trackll_basic, function(x) x$trackID[1])
    states <- lapply(trackll_basic, function(x) x$states[1])
    transit <- lapply(trackll_basic, function(x) x$transit[1])
    mask_radius_um <- lapply(trackll_basic, function(x) x$mask_radius_um[1])
    
    # Dcoef, Alpha, Rc, Drift, AngleBias, TrackLength, maxDisplacement, avgDisplacement (Basic unit: per Track)
    Temporal <- mapply(function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) cbind(trackID = a, states = b, transit = c, mask_radius_um = d,e,f,g,h,i,j,k,l,m,n,o), 
                       trackID, states, transit, mask_radius_um, Dcoef, Alpha, Rc, Drift, AngleBias, TrackLength, DispFeature, Gaussianity, Boundedness, Trappedness, MSDratio,SIMPLIFY = FALSE)
    
    # trackll, Displacement, Angle, AngleDisplacement (Basic unit: per Frame)
    Displacement_Angle <- mapply(function(a,b) inner_join(a,b, by = c("Frame", "x", "y", "z", "trackID", "states", "transit", "mask_radius_um")), trackll_basic_disp, trackll_basic_angle, SIMPLIFY = FALSE)
    
    
    trackll_analyzed[[trackll_filename[iterFile]]] <- list(Temporal = Temporal, 
                                                           Spatial = trackll[[iterFile]], 
                                                           Displacement_Angle = Displacement_Angle, 
                                                           MSD = MSD)
    
  }# End of master for loop

  return(trackll_analyzed)
}


trackll_analyzed_HMM <- track.analyze(trackll_spatial_HMM_splitTrack)

######################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################


############################################################## Messy script for parameter plotting and spatiotemporal mapping ###################################################################################################################


########################################## Step 1: Data Cleaning ######################################################################################################
# Need to define iterFile as 1. This is to prepare batch processing of mutltiple trackll_analyzed
iterFile <- 1

# Extract the Temporal data from trackll_analyzed----------------------------------------------------------------------------------------------------------------
Temporal_data <- dplyr::bind_rows(trackll_analyzed_HMM[[iterFile]]$Temporal, .id = "Trajectory")
Temporal_data$Trajectory <- as.numeric(Temporal_data$Trajectory)
Temporal_data$states <- as.numeric(Temporal_data$states)
Temporal_data$Name <- names(trackll_analyzed_HMM)[[iterFile]]
rownames(Temporal_data) <- 1:nrow(Temporal_data)

# Cleaning and imputing values
Temporal_data_clean <- Temporal_data[!is.infinite(Temporal_data$Alpha) & 
                                       !is.infinite(Temporal_data$Dcoef) &
                                       !is.na(Temporal_data$Rc) , ]

Temporal_data_clean[is.infinite(Temporal_data_clean$Rc), "Rc"] <- Temporal_data_clean[is.infinite(Temporal_data_clean$Rc), "mask_radius_um"]
Temporal_data_clean[is.infinite(Temporal_data_clean$Rc_Nuc_ratio), "Rc_Nuc_ratio"] <- 1
Temporal_data_clean[is.infinite(Temporal_data_clean$Drift), "Drift"] <- 0


# Average radius of the nucleus
mean_mask_radius_um <- mean(unique(bind_rows(trackll_analyzed_HMM[[iterFile]]$Temporal)$mask_radius_um))
########################################################################################################################################################################


########################################## Step 2: Dimension reduction and Clustering ######################################################################################################
set.seed(123)

# Check plot for hyperparameter tuning
umap <- uwot::umap(scale(Temporal_data_clean[, c(7, 11, 24, 35, 36, 37, 39, 43, 44, 47)]), n_neighbors = 13, min_dist = 0.001, spread = 1, y = Temporal_data_clean$states, target_weight = 0.7)
umap.df <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2])

# Plot data to decide number of cluster 
plot(umap.df)

# GMM
set.seed(123)
mc <- mclust::Mclust(umap.df, G = 2, modelNames = "VVV")
plot(mc, what = "classification")
plot(mc, what = "BIC")
cluster.class <- mc$classification

# Combine umap dataframe with kmean clustering data
umap.df_classified <- bind_cols(umap.df, data.frame(Class = cluster.class))
umap.df_classified$Class <- factor(umap.df_classified$Class)

# Combine Temporal data with clustering data
Temporal_data_clean_classified <- bind_cols(Temporal_data_clean, data.frame(Class = cluster.class))
Temporal_data_clean_classified$Class <- factor(Temporal_data_clean_classified$Class)

# Calculate meanParameter for classes, and use that to reorder the Class--------------------------------
meanParameter <- Temporal_data_clean_classified %>% 
  dplyr::group_by(Class) %>% 
  dplyr::mutate(meanParameter = mean(Alpha)) %>%     # Suggest using logDcoef or Alpha
  dplyr::ungroup() %>%
  dplyr::select (-Class)

class.order <- Temporal_data_clean_classified %>% 
  dplyr::group_by(Class) %>% 
  dplyr::summarize(meanParameter = mean(Alpha)) %>%    # Suggest using logDcoef or Alpha
  dplyr::arrange(meanParameter) %>% 
  dplyr::mutate(Class = 1:n())

Temporal_data_clean_classified <- dplyr::left_join(meanParameter, class.order, by = "meanParameter") %>%
  dplyr::select (-meanParameter)
Temporal_data_clean_classified$Class <- factor(Temporal_data_clean_classified$Class)

# Reorder class number and add UMAP1, UMAP2 to Temporal_data_clean_classified
Temporal_data_clean_classified <- cbind(Temporal_data_clean_classified, umap.df)

# The class number will be reorder later
ggplot2::ggplot(data = Temporal_data_clean_classified, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 2) + 
  theme_bw() +
  ggsci::scale_color_startrek() +
  coord_fixed(ratio = 1) + 
  labs(title = "Preliminary (Class number is not final)", x = "UAMP1", y = "UMAP2")

ggplot2::ggplot(data = Temporal_data_clean_classified, aes(x = UMAP1, y = UMAP2, color = logDcoef)) +           #Add "-" to UMAP1 or UMAP2 to flip the axis
  geom_point(size = 2) + 
  geom_density_2d(aes(group = Class), color = "white", alpha = 0.5) +
  #viridis::scale_color_viridis(option = "D") + 
  scale_color_gradientn(colours = pals:: parula(200), guide = "colourbar", breaks=c(seq(-2.5,0.5,1)),limits=c(-2.5,1)) +
  labs(x = "UAMP1", y = "UMAP2") +
  ggdark::dark_theme_gray() +
  coord_fixed(ratio = 1)  + xlim(-10,10) + ylim(-10,10) +
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggdark::invert_geom_defaults()    # Invert ggplot default back to normal after using dark_theme_gray

Temporal_data_clean_classified$Class <- as.numeric(as.character(Temporal_data_clean_classified$Class))
Temporal_data_clean_classified$Class <- as.factor(Temporal_data_clean_classified$Class)

ggplot2::ggplot(data = Temporal_data_clean_classified, aes(x = UMAP1, y = UMAP2, color = Class)) +  #Add "-" to UMAP1 or UMAP2 to flip the axis
  geom_point(size = 2, alpha = 0.2) + 
  geom_density_2d() + 
  theme_bw() +
  scale_colour_manual(values = c("#5C88DAFF", "#CC0C00FF")) +
  # ggsci::scale_color_startrek() +
  coord_fixed(ratio = 1) + 
  labs(x = "UAMP1", y = "UMAP2") + xlim(-10,10) + ylim(-10,10) +
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

# ----------------------------------------------------------------------------------------------------------

# Filter useful parameters for plotting (Must include "Trajectory" and "Class" columns)
# names(Temporal_data_clean_classified)
Temporal_data_clean_classified_filtered <- Temporal_data_clean_classified[, c(1,7,11,18,20,24,40,53)]

# Pivot the data into long form for facet
Temporal_data_clean_classified_filtered_long <- Temporal_data_clean_classified_filtered %>%
  tidyr::pivot_longer(!c(Trajectory, Class), names_to = "parameter", values_to = "value")

# Add median value to boxplot
meds <- plyr::ddply(Temporal_data_clean_classified_filtered_long, .(parameter, Class), summarize, med = median(value))

# Individual parameter plot----------------------------------------------------------------------------------------------------------------
ggplot2::ggplot(Temporal_data_clean_classified_filtered_long, aes(x = value, y = Class, fill = Class)) +
  geom_boxplot(color = "black", aes(fill = Class), outlier.shape = 20, outlier.alpha = 0.5) + coord_flip() +
  # geom_jitter(color = "black", size = 0.2, alpha = 0.1) + # Show dot on boxplot
  facet_wrap(~parameter,  scales = "free") +
  geom_text(data = meds, aes(x = med, , y = Class, label = round(med, 2)), size = 3, vjust = -0.5, inherit.aes = FALSE, color = "black") +
  ggsci::scale_fill_startrek() +
  theme_bw()

static.locError_um <- sqrt(mean(Temporal_data_clean_classified[Temporal_data_clean_classified$Class == 1, ]$DcoefFit_intercept) / 4)

##################################################### Step 3: Ensemble MSD temporal analysis #########################################################################################
# Individual and ensemble MSD plot----------------------------------------------------------------------------------------------------------------
dt <- 30
frameRate_sec <- 10/1000
# Extract the MSD data from trackll_analyzed
MSD_data <- lapply(trackll_analyzed_HMM[[iterFile]]$MSD, function(x) t(x))


# Organizating the data structure
MSD_data <- dplyr::bind_cols(MSD_data)
MSD_data <- t(MSD_data)
colnames(MSD_data) <- 1:dt
MSD_data <- data.frame(cbind(Trajectory = 1:nrow(MSD_data), MSD_data))

# Combine with class data
MSD_data_classified <- dplyr::left_join(MSD_data, Temporal_data_clean_classified[, c("Trajectory", "Class")], by = "Trajectory")

# Pivot longer
MSD_data_classified_long <- MSD_data_classified %>% 
  tidyr::pivot_longer(!c(Trajectory, Class), names_to = "dt", values_to = "MSD")

# Change dt number to sec
MSD_data_classified_long$dt <- readr::parse_number(MSD_data_classified_long$dt) * frameRate_sec

# Remove data with no class for a pretty plot
MSD_data_classified_long <- MSD_data_classified_long[-which(is.na(MSD_data_classified_long$Class)), ]

#MSD_data_classified_long <- MSD_data_classified_long %>% filter (Trajectory %in% static$Trajectory)


# Ensemble MSD temporal analysis
ensembleMSD_data <- MSD_data_classified_long %>%
  group_by(Class, dt) %>%
  dplyr::summarize(ensembleMSD = mean(MSD, na.rm = TRUE),
                   SEM = sd(MSD, na.rm = TRUE) / sqrt(sum(!is.na(MSD))),
                   .groups = 'drop')

ensembleMSD_data_clean <- ensembleMSD_data[complete.cases(ensembleMSD_data), ] # Remove row with no SD (just one data point)

# Converting to the format as generated by getMSD()
ensembleMSD_data_wide <- ensembleMSD_data %>% 
  dplyr::select(-SEM) %>%
  tidyr::pivot_wider(names_from = dt, values_from = ensembleMSD)

colnames(ensembleMSD_data_wide) <- c("Class", sapply(seq(1:dt), function(x) paste0("MSD_dt", x))  )
ensembleMSD_data_wide_list <- split( ensembleMSD_data_wide[-which(names(ensembleMSD_data_wide)=="Class"|names(ensembleMSD_data_wide)=="SD")], ensembleMSD_data_wide$Class )
ensembleMSD_data_wide_list <- lapply(ensembleMSD_data_wide_list, function(x) t(unlist(x)))

# Calculate ensemble behavior (Need to the getXXX() functions to environment first)----------------------------------------------------------------------------------------------------------------
Dcoef_ensemble <- get2D_Dcoef(ensembleMSD_data_wide_list, locError_um = static.locError_um)
Alpha_ensemble <- getAlpha(ensembleMSD_data_wide_list)
MSD_Alpha_ensemble <- mapply(cbind, ensembleMSD_data_wide_list, Alpha_ensemble, SIMPLIFY = FALSE)
Rc_ensemble <- getRc(MSD_Alpha_ensemble, Rc_max = mean_mask_radius_um, locError_um = static.locError_um) # mean_mask_radius_um is calculated at the beginning of this script
Drift_ensemble <- getDrift(MSD_Alpha_ensemble)

Temporal_ensembleMSD <- mapply(function(a,b,c,d) cbind(a,b,c,d), 
                               Dcoef_ensemble,Alpha_ensemble,Rc_ensemble,Drift_ensemble, SIMPLIFY = FALSE)


# This is the temporal data from ensemble MSD. Still need to append with f180.0 anisotropy data on Step 3.
Temporal_ensembleMSD <- bind_rows(Temporal_ensembleMSD, .id = "Class")
rownames(Temporal_ensembleMSD) <- NULL

# Calculate mean of other parameters and count to be put into Temporal_ensemble
mean_ensemble <- Temporal_data_clean_classified %>%
  dplyr::group_by(Class) %>%
  dplyr::summarize(meanAngleBias = mean(AngleBias),
                   meanStraightness = mean(Straightness),
                   count = n())

Temporal_ensemble <- inner_join(Temporal_ensembleMSD, mean_ensemble,by = "Class")

# Plot MSD and ensemble MSD----------------------------------------------------------------------------------------------------------------
ggplot2::ggplot() +
  #geom_line(data = MSD_data_classified_long, aes(x = dt, y = MSD, color = Class, group = Trajectory), alpha = 0.1) +   # This line plot all the individual MSD
  #geom_smooth(data = MSD_data_classified_long, aes(x = dt, y = MSD, group = Class, color = Class), alpha = 0.1) + # This line plot a smooth curve for raw MSD
  geom_point(data = ensembleMSD_data_clean, aes(x = dt, y = ensembleMSD)) +
  geom_line(data = ensembleMSD_data_clean, aes(x = dt, y = ensembleMSD, color = Class)) +
  geom_ribbon(data = ensembleMSD_data_clean, aes(x = dt, y = ensembleMSD, ymin = ensembleMSD-SEM, ymax = ensembleMSD+SEM, fill = Class), alpha = 0.1, colour = NA) +
  # scale_x_log10() +
  # scale_y_log10() +
  # ylim(0, 1) +
  xlim(0, 0.2) +
  annotation_logticks() +  # log scale looks nicer
  ggsci::scale_color_startrek() + 
  ggsci::scale_fill_startrek() + 
  labs( x = "Time interval (s)", y = expression(paste("MSD"," (", mu, m^2,")"))) +
  theme_bw() 

######################################################### Step 4: Angle ##########################################################################################
# Angle Rose plot

# Extract the Displacement_Angle data from trackll_analyzed----------------------------------------------------------------------------------------------------------------
Displacement_Angle_data <- dplyr::bind_rows(trackll_analyzed_HMM[[iterFile]]$Displacement_Angle, .id = "Trajectory")
Displacement_Angle_data$Trajectory <- as.numeric(Displacement_Angle_data$Trajectory)
Displacement_Angle_data_classified <- dplyr::full_join(Displacement_Angle_data, 
                                                       Temporal_data_clean_classified[, c("Trajectory", "Class")], 
                                                       by = "Trajectory")

# Find column with name containing "angle"
angle_column <- which(!is.na(stringr::str_match(names(Displacement_Angle_data_classified), "angle")))    

# Extract the Angle data from Displacement_Angle_data, "1" is the Trajectory column
Angle_data_classified <- Displacement_Angle_data_classified[, c(1, angle_column, which(colnames(Displacement_Angle_data_classified) == "Class"))]

# Pivot longer, based on "_" as the separator for angle_dt, angleDisp1_dt, angleDisp2_dt and angleDispAvg_dt
Angle_data_classified_long <- Angle_data_classified %>% 
  tidyr::pivot_longer(!c(Trajectory, Class), names_sep = "_", names_to = c("parameter", "dt"), values_to = "value")

Angle_data_classified_long$dt <- readr::parse_number(Angle_data_classified_long$dt)

# Create a unique identifier row for each name and then use pivot_wider
# https://stackoverflow.com/questions/58837773/pivot-wider-issue-values-in-values-from-are-not-uniquely-identified-output-w
Angle_data_classified_long_wide <- Angle_data_classified_long %>%
  dplyr::group_by(Trajectory) %>%
  dplyr::mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
  dplyr::select(-row)

# Need to move the columns for tidy data
Angle_data_classified_long_wide$angle <- dplyr::lag(Angle_data_classified_long_wide$angle, 3)
Angle_data_classified_long_wide$angleDisp1 <- dplyr::lag(Angle_data_classified_long_wide$angleDisp1, 2)
Angle_data_classified_long_wide$angleDisp2 <- dplyr::lag(Angle_data_classified_long_wide$angleDisp2, 1)

Angle_data_classified_movement <- Angle_data_classified_long_wide %>%  
  dplyr::mutate(Movement = dplyr::case_when(
    angle >=150 & angle <= 210 ~ "Backward",
    angle >=330 & angle <= 360  ~ "Forward",
    angle >=0 & angle <= 30  ~ "Forward",
    angle > 30 & angle < 180 ~ "in-between",
    angle >210 & angle < 330 ~"in-between")
  )

# Remove data with no class or no movement for a pretty plot
Angle_data_plot <- Angle_data_classified_movement[-which(is.na(Angle_data_classified_movement$Class) | is.na(Angle_data_classified_movement$Movement)), ]


# Merge f180 / 0 angle data with Temporal_ensemble generated on Step 2
f180.0 <- Angle_data_plot %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(f180.0 = sum(Movement == "Backward") / sum(Movement == "Forward"))
Temporal_ensemble_FINAL <- cbind(Temporal_ensemble, f180.0)

# facet warp does not work for coord_polar, need to use lappy and grid.arrange
# startrek palette only have 7 colors, if there are more than 7 classes, the script will break
angle_plots <- lapply(seq_along(unique(Angle_data_plot$Class)), function(x){ 
  
  ggplot2::ggplot(Angle_data_plot[as.numeric(Angle_data_plot$Class) == x, ], 
                  aes(x = angle, alpha = Movement)) +   #, fill = Movement
    geom_histogram(binwidth = 5, fill = ggsci::pal_startrek(palette = c("uniform"), alpha = 1)(7)[x], color = "black") +
    scale_alpha_manual(values = c(1, 0.3, 0)) +
    coord_polar(start = 90*pi/180) +
    scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, 30)) +
    labs(title = paste0("Class ", x))+
    theme(panel.background = element_rect(fill = NA, color = "black"),
          panel.grid.major = element_line(colour = "grey90"),
          panel.grid.minor = element_line(colour = "grey90"),
          axis.text.x = element_text(size = 6))
})

# Plot angle plots for all classes----------------------------------------------------------------------------------------------------------------
do.call(gridExtra::grid.arrange,  angle_plots)


#########################################################  Step 5: Spatiotemporal mapping ###################################################################################

track_spatial <- track.merge(trackll_spatial)
track_spatial <- dplyr::bind_rows(track_spatial, .id = "name")

# Making the error bar for nucleolus, y is 0 because it is already rotated to the x-axis
landmark_error_bar <- data.frame(x = mean(unique(track_spatial$landmark_centroid_x.normalized.centered.rotated)),
                                 y = 0,
                                 sd = sd(unique(track_spatial$landmark_centroid_x.normalized.centered.rotated)))


# Calculate the outline of a "standard nucleolus"
nucleolus_centroid <- mean(unique(track_spatial$landmark_centroid_x.normalized.centered.rotated))
i <- 0
centroid <- data.frame(x = -0.5119, y = 0)       # -0.5119, 0 is the coordinate of the centroid of a standard nucleolus
if(round(nucleolus_centroid, 4) < -0.5119){
  
  while( round(nucleolus_centroid, 4) < round(mean(centroid$x), 4) ) {
    
    i <- i + 0.0001
    
    # 1.34 and 1.67 are from Biophysical journal, 2017, Supporting Figure 6)
    d <- 1.34 + i
    n <- 1.67 + i
    
    bigCircle <- circleFun(center = c(d,0),diameter = n*2, npoints = 5000)
    smallCircle <- circleFun(center = c(0,0),diameter = 1*2, npoints = 5000)
    
    # https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    a <-  (n^2 - 1^2 + d^2) / (2*d)
    h <- sqrt(n^2 - a^2)
    P0 <- c(d,0)
    P1 <- c(0,0)
    P2 <-  P0 + a*(P1 -P0) / d
    x3 <- P2[1] + h * (0 - 0) /d
    
    b1 <- which.min((abs(bigCircle$x - x3)))
    b2 <- which.min((abs(bigCircle$x[-b1] - x3)))
    s1 <- which.min((abs(smallCircle$x - x3)))
    s2 <- which.min((abs(smallCircle$x[-s1] - x3)))
    centroid <- rbind(bigCircle[c(b2:b1), ], smallCircle[c(s2:s1), ])
  }
}
if(round(nucleolus_centroid, 4) > -0.5119){
  
  while( round(nucleolus_centroid, 4) > round(mean(centroid$x), 4) ) {
    
    i <- i - 0.0001
    
    # 1.34 and 1.67 are from Biophysical journal, 2017, Supporting Figure 6)
    d <- 1.34 + i
    n <- 1.67 + i
    
    bigCircle <- circleFun(center = c(d,0),diameter = n*2, npoints = 5000)
    smallCircle <- circleFun(center = c(0,0),diameter = 1*2, npoints = 5000)
    
    # https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    a <-  (n^2 - 1^2 + d^2) / (2*d)
    h <- sqrt(n^2 - a^2)
    P0 <- c(d,0)
    P1 <- c(0,0)
    P2 <-  P0 + a*(P1 -P0) / d
    x3 <- P2[1] + h * (0 - 0) /d
    
    b1 <- which.min((abs(bigCircle$x - x3)))
    b2 <- which.min((abs(bigCircle$x[-b1] - x3)))
    s1 <- which.min((abs(smallCircle$x - x3)))
    s2 <- which.min((abs(smallCircle$x[-s1] - x3)))
    centroid <- rbind(bigCircle[c(b2:b1), ], smallCircle[c(s2:s1), ])
  }
}


# Set plot mid point of a track or not
MID <- FALSE # (For plotting center of tracks, set to TRUE)

####################################################################################################################################################################################
# Pre-process for Class plotting analysis

# Extract the Spatial data from trackll_analyzed/trackll_analyzed_HMM----------------------------------------------------------------------------------------------------------------
Spatial_data <- dplyr::bind_rows(trackll_analyzed_HMM[[iterFile]]$Spatial, .id = "Trajectory")
Spatial_data$Trajectory <- as.numeric(Spatial_data$Trajectory)
Spatial_data$states <- as.numeric(Spatial_data$states)
Spatial_data$Name <- names(trackll_analyzed_HMM)[[iterFile]]

# Extract single-step Dcoef from trackll_analyzed/trackll_analyzed_HMM----------------------------------------------------------------------------------------------------------------
ssDcoef_data <- dplyr::bind_rows(trackll_analyzed_HMM[[iterFile]]$Displacement_Angle, .id = "Trajectory")
ssDcoef_data <- ssDcoef_data[ ,c("Trajectory", "Frame", "x", "y", "z", "SSDcoef", "logSSDcoef")]
ssDcoef_data$Trajectory <- as.numeric(Spatial_data$Trajectory)
ssDcoef_data$Name <- names(trackll_analyzed_HMM)[[iterFile]]

############### THIS IS THE ULTIMATE SPATIOTEMPORAL DATA ############### 
# Spatial_data contains tracks position after processed by HMM and temporal analysis (Before classified, and represent all "raw" detection)
# Temporal_data_clean_classified contain class data
SpatioTemporal_data <- inner_join(Spatial_data, ssDcoef_data, by = c("Trajectory", "Frame", "x", "y", "z", "Name"))
SpatioTemporal_data <- inner_join(SpatioTemporal_data, subset(Temporal_data_clean_classified, select = -c(states, transit, mask_radius_um)), by = c("Trajectory", "trackID", "Name"))
SpatioTemporal_data$cellID <- paste(SpatioTemporal_data$fileNum, SpatioTemporal_data$mask, sep = "_")    # Add Cell ID (fileNum_mask)
############### THIS IS THE ULTIMATE SPATIOTEMPORAL DATA ############### 

# Making the error bar for nucleolus, y is 0 because it is already rotated to the x-axis
landmark_error_bar <- data.frame(x = mean(unique(SpatioTemporal_data$landmark_centroid_x.normalized.centered.rotated)),
                                 y = 0,
                                 sd = sd(unique(SpatioTemporal_data$landmark_centroid_x.normalized.centered.rotated)))

# Flip and symmetric 
SpatioTemporal_data_flip <- SpatioTemporal_data %>% dplyr::mutate(y.normalized.centered.rotated = -y.normalized.centered.rotated)
SpatioTemporal_data_symmetric <- rbind(SpatioTemporal_data, SpatioTemporal_data_flip)

SpatioTemporal_data.mid <- SpatioTemporal_data %>% distinct(Trajectory, .keep_all = TRUE)
SpatioTemporal_data_flip.mid <- SpatioTemporal_data.mid %>% dplyr::mutate(y.normalized.centered.rotated.mid = -y.normalized.centered.rotated.mid)
SpatioTemporal_data_symmetric.mid <- rbind(SpatioTemporal_data.mid, SpatioTemporal_data_flip.mid)

####################################################################################################################################################################################
# Number of near neighbor
# r2: radius for nearest neighbor(um, if assuming 1 AU in the map is 1 um; check "mean_mask_radius_um")
r2 <- 0.15^2

SpatioTemporal_data_symmetric.class <- SpatioTemporal_data_symmetric %>% dplyr::group_by(Class) %>%
  dplyr::mutate(n_neighbors = count_neighbors(x.normalized.centered.rotated, y.normalized.centered.rotated, r2 = r2, xy = 1))

SpatioTemporal_data_symmetric.mid.class <- SpatioTemporal_data_symmetric.mid %>% dplyr::group_by(Class) %>%
  dplyr::mutate(n_neighbors.mid = count_neighbors(x.normalized.centered.rotated.mid, y.normalized.centered.rotated.mid, r2 = r2, xy = 1))

SpatioTemporal_data_symmetric.all <- SpatioTemporal_data_symmetric %>%
  dplyr::mutate(n_neighbors = count_neighbors(x.normalized.centered.rotated, y.normalized.centered.rotated, r2 = r2, xy = 1))

SpatioTemporal_data_symmetric.mid.all <- SpatioTemporal_data_symmetric.mid %>%
  dplyr::mutate(n_neighbors.mid = count_neighbors(x.normalized.centered.rotated.mid, y.normalized.centered.rotated.mid, r2 = r2, xy = 1))

MID <- FALSE

####################################################################################################################################################################################
# Spatiotemporal mapping (grid map)
node <- as.data.frame(tidyr::crossing(x = seq(-1,1,0.01), y= seq(-1,1,0.01)))

# Mask out node that is outside the nucleus
degree <- seq(0,360,0.001)
rad <- degree * (pi/180)
r <- seq(1.01, sqrt(2), 0.01)    # Set to 1.01 instead of 1 to preserve the nodes "on" the nuclear outline
x = round(r  * cos(rad), 2)
y = round(r * sin(rad), 2)

circle.mask <- matrix(cbind(x, y), ncol = 2)

circle.mask.matching <- paste(circle.mask[ ,1], circle.mask[ ,2])
node.matching <- paste(round(node[ ,1], 2), round(node[ ,2], 2))

node.masked <- node[!node.matching %in% circle.mask.matching, ]
node <- node.masked   # Note: Because of this line, node is the same as node.masked

#----------------------------------------------------------------------------------------------------------------------------------------------------------


# Grid map plots using lapply (It may use up all the memory)
grid.plots <- lapply(c(seq_along(unique(SpatioTemporal_data_symmetric$Class)), 2), function(x){ 
  
  # Data, Class 0: All data
  if(x == 0){
    df <- SpatioTemporal_data_symmetric
    
  }else{  
    df <- SpatioTemporal_data_symmetric[SpatioTemporal_data_symmetric$Class == x, ]
  }
  
  point <- df[, c("x.normalized.centered.rotated", "y.normalized.centered.rotated")]
  
  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  
  # Choose to plot parameter ("param") or count ("count")
  colored_by <- "count"  
  
  # Search radius = 0.3 for parameter, 0.15 for count(density)
  dist.matrix <- fields::rdist(node.masked, point)
  dist.matrix[dist.matrix > 0.15] <- NA
  dist.matrix[dist.matrix < 0.15] <- 0
  
  # Modify the number below to change parameters to plot (colnames(SpatioTemporal_data_symmetric))
  colored_by <- "param"
  param.index <- 49
  param.legend_label <- paste("Average", colnames(SpatioTemporal_data_symmetric)[param.index])
  ########  ########  ########  ########  ########  ########  ########  ########  ########  ########  
  
  
  # Put in parameter value if the point fit into the circle created by the node
  dist.matrix.param <- dist.matrix + matrix(rep(df[ ,param.index], nrow(node.masked)), nrow = nrow(node.masked), byrow = TRUE)
  
  # Calculate the mean of parameter
  param <- apply(dist.matrix.param, 1, mean, na.rm=TRUE)
  count <- apply(!is.na(dist.matrix.param), 1, sum)
  
  
  gridMap.data <- cbind(node.masked, param, count)
  
  gridMap.data$count <- (gridMap.data$count*nrow(SpatioTemporal_data_symmetric))/sum(SpatioTemporal_data_symmetric$Class != 1) 
  # Multiply nrow(SpatioTemporal_data_symmetric) is to counter-act how count is being normalized in the script at the beginning
  
  
  # Set limit !
  
  p <-  ggplot2::ggplot() +
    geom_point(data = gridMap.data, aes_string("x","y", color = colored_by), size = 2) + 
    #max(gridMap.data$count, na.rm = TRUE)
    {if(colored_by == "count") scale_color_distiller(palette = "Spectral", direction = -1, na.value = "transparent")} +
    {if(colored_by == "param") scale_color_gradientn(colors = as.character(inlmisc::GetColors(n = 512, scheme = "sunset", bias = 0.8)))} +
    
    geom_point(data = landmark_error_bar, aes(x, y), color = 'black', size = 1.5) +
    geom_errorbarh(data = landmark_error_bar, aes(xmax = x + sd, xmin = x - sd, y = 0), height = 0.01, color = 'black') + 
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 1), lty = 3, fill = "white", alpha = 0) +
    geom_path(aes(x,y), size = 1, alpha = 1, data = bigCircle[c(b2:b1), ], col = "black", lty = 3) +      # Internal outline of nucleolus
    
    {if(colored_by == "count")   
      labs(title = paste("Class", x), 
           subtitle = paste(c(0, max(gridMap.data$count, na.rm = TRUE)), collapse = " - "),
           x = "x (\u03BCm)", y = "y (\u03BCm)", 
           color = "Density")} + 
    {if(colored_by == "param")   
      labs(title = paste("Class", x), 
           subtitle = paste(c(min(gridMap.data$param, na.rm = TRUE), max(gridMap.data$param, na.rm = TRUE)), collapse = " - "),
           x = "x (\u03BCm)", y = "y (\u03BCm)", 
           color = param.legend_label)} + 
    
    xlim(-1,1) + ylim(-1, 1) +
    theme(panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA)) + 
    theme_classic(base_size = 20) +
    coord_fixed(ratio = 1)  
  # + guides(color=guide_colorbar(ticks.colour = NA)) +
  #   theme(legend.key.height=unit(9,"line")
  ggsave("plot.png", width = 15, height = 15, dpi = 300, units = "in", device='png',
         plot = egg::set_panel_size(p=p, width=unit(11, "in"), height=unit(11, "in")))
  
  
})

do.call(gridExtra::grid.arrange,  grid.plots)
