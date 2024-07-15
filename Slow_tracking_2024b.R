# Author: Yick Hin Ling, Carl Wu lab, Johns Hopkins University
# Email: yhinling@gmail.com

##################################  Step 1: Read DiaTrack .mat files ########################################################################################################################################
readDiaTrack <- function(folder){
  
  mat_files <- list.files(file.path(folder), pattern = "\\.mat$")
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

#################################### Step 2: Read mask .tif files ###########################################################################################################################################
readMask <- function(folder, resolution = c(128, 128), pixel_size_um = 16/150, MinPts = 10, eps = 3, mask.plot = TRUE, mask_outline.plot = TRUE, maskName = "MASK", landmark = FALSE){
  
  # Read _MASK.tif files and save the matrix in mask_list
  mask_files <- list.files(file.path(folder), pattern = paste0(maskName, ".tif$"))
  mask_list <- list()
  for (i in seq_along(mask_files)){
    img <- EBImage::readImage(file.path(folder, mask_files[i]))
    mat <- which(img != 0, arr.ind = TRUE, useNames = FALSE)  # Find nucleus with pixel = 1
    if(is_empty(mat)){
      mat <- matrix(c(-1,-1),1,2)
    }
    
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
mask_list.POI <- readMask(folder, resolution = c(256, 256), pixel_size_um = 16/150, maskName = "POI", mask.plot = F, mask_outline.plot = F)
mask_list.H2B <- readMask(folder, resolution = c(256, 256), pixel_size_um = 16/150, maskName = "H2B", mask.plot = F, mask_outline.plot = F)

####################################### Step 3: Mask tracks based on nuclear mask ###########################################################################################################################
track.masking <- function(trackll, mask_list, resolution = c(256, 256), pixel_size_um = 16/150, mask_track.plot = TRUE, plot_index = 1:length(trackll), maskName = "MASK", preprocess = FALSE){
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
trackll_masked.POI <- track.masking(trackll_raw, mask_list = mask_list.POI, resolution = c(256, 256), pixel_size_um = 16/150, maskName = "POI", preprocess = F, mask_track.plot = F)
trackll_masked.H2B <- track.masking(trackll_raw, mask_list = mask_list.H2B, resolution = c(256, 256), pixel_size_um = 16/150, maskName = "H2B", preprocess = F, mask_track.plot = F)

######################################## Step 4: Remove gaps generated by masking ###########################################################################################################################
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
trackll_masked_gapRemoved.POI <- track.gapRemove(trackll_masked.POI)
trackll_masked_gapRemoved.H2B <- track.gapRemove(trackll_masked.H2B)

# Trim track
track.trim <- function(trackll, minFrame = 1, maxFrame = Inf){
  Frame.trim <- lapply(trackll, function(x){
    lapply(x, function(y){
      if (nrow(y) < minFrame | nrow(y) > maxFrame){
        y <- NULL
      }else{
        y <- y
      }
    })
  })
  track.trimmed <- lapply(Frame.trim, function(x){
    x[lengths(x) != 0]
  })
  return(track.trimmed)
}
trackll.POI.trimmed <- track.trim(trackll_masked_gapRemoved.POI, minFrame = 2, maxFrame = Inf)
trackll.H2B.trimmed <- track.trim(trackll_masked_gapRemoved.H2B, minFrame = 2, maxFrame = Inf)

#########################################   Step 5: Merge data  #############################################################################################################################################
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
trackll.POI <- files.merge(trackll.POI.trimmed)
trackll.H2B <- files.merge(trackll.H2B.trimmed)

#########################################   Step 6: Linking data  ###########################################################################################################################################
track.link <- function(trackll, gap = 2, threshold = 3, threshold_unit = "pixel", pixel_size_um = 16/150){
  
  #Default unit of threshold is pixel
  if(threshold_unit == "um"){
    threshold <- threshold / pixel_size_um
  }
  
  
  
  trackll <- bind_rows(trackll[[1]], .id = "Trajectory")
  trackll <- as.data.frame(trackll)
  trackll$Trajectory <- as.numeric(trackll$Trajectory)
  trackll$cellID <- paste(trackll$fileNum, trackll$mask, sep = "_")    # Add Cell ID (fileNum_mask)
  trackll.cell <- split(trackll, trackll$cellID)
  
  
  cl <- parallel::makeCluster(detectCores()-1)
  parallel::clusterExport(cl, c("trackll.cell", "gap", "threshold"), envir=environment())
  parallel::clusterEvalQ(cl, library(tidyverse))
  
  cat(paste0("\nCommon Name:  ", common_name(trackll$filename), "\n"))
  cat(paste0("Linking trajectories within ", gap, " dark frame(s), with step size threshold of ", threshold, " pixel (", threshold*pixel_size_um, " um)", "\n"))
  
  trackll.cell.linked <- pbapply::pblapply(trackll.cell, function(x){
    
    
    if(length(unique(x$Trajectory)) > 1){     # The matrix calculation will break if there is only one track
      
      repeat{  
        
        before <- length(unique(x$Trajectory))
        
        if(before == 1){  # The matrix calculation will break if there is only one track
          break
        }
        
        mat <- matrix(0, length(unique(x$Trajectory)), length(unique(x$Trajectory)))
        rownames(mat) <- unique(x$Trajectory)
        colnames(mat) <- unique(x$Trajectory)
        
        
        for(i in unique(x$Trajectory)){
          
          last.frame <- x %>% dplyr::filter(Trajectory == i) %>% dplyr::slice(n())  # last frame of current trajectory
          
          first.frame <- x %>% dplyr::group_by(Trajectory) %>% dplyr::slice(1)      # first frame of all trajectories
          
          within.gap <- first.frame %>% filter(Frame > last.frame$Frame & Frame <= last.frame$Frame + gap + 1)
          disp <- sqrt((within.gap$x - last.frame$x)^2 + (within.gap$y - last.frame$y)^2)
          within.disp <- disp <= threshold
          
          # Row: Trajectory, Column: displacement
          
          if(length(within.gap[within.disp, ]$Trajectory) == 0){
            mat[which(as.numeric(rownames(mat)) == i), within.gap[within.disp, ]$Trajectory] <- disp[within.disp]   #It will fill the row with zeros
          }else{
            mat[which(as.numeric(rownames(mat)) == i), which(as.numeric(colnames(mat)) %in% within.gap[within.disp, ]$Trajectory)] <- disp[within.disp]
          }
        }
        # Compete with other detections from the same frame for the closest point (If you want to do this step first, change mat)
        mat4 <- apply(t(mat), 1, function(x) {ifelse(any(x == min(x[x > 0])), min(x[x > 0]), 0)})
        mat5 <- apply(t(mat), 2, FUN = function(x) {x == mat4})
        mat6 <- t(apply(mat5, 2, FUN = function(x) {x*mat4}))
        
        
        # Each detection find its closest point first (before compete with other detections from the same frame) (If you want to do this step first, change mat6)
        mat1 <- apply(mat6, 1, function(x) {ifelse(any(x == min(x[x > 0])), min(x[x > 0]), 0)})
        mat2 <- apply(mat6, 2, FUN = function(x) {x == mat1})
        mat3 <- apply(mat2, 2, FUN = function(x) {x*mat1})
        
        
        
        
        for(i in unique(x$Trajectory)){
          if(sum(mat3[which(as.numeric(rownames(mat)) == i), ] > 0) == 0){        # No point within threshold
            x$Trajectory <- x$Trajectory
          }else{
            x[x$Trajectory == i,]$Trajectory <-  names(which(mat3[which(as.numeric(rownames(mat)) == i), ] > 0))    # Change Trajectory number iteratively
          }
          
        }
        
        x <- split(x, x$Trajectory)
        names(x) <- seq_along(names(x))
        x <- dplyr::bind_rows(x, .id = "Trajectory")
        x$Trajectory <- as.numeric(x$Trajectory)
        
        
        
        after <- length(unique(x$Trajectory))
        
        if(before == after){
          break
        }
      }
    } # end of if-else testing if the mask has one track only
    
    x$trackID.linked <- paste0(x$fileNum , "_", x$mask, "_", x$Trajectory)
    return(x)
  }, cl = cl)
  
  
  parallel::stopCluster(cl)
  
  # Tidy the data
  trackll.cell.linked <-bind_rows(trackll.cell.linked)
  trackll.cell.linked <- subset(trackll.cell.linked, select = -Trajectory)   
  trackll.cell.linked <- split(trackll.cell.linked, trackll.cell.linked$trackID.linked)
  trackll.cell.linked <- trackll.cell.linked[mixedsort(names(trackll.cell.linked))]
  names(trackll.cell.linked) <- seq_along(names(trackll.cell.linked))
  trackll.cell.linked <- dplyr::bind_rows(trackll.cell.linked, .id = "Trajectory")
  trackll.cell.linked$Trajectory <- as.numeric(trackll.cell.linked$Trajectory) 
  common.name <- common_name(trackll.cell.linked$filename)
  trackll.cell.linked <- split(trackll.cell.linked, trackll.cell.linked$Trajectory)
  trackll.cell.linked <- list(trackll.cell.linked)
  names(trackll.cell.linked) <- common.name
  
  
  return(trackll.cell.linked)
  
} 
trackll.POI.linked <- track.link(trackll.POI, gap = 1, threshold = 3, threshold_unit = "pixel")
trackll.H2B.linked <- track.link(trackll.H2B, gap = 1, threshold = 3, threshold_unit = "pixel")

#########################################   Step 7a: 2-expo fit  ############################################################################################################################################
getTrackLength <- function(trackll, frameRate_sec = 250/1000){
  TrackLength <- lapply(trackll, function(x){
    TrackLength <- nrow(x) - 1
    DwellTime <- TrackLength * frameRate_sec
    
    data.frame(TrackLength, DwellTime)
    
  })
  
  return(TrackLength)
}
get2ExpoDecay <- function(folder, trackll, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 10000, save.as = "POI", weight = FALSE, cutoff = 0){
  
  TrackLength <- getTrackLength(trackll[[1]], frameRate_sec = 250/1000)
  
  fit.value.summary.fin <- c()
  plot.output.fin <- c()
  fit.value.all.fin <- c()
  for(Nmin in Nmin){
    trackLength <- dplyr::bind_rows(TrackLength)$TrackLength
    
    trackLength <- trackLength[trackLength >= Nmin & trackLength <= Nmax]
    
    cat(paste0("Number of trajectory: ", length(trackLength)), "\n")
    cat(paste0("Nmin: ", Nmin), "\n")
    
    my_ecdf <- ecdf(trackLength) #It generate another function called my_ecdf
    trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength)))))
    
    y = trackLengthCDF.$y
    x = trackLengthCDF.$x * frameRate_sec
    
    
    #raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
    raw <- data.frame(x = x, y = y) %>% distinct(x, y)
    x <- as.vector(raw$x)
    y <- as.vector(raw$y)
    
    
    # data = data.frame(P, t)
    #weights = 1/y
    # weight[which(weight == Inf)] <- max(weight[which(weight != 
    #                                                    Inf)])
    
    
    
    
    # #####################
    # three.fit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
    #                                             start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
    #                                             lower=c(0,0,0,0,0),
    #                                             upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
    #                      ,silent = TRUE)
    # 
    # 
    # if(class(three.fit.try) == "try-error" |
    #    (1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) <= 0 |
    #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) < 0.001 |
    #    abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] < 0.001) |
    #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) < 0.001){
    #   cat(paste0("Three-component exponetnial decay cannot be fitted to the data!"), "\n")
    # 
    #   if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
    # 
    #     raw <- raw %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
    #     write.csv(raw, file.path(getwd(), paste0("1-CDF.raw_", save.as, ".csv")), row.names = FALSE)
    # 
    #     write.csv(plot.output.fin, file.path(getwd(), paste0("1-CDF.BS_", save.as, ".csv")), row.names = FALSE)
    #     write.csv(fit.value.summary.fin, file.path(getwd(), paste0("fit.value.summary_", save.as, ".csv")), row.names = FALSE)
    #     write.csv(fit.value.all.fin, file.path(getwd(), paste0("fit.value.all_", save.as, ".csv")), row.names = FALSE)
    # 
    #   }
    #   stop()
    # }
    # 
    # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) +
    #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) +
    #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
    # 
    # 
    # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
    # three.k.sort <- sort(three.k, decreasing = TRUE)
    # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
    # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
    # 
    # k1 <- three.k.sort[1]
    # k2 <- three.k.sort[2]
    # k3 <- three.k.sort[3]
    # f1 <- three.fraction.sort[1]
    # f2 <- three.fraction.sort[2]
    # f3 <- three.fraction.sort[3]
    #####################
    
    plot(x,y, main= paste0("Nmin: ", Nmin), xlab="Time (s)", ylab="1-CDF", log = "xy", xlim=c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100)), ylim=c(0.001, 1))
    #lines(x,three.new_y, col = "red", log = "xy", add = TRUE)
    # legend( x = "topright",
    #         legend = c("Three"),
    #         col = c("red"), lwd = 2, lty = c(0),
    #         pch = c(19) )
    
    ####################################################################################################################################
    
    y1 <- y[y >= cutoff]
    x1 <- x[1:length(y1)]
    if(weight){
      WEIGHT = 1/y1
      WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
    }
    
    
    if(weight != FALSE){
      two.fit.try <- try(two_compFit <- minpack.lm::nlsLM(y1 ~ (F1 * exp(-k1 * (x1-(Nmin)*frameRate_sec)) + (1 - F1) * exp(-k2 * (x1-(Nmin)*frameRate_sec))), trace = F,
                                                          start=c(F1 = 0.5, k1 = 1, k2 = 0.1),
                                                          lower=c(0,0,0), 
                                                          upper=c(1, Inf, Inf), control = nls.lm.control(maxiter = 1024), weights = WEIGHT)
                         ,silent = TRUE)
    }else{
      two.fit.try <- try(two_compFit <- minpack.lm::nlsLM(y1 ~ (F1 * exp(-k1 * (x1-(Nmin)*frameRate_sec)) + (1 - F1) * exp(-k2 * (x1-(Nmin)*frameRate_sec))), trace = F,
                                                          start=c(F1 = 0.5, k1 = 1, k2 = 0.1),
                                                          lower=c(0,0,0), 
                                                          upper=c(1, Inf, Inf), control = nls.lm.control(maxiter = 1024))
                         ,silent = TRUE)            
    }
    
    
    
    
    
    
    # one.fit.try <- try(one_compFit <- minpack.lm::nlsLM(y ~ exp(-k1 * (x-Nmin*frameRate_sec)), trace = F,
    #                                                     start=c(k1 = 1),
    #                                                     lower=c(0), 
    #                                                     upper=c(Inf), control = nls.lm.control(maxiter = 1024))
    #                    ,silent = TRUE)
    # 
    # 
    # one.fit.try <- try(one_compFit <- nls(y ~ exp(-(k1 * (x-Nmin*frameRate_sec))^beta), trace = F,
    #                                                     start=c(k1 = 1, beta = 1),
    #                                                     lower=c(0,0), 
    #                                                     upper=c(Inf,1), control = nls.lm.control(maxiter = 1024))
    #                    ,silent = TRUE)
    # 
    # 
    # 
    # one.new_y <- exp(-(coef(one_compFit)[[1]]* (x-Nmin*frameRate_sec))^coef(one_compFit)[[2]])
    # lines(x,one.new_y, col = "blue")
    
    
    if(class(two.fit.try) == "try-error" | 
       abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) < 0.001){
      cat(paste0("Two-component exponetnial decay cannot be fitted to the data!"), "\n")
      
      if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
        
        raw <- raw %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
        write.csv(raw, file.path(folder, paste0("1-CDF.raw_", save.as, ".csv")), row.names = FALSE)
        
        write.csv(plot.output.fin, file.path(folder, paste0("1-CDF.BS_", save.as, ".csv")), row.names = FALSE)
        write.csv(fit.value.summary.fin, file.path(folder, paste0("fit.value.summary_", save.as, ".csv")), row.names = FALSE)
        write.csv(fit.value.all.fin, file.path(folder, paste0("fit.value.all_", save.as, ".csv")), row.names = FALSE)
        
      }
      stop()
    }
    
    two.new_y <-  (coef(two_compFit)[[1]] *exp(-coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec)))
    
    
    ktb <- max(coef(two_compFit)[[2]], coef(two_compFit)[[3]])
    ksb <- min(coef(two_compFit)[[2]], coef(two_compFit)[[3]])
    
    if(coef(two_compFit)[[2]] >= coef(two_compFit)[[3]]){
      ftb <- coef(two_compFit)[[1]]
      fsb <- 1 - coef(two_compFit)[[1]]
    }else{
      ftb <- 1- coef(two_compFit)[[1]]
      fsb <- coef(two_compFit)[[1]]
    }
    
    lines(x,two.new_y, col = "blue")
    legend( x = "topright",
            legend = c("Two"),
            col = c("blue"), lwd = 2, lty = c(0),
            pch = c(19) )
    ####################################################################################################################################
    
    
    
    plot.data <- data.frame(x = x, y = y, rep = 1)
    plot.data_BS <- plot.data
    
    # k1_BS <- c(k1)
    # k2_BS <- c(k2)
    # k3_BS <- c(k3)
    # f1_BS <- c(f1)
    # f2_BS <- c(f2)
    # f3_BS <- c(f3)
    
    ktb_BS <- c(ktb)
    ksb_BS <- c(ksb)
    ftb_BS <- c(ftb)
    fsb_BS <- c(fsb)
    
    
    set.seed(123)
    for(i in 2:Boostrap){
      
      cat('\r', paste0("Boostrap: ", i))
      
      
      
      
      repeat{
        
        trackLength_BS <- trackLength[sample(1:length(trackLength), replace = TRUE)]
        
        my_ecdf <- ecdf(trackLength_BS) #It generate another function called my_ecdf
        trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength_BS + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength_BS)))))
        
        y = trackLengthCDF.$y
        x = trackLengthCDF.$x * frameRate_sec
        
        
        #raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
        raw <- data.frame(x = x, y = y) %>% distinct(x, y)
        x <- as.vector(raw$x)
        y <- as.vector(raw$y)
        
        # three.BSfit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
        #                                               #start=c(F1 = f1, F2 = f2, k1 = k1, k2 = k2, k3 = k3),
        #                                               start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
        #                                               lower=c(0,0,0,0,0),
        #                                               upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
        #                        ,silent = TRUE)
        
        
        
        y1 <- y[y >= cutoff]
        x1 <- x[1:length(y1)]
        if(weight){
          WEIGHT = 1/y1
          WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
        }
        
        if(weight != FALSE){
          two.BSfit.try <- try(two_compFit <- minpack.lm::nlsLM(y1 ~ (F1 * exp(-k1 * (x1-(Nmin)*frameRate_sec)) + (1 - F1) * exp(-k2 * (x1-(Nmin)*frameRate_sec))), trace = F,
                                                                #start=c(F1 = ftb, k1 = ktb, k2 = ksb),
                                                                start=c(F1 = 0.5, k1 = 1, k2 = 0.1),
                                                                lower=c(0,0,0),
                                                                upper=c(1, Inf, Inf), control = nls.lm.control(maxiter = 1024), weights = WEIGHT)
                               ,silent = TRUE)
        }else{
          two.BSfit.try <- try(two_compFit <- minpack.lm::nlsLM(y1 ~ (F1 * exp(-k1 * (x1-(Nmin)*frameRate_sec)) + (1 - F1) * exp(-k2 * (x1-(Nmin)*frameRate_sec))), trace = F,
                                                                #start=c(F1 = ftb, k1 = ktb, k2 = ksb),
                                                                start=c(F1 = 0.5, k1 = 1, k2 = 0.1),
                                                                lower=c(0,0,0),
                                                                upper=c(1, Inf, Inf), control = nls.lm.control(maxiter = 1024))
                               ,silent = TRUE)
        }
        
        
        
        
        if(#class(three.BSfit.try) != "try-error" &
          #(1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) > 0 &
          # coef(three_compFit)[[1]] > 0.2 &
          # coef(three_compFit)[[1]] < 1 &
          # coef(three_compFit)[[2]] > 0.1 &
          # coef(three_compFit)[[2]] < 1 &
          # coef(three_compFit)[[3]] > 0.05 &
          # coef(three_compFit)[[3]] < 1/(frameRate_sec/2) &
          # coef(three_compFit)[[4]] > 0.05 &
          # coef(three_compFit)[[4]] < 1/(frameRate_sec/2) &
          # coef(three_compFit)[[5]] > 0.05 &
          # coef(three_compFit)[[5]] < 1/(frameRate_sec/2) &
          # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) > 0.001 &
          # abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] > 0.001) &
          # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) > 0.001 &
          ###############
          class(two.BSfit.try) != "try-error" &
          abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) > 0.001
          ###############
          
        ){
          break
        }
        
        
      }
      
      # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) + 
      #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) + 
      #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
      
      
      # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
      # three.k.sort <- sort(three.k, decreasing = TRUE)
      # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
      # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
      # 
      # k1 <- three.k.sort[1]
      # k2 <- three.k.sort[2]
      # k3 <- three.k.sort[3]
      # f1 <- three.fraction.sort[1]
      # f2 <- three.fraction.sort[2]
      # f3 <- three.fraction.sort[3]
      # 
      # k1_BS <- c(k1_BS, k1)
      # k2_BS <- c(k2_BS, k2)
      # k3_BS <- c(k3_BS, k3)
      # f1_BS <- c(f1_BS, f1)
      # f2_BS <- c(f2_BS, f2)
      # f3_BS <- c(f3_BS, f3)
      
      ###############
      ktb <- max(coef(two_compFit)[[2]], coef(two_compFit)[[3]])
      ksb <- min(coef(two_compFit)[[2]], coef(two_compFit)[[3]])
      
      if(coef(two_compFit)[[2]] >= coef(two_compFit)[[3]]){
        ftb <- coef(two_compFit)[[1]]
        fsb <- 1 - coef(two_compFit)[[1]]
      }else{
        ftb <- 1- coef(two_compFit)[[1]]
        fsb <- coef(two_compFit)[[1]]
      }
      
      ktb_BS <- c(ktb_BS, ktb)
      ksb_BS <- c(ksb_BS, ksb)
      ftb_BS <- c(ftb_BS, ftb)
      fsb_BS <- c(fsb_BS, fsb)
      ###############
      
      plot.data_BS. <- data.frame(x = x, y = y, rep = i)
      plot.data_BS <- rbind(plot.data_BS, plot.data_BS.)
      
      
    }  # end of repeat loop
    cat("\n")
    #fit.value <- data.frame(k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, rep = seq_along(k1_BS))
    fit.value <- data.frame(#k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, 
      ktb = ktb_BS, ksb = ksb_BS, ftb = ftb_BS, fsb = fsb_BS, rep = seq_along(ktb_BS))
    
    
    # Pivot longer for facet
    fit.value.long <- fit.value %>% select(-rep) %>%
      tidyr::pivot_longer(everything(), names_to = c("fit.value"), values_to = "value")
    
    fit.value.long$rep  <-  rep(fit.value$rep, each = ncol(fit.value) - 1)
    
    
    fit.value.all <- fit.value %>% dplyr::mutate(N = length(trackLength), Nmin = Nmin, Nmax = Nmax)
    fit.value.all.fin <- rbind(fit.value.all.fin, fit.value.all)
    
    fit.value.summary <- fit.value.long %>% select(-rep) %>% dplyr::group_by(fit.value) %>% dplyr::summarise(mean=mean(value), sd=sd(value), boostrap = n(), N = length(trackLength), Nmin = Nmin, Nmax = Nmax)
    fit.value.summary.fin <- rbind(fit.value.summary.fin, fit.value.summary) %>% arrange (fit.value)
    
    
    
    
    plot.data_BS <-  plot.data_BS %>% group_by(rep) %>% distinct(x, y)
    plot.output <- plot.data_BS %>% dplyr::group_by(x) %>% dplyr::summarise(mean=mean(y), sd=sd(y))
    N.range <- data.frame(x = seq(Nmin*frameRate_sec, 1000*frameRate_sec, frameRate_sec))
    plot.output <- plot.output %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
    
    # F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "f1"]
    # F2 <- fit.value.summary$mean[fit.value.summary$fit.value == "f2"]
    # F3 <- fit.value.summary$mean[fit.value.summary$fit.value == "f3"]
    # k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "k1"]
    # k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "k2"]
    # k3 <- fit.value.summary$mean[fit.value.summary$fit.value == "k3"]
    # x <- plot.output$x
    # 
    # plot.output$three <-  F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + 
    #   F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + 
    #   F3 * exp(-k3 * (x-(Nmin)*frameRate_sec))
    ##################
    F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ftb"]
    k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ktb"]
    k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "ksb"]
    x <- plot.output$x
    
    plot.output$two <-  F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + 
      (1-F1) * exp(-k2 * (x-(Nmin)*frameRate_sec)) 
    ##################
    
    
    
    
    
    plot.output.fin <- rbind(plot.output.fin, plot.output)
    
    r <- ggplot(data = plot.output) + geom_point(aes(y=mean, x=x))+
      geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=x, fill = "band"), alpha = 0.3) +
      #geom_line(aes(x=x, y = three), color = "red") +
      geom_line(aes(x=x, y = two), color = "blue") +
      scale_fill_manual("",values="blue")+ 
      scale_x_continuous(trans = "log10", limits = c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100))) + scale_y_continuous(trans = "log10", limits = c(0.001,1)) +
      labs(x = "Time (s)", y = "1-CDF", title = paste0("Nmin: ", Nmin))  +
      guides(fill = "none") 
    print(r)
    
    
    # Plot histogram of fit.value
    q <- ggplot(data = fit.value.long, aes(x = value)) +
      geom_histogram() +
      geom_vline(data = fit.value.long[fit.value.long$rep == 1, ], aes(xintercept = value)) +
      facet_wrap(~fit.value, scales = "free") +
      labs(title = paste0("Nmin: ", Nmin)) 
    
    print(q)
    
    
    
    
  } # end of for-loop for Nmin
  cat("\n")
  
  if(length(Nmin) > 1){
    s <- ggplot(data = fit.value.summary.fin, aes(x = Nmin, y = mean)) +
      geom_line()+
      geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=Nmin, fill = "band"), alpha = 0.3) +
      facet_wrap(~fit.value, scales = "free") 
    print(s)
  }
  
  raw <- raw %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
  
  write.csv(plot.output.fin, file.path(folder, paste0("1-CDF.BS_", save.as, ".csv")), row.names = FALSE)
  write.csv(fit.value.summary.fin, file.path(folder, paste0("fit.value.summary_", save.as, ".csv")), row.names = FALSE)
  write.csv(fit.value.all.fin, file.path(folder, paste0("fit.value.all_", save.as, ".csv")), row.names = FALSE)
  write.csv(raw, file.path(folder, paste0("1-CDF.raw_", save.as, ".csv")), row.names = FALSE)
  
  return(fit.value.summary.fin)
} 

get2ExpoDecay(trackll.H2B.linked, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 10000, save.as = "H2B")
get2ExpoDecay(trackll.POI.linked, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 10000, save.as = "POI")

#############################################################################################################################################################################################################


#########################################   Step 7b: 2-expo fit, with Stretched correction  #################################################################################################################
# Curve fitting was weighted by 1/y. Weights were not applied to survival probabilities below 0.01 due to increased variability and data scarcity.

get2ExpoStretchedDecay <- function(folder, trackll, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 10000, save.as = "H2B", weight = TRUE, cutoff = 0.01, ktb_min = -Inf, ktb_max = Inf, 
                                   ksb_min = -Inf, ksb_max = Inf, beta_min = 0, beta_max = 1, beta.start = 1){
    
    TrackLength <- getTrackLength(trackll[[1]], frameRate_sec = frameRate_sec)
    N.range <- data.frame(x = seq(Nmin*frameRate_sec, 1000*frameRate_sec, frameRate_sec))
    
    fit.value.summary.fin <- c()
    plot.output.fin <- c()
    fit.value.all.fin <- c()
    for(Nmin in Nmin){
        trackLength <- dplyr::bind_rows(TrackLength)$TrackLength
        
        trackLength <- trackLength[trackLength >= Nmin & trackLength <= Nmax]
        
        cat(paste0("Number of trajectory: ", length(trackLength)), "\n")
        cat(paste0("Nmin: ", Nmin), "\n")
        
        my_ecdf <- ecdf(trackLength) #It generate another function called my_ecdf
        trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength)))))
        
        y = trackLengthCDF.$y
        x = trackLengthCDF.$x * frameRate_sec
        
        
        # raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
        RAW <- data.frame(x = x, y = y) %>% distinct(x, y)
        x <- as.vector(RAW$x)
        y <- as.vector(RAW$y) #/  (exp(- (0.206147945 * (x-(Nmin)*frameRate_sec))^0.749468478 ) )
        
        
        
        # data = data.frame(P, t)
        #weights = 1/y
        # weight[which(weight == Inf)] <- max(weight[which(weight != 
        #                                                    Inf)])
        
        
        
        
        ######################
        # three.fit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
        #                                             start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
        #                                             lower=c(0,0,0,0,0),
        #                                             upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
        #                      ,silent = TRUE)
        # 
        # 
        # if(class(three.fit.try) == "try-error" | 
        #    (1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) <= 0 | 
        #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) < 0.001 |
        #    abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] < 0.001) |
        #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) < 0.001){
        #   cat(paste0("Three-component exponetnial decay cannot be fitted to the data!"), "\n")
        #   
        #   if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
        #     
        #     raw <- raw %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
        #     write.csv(raw, file.path(getwd(), paste0("1-CDF.raw_", save.as, ".csv")), row.names = FALSE)
        #     
        #     write.csv(plot.output.fin, file.path(getwd(), paste0("1-CDF.BS_", save.as, ".csv")), row.names = FALSE)
        #     write.csv(fit.value.summary.fin, file.path(getwd(), paste0("fit.value.summary_", save.as, ".csv")), row.names = FALSE)
        #     write.csv(fit.value.all.fin, file.path(getwd(), paste0("fit.value.all_", save.as, ".csv")), row.names = FALSE)
        #     
        #   }
        #   stop()
        # }
        # 
        # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) + 
        #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) + 
        #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
        # 
        # 
        # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
        # three.k.sort <- sort(three.k, decreasing = TRUE)
        # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
        # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
        # 
        # k1 <- three.k.sort[1]
        # k2 <- three.k.sort[2]
        # k3 <- three.k.sort[3]
        # f1 <- three.fraction.sort[1]
        # f2 <- three.fraction.sort[2]
        # f3 <- three.fraction.sort[3]
        ######################
        
        plot(x,y, main= paste0("Nmin: ", Nmin), xlab="Time (s)", ylab="1-CDF", log = "xy", xlim=c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100)), ylim=c(0.001, 1))
        #lines(x,three.new_y, col = "red", log = "xy", add = TRUE)
        # legend( x = "topright",
        #         legend = c("Three"),
        #         col = c("red"), lwd = 2, lty = c(0),
        #         pch = c(19) )
        
        ####################################################################################################################################
        # two.fit.try <- try(two_compFit <- nlsLM(y ~ (   F1 * exp(-(k1 * (x-(Nmin)*frameRate_sec) + Kh2b) ) + (1 - F1) * exp(-(k2 * (x-(Nmin)*frameRate_sec) + Kh2b   )   )), trace = F,
        #                                         start=c(F1 = 0.3, k1 = 1, k2 = 0.008),
        #                                         lower=c(0,0,0),
        #                                         upper=c(1, Inf, Inf), control = nls.lm.control(maxiter = 1024))
        #                    ,silent = TRUE)
        
        
        y1 <- y#[y >= cutoff]
        x1 <- x[1:length(y1)]
        if(weight){
            WEIGHT = ifelse(y1 > cutoff, 1/y1, 1)
            WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
            WEIGHTS = "TRUE"
        }else{
            WEIGHTS = "FALSE" 
        }
        
        if(weight != FALSE){
            two.fit.try <- try(two_compFit <- nls(y1 ~ F1 * exp(-k1 * (x1-Nmin*frameRate_sec)) + (1 - F1) * exp(- (k2 * (x1-Nmin*frameRate_sec))^beta), trace = F,
                                                  start=list(F1 = 0.3, k1 = 1, k2 = 0.2, beta = beta.start),
                                                  lower=list(0,0,0, 0),
                                                  upper=list(1, Inf, Inf, 1),algorithm = "port", weights = WEIGHT)
                               ,silent = TRUE)
        }else{
            two.fit.try <- try(two_compFit <- nls(y1 ~ F1 * exp(-k1 * (x1-Nmin*frameRate_sec)) + (1 - F1) * exp(- (k2 * (x1-Nmin*frameRate_sec))^beta), trace = F,
                                                  start=list(F1 = 0.3, k1 = 1, k2 = 0.2, beta = beta.start),
                                                  lower=list(0,0,0, 0),
                                                  upper=list(1, Inf, Inf, 1),algorithm = "port")
                               ,silent = TRUE)
        }
        
        
        
        if(class(two.fit.try) == "try-error" | 
           abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) < 0.001){
            cat(paste0("Two-component exponetnial decay cannot be fitted to the data!"), "\n")
            
            if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
                
                write.csv(plot.output.fin, file.path(folder, paste0(save.as, "_1-CDF.BS", ".csv")), row.names = FALSE)
                write.csv(fit.value.summary.fin, file.path(folder, paste0(save.as, "_fit.value.summary", ".csv")), row.names = FALSE)
                write.csv(fit.value.all.fin, file.path(folder, paste0(save.as, "_fit.value.all", ".csv")), row.names = FALSE)
                write.csv(RAW, file.path(folder, paste0(save.as, "_1-CDF.raw", ".csv")), row.names = FALSE)
                
            }
            stop()
        }
        
        two.new_y <-  (coef(two_compFit)[[1]] *exp(-coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-(coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec))^coef(two_compFit)[[4]]))
        # two.new_y <-  coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec)) + (1 - coef(two_compFit)[[1]]) * exp(-  ((coef(two_compFit)[[3]]+Kh2b) * (x-(Nmin)*frameRate_sec))   )   
        # 
        #two.new_y <- coef(two_compFit)[[1]] * exp(-(coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec) + Kh2b) ) + (1 - coef(two_compFit)[[1]]) * exp(-(coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec) + Kh2b   )   )
        
        
        
        ktb <- coef(two_compFit)[[2]]
        ksb <- coef(two_compFit)[[3]]
        beta <- coef(two_compFit)[[4]]
        ftb <- coef(two_compFit)[[1]]
        fsb <- 1 - coef(two_compFit)[[1]]
        
        lines(x,two.new_y, col = "blue")
        legend( x = "topright",
                legend = c("Two"),
                col = c("blue"), lwd = 2, lty = c(0),
                pch = c(19) )
        
        RAW <- RAW %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax, cutoff = cutoff)
        RAW$fit <- (coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (RAW$x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-  (((coef(two_compFit)[[3]] * (RAW$x-(Nmin)*frameRate_sec))^coef(two_compFit)[[4]])))  )
        RAW$AIC <- AIC(two_compFit)
        RAW$BIC <- BIC(two_compFit)
        
        ####################################################################################################################################
        
        # Calculate AIC and BIC
        aic <- AIC(two_compFit)
        bic <- BIC(two_compFit)
        # residuals_info <- nlstools::nlsResiduals(two_compFit)
        
        plot.data <- data.frame(x = x, y = y)
        plot.data <- plot.data %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(rep = 1, fit = RAW$fit)
        plot.data_BS <- plot.data
        
        # k1_BS <- c(k1)
        # k2_BS <- c(k2)
        # k3_BS <- c(k3)
        # f1_BS <- c(f1)
        # f2_BS <- c(f2)
        # f3_BS <- c(f3)
        
        ktb_BS <- c(ktb)
        ksb_BS <- c(ksb)
        ftb_BS <- c(ftb)
        fsb_BS <- c(fsb)
        beta_BS <- c(beta)
        
        set.seed(123)
        for(i in 2:Boostrap){
            
            cat('\r', paste0("Boostrap: ", i))
            
            
            
            
            repeat{
                
                trackLength_BS <- trackLength[sample(1:length(trackLength), replace = TRUE)]
                
                my_ecdf <- ecdf(trackLength_BS) #It generate another function called my_ecdf
                trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength_BS + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength_BS)))))
                
                y = trackLengthCDF.$y
                x = trackLengthCDF.$x * frameRate_sec
                
                
                #raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
                raw <- data.frame(x = x, y = y) %>% distinct(x, y)
                x <- as.vector(raw$x)
                y <- as.vector(raw$y)
                
                # three.BSfit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
                #                                               #start=c(F1 = f1, F2 = f2, k1 = k1, k2 = k2, k3 = k3),
                #                                               start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
                #                                               lower=c(0,0,0,0,0),
                #                                               upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
                #                        ,silent = TRUE)
                
                
                y1 <- y#[y >= cutoff]
                x1 <- x[1:length(y1)]
                if(weight){
                    WEIGHT = ifelse(y1 > cutoff, 1/y1, 1)
                    WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
                }
                
                if(weight != FALSE){
                    two.BSfit.try <- try(two_compFit <- nls(y1 ~ F1 * exp(-k1 * (x1-Nmin*frameRate_sec)) + (1 - F1) * exp(- (k2 * (x1-Nmin*frameRate_sec))^beta), trace = F,
                                                            start=list(F1 = 0.3, k1 = 1, k2 = 0.2, beta = beta.start),
                                                            lower=list(0,0,0, 0),
                                                            upper=list(1, Inf, Inf, 1),algorithm = "port", weights = WEIGHT)
                                         ,silent = TRUE)
                }else{
                    two.BSfit.try <- try(two_compFit <- nls(y1 ~ F1 * exp(-k1 * (x1-Nmin*frameRate_sec)) + (1 - F1) * exp(- (k2 * (x1-Nmin*frameRate_sec))^beta), trace = F,
                                                            start=list(F1 = 0.3, k1 = 1, k2 = 0.2, beta = beta.start),
                                                            lower=list(0,0,0, 0),
                                                            upper=list(1, Inf, Inf, 1),algorithm = "port")
                                         ,silent = TRUE)
                }
                
                
                
                if(#class(three.BSfit.try) != "try-error" &
                    #(1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) > 0 &
                    # coef(three_compFit)[[1]] > 0.2 &
                    # coef(three_compFit)[[1]] < 1 &
                    # coef(three_compFit)[[2]] > 0.1 &
                    # coef(three_compFit)[[2]] < 1 &
                    # coef(three_compFit)[[3]] > 0.05 &
                    # coef(three_compFit)[[3]] < 1/(frameRate_sec/2) &
                    # coef(three_compFit)[[4]] > 0.05 &
                    # coef(three_compFit)[[4]] < 1/(frameRate_sec/2) &
                    # coef(three_compFit)[[5]] > 0.05 &
                    # coef(three_compFit)[[5]] < 1/(frameRate_sec/2) &
                    # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) > 0.001 &
                    # abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] > 0.001) &
                    # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) > 0.001 &
                    ###############
                    class(two.BSfit.try) != "try-error" &
                    abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) > 0.001 &
                    coef(two_compFit)[[2]] >= ktb_min & 
                    coef(two_compFit)[[2]] <= ktb_max & 
                    coef(two_compFit)[[3]] >= ksb_min &
                    coef(two_compFit)[[3]] <= ksb_max &
                    coef(two_compFit)[[4]] >= beta_min &
                    coef(two_compFit)[[4]] <= beta_max 
                    
                    # coef(two_compFit)[[2]] > 0.5 & 
                    # coef(two_compFit)[[2]] < 1/(frameRate_sec) & 
                    # coef(two_compFit)[[3]] > 0 &
                    # coef(two_compFit)[[4]] < 1 &
                    # coef(two_compFit)[[4]] > 0
                    
                    ###############
                    
                ){
                    break
                }
                
                
            }
            
            # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) + 
            #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) + 
            #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
            
            
            # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
            # three.k.sort <- sort(three.k, decreasing = TRUE)
            # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
            # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
            # 
            # k1 <- three.k.sort[1]
            # k2 <- three.k.sort[2]
            # k3 <- three.k.sort[3]
            # f1 <- three.fraction.sort[1]
            # f2 <- three.fraction.sort[2]
            # f3 <- three.fraction.sort[3]
            # 
            # k1_BS <- c(k1_BS, k1)
            # k2_BS <- c(k2_BS, k2)
            # k3_BS <- c(k3_BS, k3)
            # f1_BS <- c(f1_BS, f1)
            # f2_BS <- c(f2_BS, f2)
            # f3_BS <- c(f3_BS, f3)
            
            ###############
            ktb <- coef(two_compFit)[[2]]
            ksb <- coef(two_compFit)[[3]]
            beta <- coef(two_compFit)[[4]]
            
            
            ftb <- coef(two_compFit)[[1]]
            fsb <- 1 - coef(two_compFit)[[1]]
            
            
            ktb_BS <- c(ktb_BS, ktb)
            ksb_BS <- c(ksb_BS, ksb)
            ftb_BS <- c(ftb_BS, ftb)
            fsb_BS <- c(fsb_BS, fsb)
            beta_BS <- c(beta_BS, beta)
            
            ###############
            
            plot.data_BS. <- data.frame(x = x, y = y, rep = i)
            
            plot.data_BS. <- plot.data_BS. %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(rep = i)
            
            plot.data_BS.$fit <- (coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (plot.data_BS.$x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-  (((coef(two_compFit)[[3]] * (plot.data_BS.$x-(Nmin)*frameRate_sec))^coef(two_compFit)[[4]])))  )
            
            plot.data_BS <- rbind(plot.data_BS, plot.data_BS.)
            
            
        }  # end of repeat loop
        cat("\n")
        #fit.value <- data.frame(k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, rep = seq_along(k1_BS))
        fit.value <- data.frame(#k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, 
            ktb = ktb_BS, ksb = ksb_BS, ftb = ftb_BS, fsb = fsb_BS, beta = beta_BS, rep = seq_along(ktb_BS))
        
        
        # Pivot longer for facet
        fit.value.long <- fit.value %>% dplyr::select(-rep) %>%
            tidyr::pivot_longer(everything(), names_to = c("fit.value"), values_to = "value")
        
        fit.value.long$rep  <-  rep(fit.value$rep, each = ncol(fit.value) - 1)
        
        
        fit.value.all <- fit.value %>% dplyr::mutate(N = length(trackLength), Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS)
        fit.value.all.fin <- rbind(fit.value.all.fin, fit.value.all)
        
        fit.value.summary <- fit.value.long %>% dplyr::select(-rep) %>% dplyr::group_by(fit.value) %>% dplyr::summarise(mean=mean(value), sd=sd(value), boostrap = n(), N = length(trackLength), Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS)
        fit.value.summary.fin <- rbind(fit.value.summary.fin, fit.value.summary) %>% arrange (fit.value)
        
        
        
        plot.data_BS <-  plot.data_BS %>% dplyr::group_by(rep) 
        plot.output <- plot.data_BS %>% dplyr::group_by(x) %>% dplyr::summarise(mean=mean(y, na.rm = TRUE), sd=sd(y, na.rm = TRUE), fit.mean=mean(fit, na.rm = TRUE), fit.sd=sd(fit, na.rm = TRUE))
        N.range <- data.frame(x = seq(Nmin*frameRate_sec, 1000*frameRate_sec, frameRate_sec))
        plot.output <- plot.output %>% dplyr::mutate(Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS)
        
        
        # F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "f1"]
        # F2 <- fit.value.summary$mean[fit.value.summary$fit.value == "f2"]
        # F3 <- fit.value.summary$mean[fit.value.summary$fit.value == "f3"]
        # k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "k1"]
        # k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "k2"]
        # k3 <- fit.value.summary$mean[fit.value.summary$fit.value == "k3"]
        # x <- plot.output$x
        # 
        # plot.output$three <-  F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + 
        #   F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + 
        #   F3 * exp(-k3 * (x-(Nmin)*frameRate_sec))
        ##################
        F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ftb"]
        k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ktb"]
        k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "ksb"]
        beta <- fit.value.summary$mean[fit.value.summary$fit.value == "beta"]
        x <- plot.output$x
        #Kh2b <- (0.194058883 * ((x-(Nmin)*frameRate_sec))        )         ^0.785879326
        plot.output$fit.avg <-  F1 * exp(-k1 * (x-Nmin*frameRate_sec)) + (1 - F1) * exp(- (k2 * (x-Nmin*frameRate_sec))^beta)
        ##################
        
        plot.output$AIC <- aic
        plot.output$BIC <- bic
        
        
        plot.output.fin <- rbind(plot.output.fin, plot.output)
        
        r <- ggplot(data = plot.output) + geom_point(aes(y=mean, x=x))+
            geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=x, fill = "band"), alpha = 0.3) +
            #geom_line(aes(x=x, y = three), color = "red") +
            geom_line(aes(x=x, y = fit.mean), color = "blue") +
            scale_fill_manual("",values="blue")+ 
            scale_x_continuous(trans = "log10", limits = c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100))) + scale_y_continuous(trans = "log10", limits = c(0.001,1)) +
            labs(x = "Time (s)", y = "1-CDF", title = paste0("Nmin: ", Nmin))  +
            guides(fill = "none") 
        print(r)
        
        
        # Plot histogram of fit.value
        q <- ggplot(data = fit.value.long, aes(x = value)) +
            geom_histogram() +
            geom_vline(data = fit.value.long[fit.value.long$rep == 1, ], aes(xintercept = value)) +
            geom_vline(data = fit.value.summary.fin, aes(xintercept = mean), color = "red") +
            facet_wrap(~fit.value, scales = "free") +
            labs(title = paste0("Nmin: ", Nmin)) 
        
        print(q)
        
        
        
        
    } # end of for-loop for Nmin
    cat("\n")
    
    if(length(Nmin) > 1){
        s <- ggplot(data = fit.value.summary.fin, aes(x = Nmin, y = mean)) +
            geom_line()+
            geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=Nmin, fill = "band"), alpha = 0.3) +
            facet_wrap(~fit.value, scales = "free") 
        print(s)
    }
    
    
    write.csv(plot.output.fin, file.path(folder, paste0(save.as, "_1-CDF.BS", ".csv")), row.names = FALSE)
    write.csv(fit.value.summary.fin, file.path(folder, paste0(save.as, "_fit.value.summary", ".csv")), row.names = FALSE)
    write.csv(fit.value.all.fin, file.path(folder, paste0(save.as, "_fit.value.all", ".csv")), row.names = FALSE)
    write.csv(RAW, file.path(folder, paste0(save.as, "_1-CDF.raw", ".csv")), row.names = FALSE)
    
    return(fit.value.summary.fin)
    
}

get2ExpoCorrectedDecay <- function(folder, trackll, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 10000, save.as = "POI", weight = TRUE, a,b, cutoff = 0.01, 
                                   ktb_min = -Inf, ktb_max = Inf, ksb_min = -Inf, ksb_max = Inf){
    
    TrackLength <- getTrackLength(trackll[[1]], frameRate_sec = 250/1000)
    N.range <- data.frame(x = seq(Nmin*frameRate_sec, 1000*frameRate_sec, frameRate_sec))
    
    fit.value.summary.fin <- c()
    plot.output.fin <- c()
    fit.value.all.fin <- c()
    for(Nmin in Nmin){
        trackLength <- dplyr::bind_rows(TrackLength)$TrackLength
        
        trackLength <- trackLength[trackLength >= Nmin & trackLength <= Nmax]
        
        cat(paste0("Number of trajectory: ", length(trackLength)), "\n")
        cat(paste0("Nmin: ", Nmin), "\n")
        
        my_ecdf <- ecdf(trackLength) #It generate another function called my_ecdf
        trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength)))))
        
        y = trackLengthCDF.$y
        x = trackLengthCDF.$x * frameRate_sec
        
        
        # raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
        RAW <- data.frame(x = x, y = y) %>% distinct(x, y)
        x <- as.vector(RAW$x)
        y <- as.vector(RAW$y) #/  (exp(- (0.206147945 * (x-(Nmin)*frameRate_sec))^0.749468478 ) )
        
        # data = data.frame(P, t)
        #weights = 1/y
        # weight[which(weight == Inf)] <- max(weight[which(weight != 
        #                                                    Inf)])
        
        
        
        
        ######################
        # three.fit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
        #                                             start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
        #                                             lower=c(0,0,0,0,0),
        #                                             upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
        #                      ,silent = TRUE)
        # 
        # 
        # if(class(three.fit.try) == "try-error" | 
        #    (1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) <= 0 | 
        #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) < 0.001 |
        #    abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] < 0.001) |
        #    abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) < 0.001){
        #   cat(paste0("Three-component exponetnial decay cannot be fitted to the data!"), "\n")
        #   
        #   if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
        #     
        #     raw <- raw %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
        #     write.csv(raw, file.path(getwd(), paste0("1-CDF.raw_", save.as, ".csv")), row.names = FALSE)
        #     
        #     write.csv(plot.output.fin, file.path(getwd(), paste0("1-CDF.BS_", save.as, ".csv")), row.names = FALSE)
        #     write.csv(fit.value.summary.fin, file.path(getwd(), paste0("fit.value.summary_", save.as, ".csv")), row.names = FALSE)
        #     write.csv(fit.value.all.fin, file.path(getwd(), paste0("fit.value.all_", save.as, ".csv")), row.names = FALSE)
        #     
        #   }
        #   stop()
        # }
        # 
        # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) + 
        #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) + 
        #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
        # 
        # 
        # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
        # three.k.sort <- sort(three.k, decreasing = TRUE)
        # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
        # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
        # 
        # k1 <- three.k.sort[1]
        # k2 <- three.k.sort[2]
        # k3 <- three.k.sort[3]
        # f1 <- three.fraction.sort[1]
        # f2 <- three.fraction.sort[2]
        # f3 <- three.fraction.sort[3]
        ######################
        
        plot(x,y, main= paste0("Nmin: ", Nmin), xlab="Time (s)", ylab="1-CDF", log = "xy", xlim=c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100)), ylim=c(0.001, 1))
        #lines(x,three.new_y, col = "red", log = "xy", add = TRUE)
        # legend( x = "topright",
        #         legend = c("Three"),
        #         col = c("red"), lwd = 2, lty = c(0),
        #         pch = c(19) )
        Kh2b <- (a * ((x-(Nmin)*frameRate_sec))        )         ^b
        ####################################################################################################################################
        y1 <- y#[y >= cutoff]
        x1 <- x[1:length(y1)]
        Kh2b1 <- (a * ((x1-(Nmin)*frameRate_sec))        )         ^b
        
        if(weight){
            WEIGHT = ifelse(y1 > cutoff, 1/y1, 1)
            WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
            WEIGHTS = "TRUE"
        }else{
            WEIGHTS = "FALSE"
        }
        
        if(weight != FALSE){
            two.fit.try <- try(two_compFit <- nls(y1 ~ (   F1 * exp(-(k1 * (x1-(Nmin)*frameRate_sec) ) ) + (1 - F1) * exp(-(k2 * (x1-(Nmin)*frameRate_sec) + Kh2b1   )   )), trace = F,
                                                  start=c(F1 = 0.3, k1 = 1, k2 = 0.008),
                                                  lower=c(0,0,-Inf), 
                                                  upper=c(1, Inf, Inf), ,algorithm = "port", weights = WEIGHT)
                               ,silent = TRUE)
        }else{
            two.fit.try <- try(two_compFit <- nls(y1 ~ (   F1 * exp(-(k1 * (x1-(Nmin)*frameRate_sec) ) ) + (1 - F1) * exp(-(k2 * (x1-(Nmin)*frameRate_sec) + Kh2b1   )   )), trace = F,
                                                  start=c(F1 = 0.3, k1 = 1, k2 = 0.008),
                                                  lower=c(0,0,-Inf), 
                                                  upper=c(1, Inf, Inf), ,algorithm = "port")
                               ,silent = TRUE)
            
            
        }
        
        
        
        
        # two.fit.try <- try(two_compFit <- nlsLM(y ~ A + ( exp(-  (k1 * (x-(Nmin)*frameRate_sec))^beta   )   ), trace = F,
        #                                                     start=c(k1 = 0.1, beta = 1, A = 0.6),
        #                                                     lower=c(0,0, 0),
        #                                                     upper=c(Inf,  1, 1), control = nls.lm.control(maxiter = 1024))
        #                    ,silent = TRUE)
        
        
        
        
        if(class(two.fit.try) == "try-error" | 
           abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) < 0.001){
            cat(paste0("Two-component exponetnial decay cannot be fitted to the data!"), "\n")
            
            if(exists("plot.output.fin") & exists("fit.value.summary.fin") & exists("fit.value.all.fin") & exists("raw")){
                
                RAW <- RAW %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax)
                write.csv(plot.output.fin, file.path(folder, paste0(save.as, "_1-CDF.BS", ".csv")), row.names = FALSE)
                write.csv(fit.value.summary.fin, file.path(folder, paste0(save.as, "_fit.value.summary", ".csv")), row.names = FALSE)
                write.csv(fit.value.all.fin, file.path(folder, paste0(save.as, "_fit.value.all", ".csv")), row.names = FALSE)
                write.csv(RAW, file.path(folder, paste0(save.as, "_1-CDF.raw", ".csv")), row.names = FALSE)
                
            }
            stop()
        }
        
        # two.new_y <-  (coef(two_compFit)[[1]] *exp(-coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-(coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec))^coef(two_compFit)[[4]]))
        # two.new_y <-  coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec)) + (1 - coef(two_compFit)[[1]]) * exp(-  ((coef(two_compFit)[[3]]+Kh2b) * (x-(Nmin)*frameRate_sec))   )   
        # 
        two.new_y <- coef(two_compFit)[[1]] * exp(-(coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec) ) ) + (1 - coef(two_compFit)[[1]]) * exp(-(coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec) + Kh2b   )   )
        
        
        
        ktb <- coef(two_compFit)[[2]]
        ksb <- coef(two_compFit)[[3]]
        #beta <- coef(two_compFit)[[4]]
        
        ftb <- coef(two_compFit)[[1]]
        fsb <- 1 - coef(two_compFit)[[1]]
        
        lines(x,two.new_y, col = "blue")
        legend( x = "topright",
                legend = c("Two"),
                col = c("blue"), lwd = 2, lty = c(0),
                pch = c(19) )
        
        RAW <- RAW %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(Nmin = Nmin, Nmax = Nmax, cutoff = cutoff)
        KH2B <- (a * ((RAW$x-(Nmin)*frameRate_sec))        )         ^b
        RAW$fit <- (coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (RAW$x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-  (((coef(two_compFit)[[3]] * (RAW$x-(Nmin)*frameRate_sec))) + KH2B))  )
        RAW$AIC <- AIC(two_compFit)
        RAW$BIC <- BIC(two_compFit)
        
        ####################################################################################################################################
        
        # Calculate AIC and BIC
        aic <- AIC(two_compFit)
        bic <- BIC(two_compFit)
        # residuals_info <- nlstools::nlsResiduals(two_compFit)
        
        
        plot.data <- data.frame(x = x, y = y)
        plot.data <- plot.data %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(rep = 1, fit = RAW$fit)
        plot.data_BS <- plot.data
        
        # k1_BS <- c(k1)
        # k2_BS <- c(k2)
        # k3_BS <- c(k3)
        # f1_BS <- c(f1)
        # f2_BS <- c(f2)
        # f3_BS <- c(f3)
        
        ktb_BS <- c(ktb)
        ksb_BS <- c(ksb)
        ftb_BS <- c(ftb)
        fsb_BS <- c(fsb)
        #beta_BS <- c(beta)
        
        set.seed(123)
        for(i in 2:Boostrap){
            
            cat('\r', paste0("Boostrap: ", i))
            
            
            
            
            repeat{
                
                trackLength_BS <- trackLength[sample(1:length(trackLength), replace = TRUE)]
                
                my_ecdf <- ecdf(trackLength_BS) #It generate another function called my_ecdf
                trackLengthCDF. <- data.frame(x = c(sort(c(Nmin,trackLength_BS + 1))), y = c(1 - my_ecdf(sort(c(Nmin-1,trackLength_BS)))))
                
                y = trackLengthCDF.$y
                x = trackLengthCDF.$x * frameRate_sec
                
                
                #raw <- data.frame(x = x, y = y) %>% filter(x <= Nmax* frameRate_sec) %>% distinct(x, y)
                raw <- data.frame(x = x, y = y) %>% distinct(x, y)
                x <- as.vector(raw$x)
                y <- as.vector(raw$y)
                Kh2b <- (a * ((x-(Nmin)*frameRate_sec))        )         ^b
                
                # three.BSfit.try <- try(three_compFit <- nlsLM(y ~ (F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + (1 - F1 - F2) * exp(-k3 * (x-(Nmin)*frameRate_sec))), trace = F,
                #                                               #start=c(F1 = f1, F2 = f2, k1 = k1, k2 = k2, k3 = k3),
                #                                               start=c(F1 = 0.3, F2 = 0.5, k1 = 10, k2 = 1, k3 = 0.1),
                #                                               lower=c(0,0,0,0,0),
                #                                               upper=c(1, 1, Inf,Inf, Inf), control = nls.lm.control(maxiter = 1024))
                #                        ,silent = TRUE)
                
                y1 <- y#[y >= cutoff]
                x1 <- x[1:length(y1)]
                Kh2b1 <- (a * ((x1-(Nmin)*frameRate_sec))        )         ^b
                
                if(weight){
                    WEIGHT = ifelse(y1 > cutoff, 1/y1, 1)
                    WEIGHT[which(WEIGHT == Inf)] <- max(WEIGHT[which(WEIGHT !=Inf)])
                }
                
                if(weight != FALSE){
                    
                    two.BSfit.try <- try(two_compFit <- nls(y1 ~ (   F1 * exp(-(k1 * (x1-(Nmin)*frameRate_sec) ) ) + (1 - F1) * exp(-(k2 * (x1-(Nmin)*frameRate_sec) + Kh2b1   )   )), trace = F,
                                                            #start=c(F1 = ftb, k1 = ktb, k2 = ksb),
                                                            start=c(F1 = 0.3, k1 = 1, k2 = 0.008),
                                                            lower=c(0,0,-Inf),
                                                            upper=c(1, Inf, Inf), algorithm = "port", weight = WEIGHT)
                                         ,silent = TRUE)
                }else{
                    two.BSfit.try <- try(two_compFit <- nls(y1 ~ (   F1 * exp(-(k1 * (x1-(Nmin)*frameRate_sec) ) ) + (1 - F1) * exp(-(k2 * (x1-(Nmin)*frameRate_sec) + Kh2b1   )   )), trace = F,
                                                            #start=c(F1 = ftb, k1 = ktb, k2 = ksb),
                                                            start=c(F1 = 0.3, k1 = 1, k2 = 0.008),
                                                            lower=c(0,0,-Inf),
                                                            upper=c(1, Inf, Inf), algorithm = "port")
                                         ,silent = TRUE)
                }
                
                
                
                if(#class(three.BSfit.try) != "try-error" &
                    #(1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]]) > 0 &
                    # coef(three_compFit)[[1]] > 0.2 &
                    # coef(three_compFit)[[1]] < 1 &
                    # coef(three_compFit)[[2]] > 0.1 &
                    # coef(three_compFit)[[2]] < 1 &
                    # coef(three_compFit)[[3]] > 0.05 &
                    # coef(three_compFit)[[3]] < 1/(frameRate_sec/2) &
                    # coef(three_compFit)[[4]] > 0.05 &
                    # coef(three_compFit)[[4]] < 1/(frameRate_sec/2) &
                    # coef(three_compFit)[[5]] > 0.05 &
                    # coef(three_compFit)[[5]] < 1/(frameRate_sec/2) &
                    # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[4]]) > 0.001 &
                    # abs(coef(three_compFit)[[4]] - coef(three_compFit)[[5]] > 0.001) &
                    # abs(coef(three_compFit)[[3]] - coef(three_compFit)[[5]]) > 0.001 &
                    ###############
                    class(two.BSfit.try) != "try-error" &
                    abs(coef(two_compFit)[[2]] - coef(two_compFit)[[3]]) > 0.001 &
                    coef(two_compFit)[[2]] >= ktb_min & 
                    coef(two_compFit)[[2]] <= ktb_max & 
                    coef(two_compFit)[[3]] >= ksb_min &
                    coef(two_compFit)[[3]] <= ksb_max 
                    ###############
                    
                ){
                    break
                }
                
                
            }
            # Kh2b <- (0.194058883 * ((x-(Nmin)*frameRate_sec))        )         ^0.785879326
            # plot(x,y, main= paste0("Nmin: ", Nmin), xlab="Time (s)", ylab="1-CDF", log = "xy", xlim=c(0.25, 100), ylim=c(0.001, 1))
            # two.new_y <- coef(two_compFit)[[1]] * exp(-(coef(two_compFit)[[2]] * (x-(Nmin)*frameRate_sec) + Kh2b) ) + (1 - coef(two_compFit)[[1]]) * exp(-(coef(two_compFit)[[3]] * (x-(Nmin)*frameRate_sec) + Kh2b   )   )
            # lines(x, two.new_y)
            
            
            # three.new_y <-  coef(three_compFit)[[1]] * exp(-coef(three_compFit)[[3]] * (x-(Nmin)*frameRate_sec)) + 
            #   coef(three_compFit)[[2]] * exp(-coef(three_compFit)[[4]] * (x-(Nmin)*frameRate_sec)) + 
            #   (1-coef(three_compFit)[[1]]-coef(three_compFit)[[2]]) * exp(-coef(three_compFit)[[5]] * (x-(Nmin)*frameRate_sec))
            
            
            # three.k <- c(coef(three_compFit)[[3]], coef(three_compFit)[[4]], coef(three_compFit)[[5]])
            # three.k.sort <- sort(three.k, decreasing = TRUE)
            # three.fraction <- c(coef(three_compFit)[[1]], coef(three_compFit)[[2]], 1 - coef(three_compFit)[[1]] - coef(three_compFit)[[2]])
            # three.fraction.sort <- three.fraction[order(three.k, decreasing = TRUE)]
            # 
            # k1 <- three.k.sort[1]
            # k2 <- three.k.sort[2]
            # k3 <- three.k.sort[3]
            # f1 <- three.fraction.sort[1]
            # f2 <- three.fraction.sort[2]
            # f3 <- three.fraction.sort[3]
            # 
            # k1_BS <- c(k1_BS, k1)
            # k2_BS <- c(k2_BS, k2)
            # k3_BS <- c(k3_BS, k3)
            # f1_BS <- c(f1_BS, f1)
            # f2_BS <- c(f2_BS, f2)
            # f3_BS <- c(f3_BS, f3)
            
            ###############
            ktb <- coef(two_compFit)[[2]]
            ksb <- coef(two_compFit)[[3]]
            #beta <- coef(two_compFit)[[4]]
            
            
            ftb <- coef(two_compFit)[[1]]
            fsb <- 1 - coef(two_compFit)[[1]]
            
            
            ktb_BS <- c(ktb_BS, ktb)
            ksb_BS <- c(ksb_BS, ksb)
            ftb_BS <- c(ftb_BS, ftb)
            fsb_BS <- c(fsb_BS, fsb)
            #beta_BS <- c(beta_BS, beta)
            
            ###############
            
            plot.data_BS. <- data.frame(x = x, y = y, rep = i)
            
            plot.data_BS. <- plot.data_BS. %>% dplyr::full_join(N.range, by = "x") %>% arrange(x) %>% mutate(rep = i)
            
            KH2B <- (a * ((plot.data_BS.$x-(Nmin)*frameRate_sec))        )         ^b
            
            plot.data_BS.$fit <- (coef(two_compFit)[[1]] * exp(-coef(two_compFit)[[2]] * (plot.data_BS.$x-(Nmin)*frameRate_sec)) + (1-coef(two_compFit)[[1]]) * exp(-  (((coef(two_compFit)[[3]] * (plot.data_BS.$x-(Nmin)*frameRate_sec))) + KH2B))  )
            
            plot.data_BS <- rbind(plot.data_BS, plot.data_BS.)
            
            
            
            
            
        }  # end of for loop
        cat("\n")
        #fit.value <- data.frame(k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, rep = seq_along(k1_BS))
        fit.value <- data.frame(#k1 = k1_BS, k2 = k2_BS, k3 = k3_BS, f1 = f1_BS, f2 = f2_BS, f3 = f3_BS, 
            ktb = ktb_BS, ksb = ksb_BS, ftb = ftb_BS, fsb = fsb_BS,rep = seq_along(ktb_BS))   #beta = beta_BS, 
        
        
        # Pivot longer for facet
        fit.value.long <- fit.value %>% dplyr::select(-rep) %>%
            tidyr::pivot_longer(everything(), names_to = c("fit.value"), values_to = "value")
        
        fit.value.long$rep  <-  rep(fit.value$rep, each = ncol(fit.value) - 1)
        
        
        fit.value.all <- fit.value %>% dplyr::mutate(N = length(trackLength), Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS)
        fit.value.all.fin <- rbind(fit.value.all.fin, fit.value.all)
        
        fit.value.summary <- fit.value.long %>% 
            dplyr::select(-rep) %>% 
            dplyr::group_by(fit.value) %>% 
            dplyr::summarise(mean=mean(value), mean.sd=sd(value), boostrap = n(), N = length(trackLength), Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS)
        
        
        fit.value.summary.fin <- rbind(fit.value.summary.fin, fit.value.summary) %>% arrange (fit.value)
        
        
        KSB <- c(fit.value.summary.fin[[3,2]], fit.value.summary.fin[[3,3]])
        EXPR1 <- expression(1/KSB)
        DF1 <- cbind(KSB)
        RES1 <- propagate::propagate(expr = EXPR1, data = DF1, type = "stat", do.sim = F, verbose = F)
        KTB <- c(fit.value.summary.fin[[4,2]], fit.value.summary.fin[[4,3]])
        EXPR2 <- expression(1/KTB)
        DF2 <- cbind(KTB)
        RES2 <- propagate::propagate(expr = EXPR2, data = DF2, type = "stat", do.sim = F, verbose = F)
        
        Tau_sb <- c("Tau_sb", RES1$prop[[1]], RES1$prop[[3]], Boostrap, length(trackLength), Nmin, Nmax, cutoff, WEIGHTS)
        Tau_tb <- c("Tau_tb", RES2$prop[[1]], RES2$prop[[3]], Boostrap, length(trackLength), Nmin, Nmax, cutoff, WEIGHTS)
        
        #Tau_sb <- c("Tau_sb", 1/fit.value.summary.fin[3,2], msm::deltamethod(g = ~1/x1, mean= fit.value.summary.fin[3,2], cov=fit.value.summary.fin[3,3]^2), Boostrap, length(trackLength), Nmin, Nmax)
        #Tau_tb <- c("Tau_tb", 1/fit.value.summary.fin[4,2], msm::deltamethod(g = ~1/x1, mean= fit.value.summary.fin[4,2], cov=fit.value.summary.fin[4,3]^2), Boostrap, length(trackLength), Nmin, Nmax)
        
        fit.value.summary.fin <- rbind(fit.value.summary.fin, unname(Tau_sb), unname(Tau_tb))
        
        
        plot.data_BS <-  plot.data_BS %>% dplyr::group_by(rep) 
        plot.output <- plot.data_BS %>% dplyr::group_by(x) %>% dplyr::summarise(mean=mean(y, na.rm = TRUE), sd=sd(y, na.rm = TRUE), fit.mean=mean(fit, na.rm = TRUE), fit.sd=sd(fit, na.rm = TRUE))
        N.range <- data.frame(x = seq(Nmin*frameRate_sec, 1000*frameRate_sec, frameRate_sec))
        plot.output <- plot.output %>% dplyr::full_join(N.range, by = "x") %>% dplyr::arrange(x) %>% dplyr::mutate(Nmin = Nmin, Nmax = Nmax, cutoff = cutoff, WEIGHTS = WEIGHTS,
                                                                                                                   ktb_min = ktb_min, ktb_max = ktb_max, ksb_min = -ksb_min, ksb_max = ksb_max, a = a, b = b)
        
        plot.output$ktb_min[plot.output$ktb_min == -Inf | plot.output$ktb_min == Inf] <- NA
        plot.output$ktb_max[plot.output$ktb_max == -Inf | plot.output$ktb_max == Inf] <- NA
        plot.output$ksb_min[plot.output$ksb_min == -Inf | plot.output$ksb_min == Inf] <- NA
        plot.output$ksb_max[plot.output$ksb_max == -Inf | plot.output$ksb_max == Inf] <- NA
        
        
        # F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "f1"]
        # F2 <- fit.value.summary$mean[fit.value.summary$fit.value == "f2"]
        # F3 <- fit.value.summary$mean[fit.value.summary$fit.value == "f3"]
        # k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "k1"]
        # k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "k2"]
        # k3 <- fit.value.summary$mean[fit.value.summary$fit.value == "k3"]
        # x <- plot.output$x
        # 
        # plot.output$three <-  F1 * exp(-k1 * (x-(Nmin)*frameRate_sec)) + 
        #   F2 * exp(-k2 * (x-(Nmin)*frameRate_sec)) + 
        #   F3 * exp(-k3 * (x-(Nmin)*frameRate_sec))
        ##################
        F1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ftb"]
        k1 <- fit.value.summary$mean[fit.value.summary$fit.value == "ktb"]
        k2 <- fit.value.summary$mean[fit.value.summary$fit.value == "ksb"]
        #beta <- fit.value.summary$mean[fit.value.summary$fit.value == "beta"]
        x <- plot.output$x
        Kh2b <- (a * ((x-(Nmin)*frameRate_sec))        )         ^b
        plot.output$fit.avg <-  F1 * exp(-(k1 * (x-(Nmin)*frameRate_sec)) ) + 
            (1-F1) * exp(- (k2 * (x-(Nmin)*frameRate_sec) + Kh2b) ) 
        ##################
        
        
        plot.output$AIC <- aic
        plot.output$BIC <- bic
        
        plot.output.fin <- rbind(plot.output.fin, plot.output)
        
        
        r <- ggplot(data = plot.output) + geom_point(aes(y=mean, x=x))+
            geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=x, fill = "band"), alpha = 0.3) +
            #geom_line(aes(x=x, y = three), color = "red") +
            geom_line(aes(x=x, y = fit.mean), color = "blue") +
            scale_fill_manual("",values="blue")+ 
            scale_x_continuous(trans = "log10", limits = c(Nmin*frameRate_sec, min(Nmax*frameRate_sec,100))) + scale_y_continuous(trans = "log10", limits = c(0.001,1)) +
            labs(x = "Time (s)", y = "1-CDF", title = paste0("Nmin: ", Nmin))  +
            guides(fill = "none") 
        print(r)
        
        
        # Plot histogram of fit.value
        q <- ggplot(data = fit.value.long, aes(x = value)) +
            geom_histogram() +
            geom_vline(data = fit.value.long[fit.value.long$rep == 1, ], aes(xintercept = value)) +
            facet_wrap(~fit.value, scales = "free") +
            labs(title = paste0("Nmin: ", Nmin)) 
        
        print(q)
        
        
        
        
    } # end of for-loop for Nmin
    cat("\n")
    
    if(length(Nmin) > 1){
        s <- ggplot(data = fit.value.summary.fin, aes(x = Nmin, y = mean)) +
            geom_line()+
            geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, x=Nmin, fill = "band"), alpha = 0.3) +
            facet_wrap(~fit.value, scales = "free") 
        print(s)
    }
    
    write.csv(plot.output.fin, file.path(folder, paste0(save.as, "_1-CDF.BS", ".csv")), row.names = FALSE)
    write.csv(fit.value.summary.fin, file.path(folder, paste0(save.as, "_fit.value.summary", ".csv")), row.names = FALSE)
    write.csv(fit.value.all.fin, file.path(folder, paste0(save.as, "_fit.value.all", ".csv")), row.names = FALSE)
    write.csv(RAW, file.path(folder, paste0(save.as, "_1-CDF.raw", ".csv")), row.names = FALSE)
    return(fit.value.summary.fin)
    
}


fit.value.H2B <- get2ExpoStretchedDecay(trackll.H2B.linked, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 5000)
fit.value.POI <- get2ExpoCorrectedDecay(trackll.POI.linked, frameRate_sec = 250/1000, Nmin = 1, Nmax = Inf, Boostrap = 5000, a = unlist(fit.value.H2B[4,2]), b = unlist(fit.value.H2B[1,2]))

#############################################################################################################################################################################################################






