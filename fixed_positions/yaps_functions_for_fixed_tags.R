#' Prepare and run yaps on a piece of toa-dataframe
prepare_and_run <- function(toa_list_el, hydro_data, pingType, rbi_min, rbi_max, X_shift, Y_shift, results_path, chunk_len, time_col='SyncTime', min_obs=3){
  
  library(dplyr)
  
  toa_on_recs <- toa_list_el %>% dplyr:: select(starts_with("ST"))
  receivers <- colnames(toa_on_recs)
  hydros <- hydro_data[rownames(hydro_data) %in% receivers, c('X','Y')]
  colnames(hydros) <- c('hx', 'hy')
  colnames(toa_on_recs) == rownames(hydros)
  
  track <- data.frame() # return empty if try runs into error
  try({
    toa_df <- toa_list_el[,colnames(toa_list_el)%in% receivers]
    toa_info <- toa_list_el[,colnames(toa_list_el)%in% receivers==FALSE]
    toa <- t(as.matrix(toa_df))
    track_and_residuals <- runYAPS(toa, toa_info, hydros, pingType, rbi_min, rbi_max, X_shift, Y_shift, time_col, min_obs)
    track <- track_and_residuals[[1]]
    residuals_df <- track_and_residuals[[2]]
    colnames(residuals_df) <- receivers
  })
  
  chunk_nb <- toa_list_el$chunks[[1]]
  groups_pas <- toa_list_el$groups_pas[[1]]
  ID <- toa_list_el$ID[[1]]
  
  write.csv(track,paste0(results_path,'/track_',toString(ID),'_pass_', toString(groups_pas), '_chunklen_', toString(chunk_len), '_chunk_', toString(chunk_nb), '.csv'))
  write.csv(residuals_df,paste0(results_path,'/residuals_',toString(ID),'_pass_', toString(groups_pas), '_chunklen_', toString(chunk_len), '_chunk_', toString(chunk_nb),'.csv'))
  
  if (nrow(track)==0){
    DATETIME <- toa_info[,time_col]
    write.csv(DATETIME, paste0(results_path, '/failed_period_',toString(ID),'_pass_', toString(groups_pas), '_chunklen_', toString(chunk_len), '_chunk_', toString(chunk_nb), '.csv'))
  }
  
  return(list(track, residuals_df))
}


#' Function to run yaps on a given toa matrix
#' @param toa TOA-matrix: matrix with receivers in rows and detections in columns.
#' @param toa_info Dataframe with datetime, soundspeed and groups_pas columns corresponding to toa
#' @param hydros Dataframe with columns hx and hy containing positions of the receivers. 
#' @param pingType Type of transmitter to simulate - either stable burst interval ('sbi') or random burst interval ('rbi')
#' @param rbi_min,rbi_max Minimum and maximum BI for random burst interval transmitters
#' @param X_shift,Y_shift int: shift of X,Y coordinates to the center for easier running yaps

runYAPS <- function(toa, toa_info, hydros, pingType, rbi_min, rbi_max, X_shift, Y_shift, time_col, min_obs=3){
  # Make the grid centred close to 0;0 (not absolutely necessary, but don't work with UTM coordinates)
  hydros[,'hx'] <- hydros[,'hx'] - X_shift
  hydros[,'hy'] <- hydros[,'hy'] - Y_shift
  
  DATETIME <- toa_info[,time_col]
  ss_measured <- toa_info$soundspeed
  chunk_nb <- toa_info$chunks[[1]]
  groups_pas <- toa_info$groups_pas[[1]]
  
  track <- data.frame() # return empty if try runs into error
  residuals_df <- data.frame() # return empty if try runs into error
  
  try({
    # YAPS is most sensible to few observation in head and tail
    # throw away observations in head and tail where nobs < 3 and nobs10 < 2
    # If nobs is overall too small, code below will return error
    nobs <- apply(toa, 2, function(k) sum(!is.na(k))) # notNA sum over columns
    nobs10 <- caTools::runmean(nobs, 10) # running mean over 10 elements
    first <- which(nobs10 >= 2 & nobs >= min_obs)[1]
    last <- rev(which(nobs10 >= 2 & nobs >= min_obs))[1]
    toa <- toa[,first:last]
    T0 <- min(toa, na.rm=TRUE) #reshape so first TOA is 0
    toa <- toa - T0
    DATETIME <- DATETIME[first:last]
    ss_measured <- ss_measured[first:last]
    
    attempt <- 1
    
    while(nrow(track)==0 && attempt <=5){
      attempt <- attempt+1
      try({
        outTmb_list <- runTrack(toa, hydros, pingType, rbi_min, rbi_max, X_shift, Y_shift, ss_measured)
        outTmb <- outTmb_list[[1]]
        outTmb$DATETIME <- DATETIME
        outTmb$ss_measured <- ss_measured
        outTmb$pas_nr <- groups_pas
        outTmb$chunk_nr <- chunk_nb
        track <- rbind(track,outTmb)
        
        resids <- outTmb_list[[2]]
        residuals_t <- t(data.frame(resids))
        residuals_df <- rbind(residuals_df, residuals_t)
      })
    }

  })# return empty track if error occurs
  
  return(list(track, residuals_df))
}

#' Proper yaps running
runTrack <- function(toa, hydros, pingType, rbi_min, rbi_max, X_shift, Y_shift, ss_data){
  # n_ss: number of soundspeed estimates. One estimate of SS per hour or so is enough (not needed to have one per ping).
  # sdInits=1 : use new inits each time
  if(pingType == 'sbi'){
    inp <- getInp(hydros, toa, E_dist="Mixture", n_ss=10, pingType=pingType, sdInits=1)
  } else if(pingType == 'rbi'){
    inp <- getInp(hydros, toa, E_dist="Mixture", n_ss=10, pingType=pingType, sdInits=1, rbi_min=rbi_min, rbi_max=rbi_max, ss_data_what = "data", ss_data)
  }
  
  pl <- c()
  maxIter <- ifelse(pingType=="sbi", 500, 5000)
  outTmb <- runTmb(inp, maxIter=maxIter, getPlsd=TRUE, getRep=TRUE)
  
  # Estimates in pl
  pl <- outTmb$pl
  
  # Error estimates in plsd
  plsd <- outTmb$plsd
  
  # add again center correction
  X <- pl$X + X_shift
  Y <- pl$Y + Y_shift
  
  Xe <- plsd$X
  Ye <- plsd$Y
  
  est_error <- sqrt(Xe**2+Ye**2)
  
  resids <- outTmb$rep$mu_toa - inp$datTmb$toa # difference between final estimation of toa (from estimated track), and observed toa
  
  track <- data.frame(X, Y, Xe, Ye, stringsAsFactors=FALSE) 
  if (mean(est_error)>2){ # return empty dataframe if estimated error is too large
    track <- data.frame()
  }
  
  return(list(track, resids))
}



#' @param toa_df TOA dataframe with receivers in columns and detections in rows (NO OTHER COLUMNS!)
#' @param window window of moving average 
plot_nobs <- function(toa_df, window=10){
  toa <- t(as.matrix(toa_df))
  T0 <- min(toa, na.rm=TRUE) #reshape so first TOA is 0
  toa <- toa - T0
  nobs <- apply(toa, 2, function(k) sum(!is.na(k))) # notNA sum over columns
  nobs_ma <- caTools::runmean(nobs, window) # running mean over 10 elements
  first <- which(nobs_ma >= 2 & nobs >= 3)[1]
  last <- rev(which(nobs_ma >= 2 & nobs >= 3))[1]
  plot(nobs)
  lines(nobs_ma, col="red")
  abline(v=c(first,last), col="blue")
}

#' split data in smaller chunks manageable by yaps
#' @param toa_data TOA dataframe with receivers and groups_pas in columns.
#' @param chunklen length of each chunk witin each groups_pas 
chunk_toa <- function(toa_data, chunklen){
  toa_data$chunks <- NA
  chunk_counter <- 1
  for (nr in as.list(unique(toa_data$groups_pas))){
    len <- length(toa_data[toa_data$groups_pas==nr,]$groups_pas)
    # len is minimal 120, since this was min_track_length in creation of TOA-data
    n <- ceiling(len/chunklen)
    # in case that it is needed to split up
    if (len > chunklen){
      chunks <- c()
      for (i in as.list(1:(n-1))){ # going until the before-last chunk, because the last chunk gets special treatment
        chunks <- append(chunks, c(rep(chunk_counter,chunklen))) # changed i by chunk_couter
        chunk_counter <- chunk_counter+1
      }
      if (len-chunklen*(n-1) > 0.5*chunklen){ # take the last chunk separate only if it's at least 1/2 of chunklen
        chunks <- append(chunks, c(rep(chunk_counter,chunklen))) # changed n by chunk_couter
        chunk_counter <- chunk_counter+1
      }
      else{ # add last chunk to the before-last is the last one is too short
        chunks <- append(chunks, c(rep(chunk_counter-1,chunklen))) # changed n by chunk_couter
        # here don't increase chunk_counter, because previous chunk_counter is used
      }

      chunks <- chunks[1:len]
      toa_data[toa_data$groups_pas==nr,]$chunks <- chunks
    }
    # when not needed to split up within the groups_pas
    else{
      toa_data[toa_data$groups_pas==nr,]$chunks <- chunk_counter
      chunk_counter <- chunk_counter+1
    }
  }
  return(toa_data)
}