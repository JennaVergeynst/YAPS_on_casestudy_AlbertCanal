rm(list=ls())
graphics.off()

## To update YAPS
#devtools::install_github("baktoft/yaps")

## For R-studio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set file directory as working directory

## Outside R-studio:
#setwd(getSrcDirectory()[1])

path <- getwd()
results_path <- path

library(yaps)
library("parallel")
library(dplyr)
library(data.table)
source('yaps_functions_for_fish_tags.R')

## Parameters
pingType = 'rbi' # min and max bi is checked inside function
hydro_data<- read.csv(paste0(path, '/../Receiver_location.csv'),row.names = 'station_name')
# Make the grid centred close to 0;0 (not absolutely necessary, but don't work with UTM coordinates)
X_shift <- 0
Y_shift <- 0

## Data
ID = '100'
data <- read.csv(paste(path,'/TOA_data_',toString(ID),'.csv',sep = ''),check.names = FALSE)
data$ID <- ID

## split along passings
passings_list <- split(data, data$groups_pas)

## Run one track
toa_list_el = passings_list$'0'
track_res_list <- prepare_and_run(toa_list_el, hydro_data, pingType, X_shift, Y_shift, results_path, time_col='synced_time')
track <- track_res_list[[1]]
residuals <- track_res_list[[2]]
colMeans(residuals, na.rm = TRUE) # average residual per receiver
rowMeans(residuals, na.rm = TRUE) # average residual per calculated position
X <- track$X
Y <- track$Y
plot(X,Y)
est_error <- sqrt(track$Xe**2+track$Ye**2)
mean(est_error)

## Run all tracks in parallel
ptm <- proc.time()
cl = makeCluster(4, outfile="") # Adapt to number of cores on your system
clusterExport(cl, list("prepare_and_run", "runYAPS", "runTrack", "getInp", "runTmb"))
track <- clusterApplyLB(cl, passings_list, prepare_and_run, hydro_data, pingType, 
                        X_shift, Y_shift, results_path, time_col='synced_time')
proc.time() - ptm
