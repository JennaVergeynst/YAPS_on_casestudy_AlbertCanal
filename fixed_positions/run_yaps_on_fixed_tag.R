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
source('yaps_functions_for_fixed_tags.R')

## Parameters
pingType = 'rbi'
mean_bi = 600
rbi_min = 540
rbi_max = 660
ID = '65043'
self = 'ST15'
chunk_len <- 2000

## Data
data <- read.csv(paste('TOA_data_sample_',toString(ID),'.csv',sep = ''),check.names = FALSE)

data$ID <- ID
hydro_data<- read.csv(paste0(path, '/../Receiver_location.csv'),row.names = 'station_name')

toa_df <- data %>% dplyr:: select(starts_with("ST"))
receivers <- colnames(toa_df)
hydros <- hydro_data[rownames(hydro_data) %in% receivers, c('X','Y')] 
colnames(hydros) <- c('hx', 'hy')
colnames(toa_df) == rownames(hydros)

toa_all <- as.matrix(toa_df)
ss_all <- data$soundspeed

## Sync tag is co-located with self receiver
focal_hydro   <- which(rownames(hydros) == self)
focal_toa_col <- which(colnames(toa_all) == self)

## estimating distances from the self based on hydro positions and data in TOA
hydros_dt <- data.table(hydros)
hydros_dt[, name:=rownames(hydros)]
hydros_dt[, dist2focal:=sqrt((hx-hx[focal_hydro])^2 + (hy-hy[focal_hydro])^2)]
hydros_dt[, dist_from_toa:=colMeans((toa_all- toa_all[,focal_toa_col])*1460, na.rm=TRUE)]
hydros_dt[, dist_diffs:= abs(dist2focal - dist_from_toa)]
hydros_dt

## Remove behind-the-corner receivers (in this case ST7, ST10, ST12) and self-receiver
corner_idx <- which(hydros_dt$dist_diffs>20)
corner_hydros <- hydros_dt$name[corner_idx]
data[, which(colnames(data) %in% c(corner_hydros,self))] <- NA

## Split data according to passings, and then each passing in chunks
data$chunks <- 0
## split along passings
passings_list <- split(data, data$groups_pas)
## add chunks to passings that are too long
for (i in 1:length(passings_list)){
  if (dim(passings_list[[i]])[1]>chunk_len){
    passings_list[[i]] <- chunk_toa(passings_list[[i]], chunk_len)
  }
}
## remake the dataframe
data_rebind <- bind_rows(passings_list)
## split according to the combo passing-chunk
chunk_list <- split(data_rebind, data_rebind[,c('groups_pas', 'chunks')])
## remove empty dataframes from list
remove <- sapply(chunk_list, function(i) dim(i)[1]==0)
chunk_list <- chunk_list[!remove]


# Make the grid centred close to 0;0 (not absolutely necessary, but don't work with UTM coordinates)
X_shift <- 0
Y_shift <- 0

toa_list_el = chunk_list$`0.1`

# Run yaps and write out the result and residuals (or the timeperiod if the run failed)
track_res_list <- prepare_and_run(toa_list_el, hydro_data, pingType, rbi_min, rbi_max, 
                                  X_shift, Y_shift, results_path, chunk_len, time_col='synced_time')

# Check track
track <- track_res_list[[1]]
residuals <- track_res_list[[2]]
colMeans(residuals, na.rm = TRUE) # average residual per receiver
rowMeans(residuals, na.rm = TRUE) # average residual per calculated position


## Run all tracks in parallel

ptm <- proc.time()
cl = makeCluster(3, outfile="")
clusterExport(cl, list("prepare_and_run", "runYAPS", "runTrack", "getInp", "runTmb"))
track <- clusterApplyLB(cl, chunk_list, prepare_and_run, hydro_data, pingType, 
                        rbi_min, rbi_max, X_shift, Y_shift, results_path, chunk_len, time_col='synced_time')
proc.time() - ptm
