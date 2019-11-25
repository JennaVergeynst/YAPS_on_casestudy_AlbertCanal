rm(list=ls())
graphics.off()

## To update YAPS
#devtools::install_github("baktoft/yaps")

## For R-studio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set file directory as working directory

## Outside R-studio:
#setwd(getSrcDirectory()[1])

path <- getwd()

library(yaps)
library(dplyr)
library(data.table)

ID ='255' # or '54329' or '16200'

data <- read.csv(paste('TOA_data_',toString(ID),'.csv',sep = ''),check.names = FALSE)
data$ID <- ID
hydro_data<- read.csv(paste0(path, '/../Receiver_location.csv'),row.names = 'station_name')

toa_df <- data %>% dplyr:: select(starts_with("ST"))
receivers <- colnames(toa_df)
hydros <- hydro_data[rownames(hydro_data) %in% receivers, c('X','Y')] 
colnames(hydros) <- c('hx', 'hy')
colnames(toa_df) == rownames(hydros)

toa_all <- as.matrix(toa_df)
ss_all <- data$soundspeed

## Check where diff is too large
plot(diff(rowMeans(toa_all, na.rm=TRUE)))

nobs <- apply(toa_all, 1, function(k) sum(!is.na(k)))
plot(nobs)

hydros_yaps <- hydros[, c('hx','hy')]
# Make sure to have enough observations in the beginning and at the end
first_idx <- which(nobs >= 3)[1]
last_idx <- rev(which(nobs >= 3))[1]
toa_yaps <- toa_all[first_idx:last_idx, ]
dimnames(toa_yaps) <- NULL
ss_yaps <- ss_all[first_idx:last_idx]

T0 <- min(toa_yaps, na.rm=TRUE)
toa_yaps <- toa_yaps - T0

bi_range <- range(diff(rowMeans(toa_yaps, na.rm=TRUE)), na.rm=TRUE)
bi_range

inp <- yaps::getInp(hydros=hydros_yaps, toa=t(toa_yaps), E_dist="Mixture", n_ss=10, pingType='rbi', sdInits = 1, rbi_min =bi_range[1]-2 ,  rbi_max = bi_range[2]+2, ss_data_what = "data", ss_data = ss_yaps)
str(inp)

# might require a few attempts but should get there - remember to rerun getInp as well to get new random inits
outTmb <- runTmb(inp, maxIter = 1000, getPlsd = TRUE, getRep = TRUE)
X <- outTmb$pl$X
Y <- outTmb$pl$Y
plot(X,Y)

Xe <- outTmb$plsd$X
Ye <- outTmb$plsd$Y
error <- sqrt(Xe**2+Ye**2)
print('If error too big, run again with new initial parameters (inp).')
mean(error)

DATETIME <- data$synced_time[first_idx:last_idx]
ss_measured <- data$soundspeed[first_idx:last_idx]

track <- data.frame(X, Y, Xe, Ye, DATETIME, ss_measured, stringsAsFactors=FALSE) 
track$ID <- ID

residuals <- outTmb$rep$mu_toa - inp$datTmb$toa
residuals <- transpose(data.frame(residuals))
# colnames(residuals) <- receivers
# colMeans(residuals, na.rm=TRUE)

write.csv(track,paste0('track_',toString(ID),'.csv'))
write.csv(residuals,paste0('residuals_',toString(ID),'.csv'))


