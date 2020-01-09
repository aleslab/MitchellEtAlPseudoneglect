## Analysis of tactile line responses from Mitchell, Ales & Harris Pseudoneglect reliability study
## Code: A.G. Mitchell, 09.02.20

###### libraries and loading ######
library(readr)
library(ggplot2)
library(reshape)

# load data
dataPath <- 'M:/Alex_Files/Experiments/Bias/Data/tactile-rod/'
anaPath <- 'M:/Alex_Files/Experiments/Bias/Analysis/tactile-rod_Abi/'

setwd(dataPath)
filenames <- dir(
  dataPath, 
  recursive = TRUE, 
  full.names = FALSE, 
  pattern = '.csv') #filelist
# results file
res <- read.csv(text = 'size,shift,start,hand,response,SUB,SESS')

##### extracting data ######
# get key data from .csv (adding sub and session)
for (file in filenames){
  tmp <- read.csv(file)[, c(1:5)]
  tmp$SUB <- substr(file ,1,3)
  tmp$SESS <- substr(file ,6,6)
  res <- rbind(res,tmp)
}

# extracting rod physical midpoint - adding
res$MID <- factor(cut(res$size, 3), labels = c('50', '100', '150'))
res$MID <- as.numeric(as.character(res$MID))
res$LEN <- factor(cut(res$size, 3), labels = c('100', '200', '300'))
# column names for clarity
colnames(res)[colnames(res) == 'size'] <- 'LEN'
colnames(res)[colnames(res) == 'shift'] <- 'POS'
colnames(res)[colnames(res) == 'start'] <- 'START'
colnames(res)[colnames(res) == 'hand'] <- 'HAND'
colnames(res)[colnames(res) == 'response'] <- 'pMID'
# reorder data-frame
res <- res[, c(6,7,9,8,2:5)]
# calculating error
res$ERR <- res$pMID - res$MID #how to do this with mid as a factor

### continue from here! collapsing, calculating and plotting