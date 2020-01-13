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

##### Averaging data across factors #####
# aggregate by line length and position
res_means_pos <- aggregate(ERR ~ SESS*LEN*POS*SUB, mean, data = res)
res_means_pos <- res_means_pos[, c(4,1:3,5)] #putting subject first again
# mean across position
res_means <- aggregate(ERR ~ SESS*LEN*SUB, mean, data = res)
res_means <- res_means[, c(3,1,2,4)] #putting subject first again
# saving
setwd(anaPath)
write.csv(res_means_pos, 'tactile-rod_data.csv', row.names = FALSE)


# plot to have a look :)
ggplot(res_means, aes(x = ERR, y = LEN), colour = SESS) +
  geom_point(aes(colour = SESS)) + facet_wrap(~SUB) +
  theme_bw()

# mean of all participants
res_length <- aggregate(ERR ~ LEN*SUB, mean, data = res_means)
# saving
setwd(anaPath)
write.csv(res_length, 'tactile-rod_meanbylength.csv', row.names = FALSE)

##### plotting the data #####
dodge = position_dodge(0.2)
ggplot(res_length, aes(x = LEN, y = ERR)) +
  geom_point(shape = 1, size = 4.5) +
  geom_line(aes(group = SUB), alpha = .3, size = .85) +
  labs(title = 'Tactile rod bisection', x = 'Line length (mm)', 
       y = 'Bisection error (mm)', element_text(size = 12)) +
  geom_hline(yintercept = 0, size = 0.7) + ylim(-20,15) +
  theme_bw() 
# saving plot


