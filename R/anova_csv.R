# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())
# Set working directory
setwd("C:/Users/bdour/Documents/Data/Processing/Biomechanical variables/Comparison/Vicon vs. Kinect")

input_file <- "raw_results.csv" 

# read input
data <- read.csv(file = input_file, header = TRUE)

# ICC
library(irr)
#   Trunk flexion
vicon_trunk_fl <- data[-1,3]
kinect_trunk_fl <- data[-1,4]
cbind(vicon_trunk_fl, kinect_trunk_fl)
icc_trunk <- icc(cbind(vicon_trunk_fl, kinect_trunk_fl), "twoway", "agreement", "single", conf.level = 0.95)
icc_trunk
# ANOVA
#fit <- aov(data[,1] ~ data[,2], data=data)




