# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

# Vicon vs. Kinect (ICC on critical values)

#   Set working directory
setwd("C:/Users/bdour/Documents/Work/Toronto/Sunnybrook/ACL Injury Screening/Data/Processed/Biomechanical variables/Comparison/Vicon vs. Kinect")

#   Data import
input_file <- "V-vs-K_critical_values.csv" 
data <- read.csv(file = input_file, header = TRUE)
length_data = dim(data)[1]-2
data <- data[1:length_data,]

#   Library import
library(irr)

#   Trunk flexion
#     ICC
vicon  <- as.numeric(matrix(data[-1,3]))
kinect <- as.numeric(matrix(data[-1,4]))
icc_trunk_fl <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_trunk_fl
#     R2
trunk_fl.lm <- lm(vicon ~ kinect)
r2_trunk_fl <- summary(trunk_fl.lm)$r.squared
r2_trunk_fl

#   Hip flexion

#   Right
#     ICC
vicon  <- as.numeric(matrix(data[-1,5]))
kinect <- as.numeric(matrix(data[-1,6]))
icc_right_hip_fl <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_right_hip_fl
#     R2
right_hip_fl.lm <- lm(vicon ~ kinect)
r2_right_hip_fl <- summary(right_hip_fl.lm)$r.squared
r2_right_hip_fl

#   Left
#     ICC
vicon  <- as.numeric(matrix(data[-1,7]))
kinect <- as.numeric(matrix(data[-1,8]))
icc_left_hip_fl <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_left_hip_fl
#     R2
left_hip_fl.lm <- lm(vicon ~ kinect)
r2_left_hip_fl <- summary(left_hip_fl.lm)$r.squared
r2_left_hip_fl

#   Knee flexion

#   Right
#     ICC
vicon  <- as.numeric(matrix(data[-1,9]))
kinect <- as.numeric(matrix(data[-1,10]))
icc_right_knee_fl <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_right_knee_fl
#     R2
right_knee_fl.lm <- lm(vicon ~ kinect)
r2_right_knee_fl <- summary(right_knee_fl.lm)$r.squared
r2_right_knee_fl

#   Left
#     ICC
vicon  <- as.numeric(matrix(data[-1,11]))
kinect <- as.numeric(matrix(data[-1,12]))
icc_left_knee_fl <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_left_knee_fl
#     R2
left_knee_fl.lm <- lm(vicon ~ kinect)
r2_left_knee_fl <- summary(left_knee_fl.lm)$r.squared
r2_left_knee_fl

#   Knee abduction

#   Right
#     ICC
vicon  <- as.numeric(matrix(data[-1,13]))
kinect <- as.numeric(matrix(data[-1,14]))
icc_right_knee_ab <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_right_knee_ab
#     R2
right_knee_ab.lm <- lm(vicon ~ kinect)
r2_right_knee_ab <- summary(right_knee_ab.lm)$r.squared
r2_right_knee_ab

#   Left
#     ICC
vicon  <- as.numeric(matrix(data[-1,15]))
kinect <- as.numeric(matrix(data[-1,16]))
icc_left_knee_ab <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_left_knee_ab
#     R2
left_knee_ab.lm <- lm(vicon ~ kinect)
r2_left_knee_ab <- summary(left_knee_ab.lm)$r.squared
r2_left_knee_ab

#   KASR
#     ICC
vicon  <- as.numeric(matrix(data[-1,17]))
kinect <- as.numeric(matrix(data[-1,18]))
vicon <- vicon[!is.na(vicon)]
kinect <- kinect[!is.na(kinect)]
icc_kasr <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_kasr
#     R2
kasr.lm <- lm(vicon ~ kinect)
r2_kasr <- summary(kasr.lm)$r.squared
r2_kasr

#   Performance variable
#     ICC
vicon  <- as.numeric(matrix(data[-1,19]))
kinect <- as.numeric(matrix(data[-1,20]))
icc_perf <- icc(cbind(vicon, kinect), "twoway", "agreement", "single", conf.level = 0.95)
icc_perf
#     R2
perf.lm <- lm(vicon ~ kinect)
r2_perf <- summary(perf.lm)$r.squared
r2_perf



# Vicon vs. Curv (ICC on critical values)

#   Set working directory
setwd("C:/Users/bdour/Documents/Work/Toronto/Sunnybrook/ACL Injury Screening/Data/Processed/Biomechanical variables/Comparison/Vicon vs. Curv")

#   Data import
input_file <- "V-vs-C_critical_values.csv" 
data <- read.csv(file = input_file, header = TRUE)
length_data = dim(data)[1]-2
data <- data[1:length_data,]

#   Library import
library(irr)

#   Knee abduction

#   Right
#     ICC
vicon  <- as.numeric(matrix(data[-1,3]))
curv <- as.numeric(matrix(data[-1,4]))
icc_right_knee_ab <- icc(cbind(vicon, curv), "twoway", "agreement", "single", conf.level = 0.95)
icc_right_knee_ab
#     R2
right_knee_ab.lm <- lm(vicon ~ curv)
r2_right_knee_ab <- summary(right_knee_ab.lm)$r.squared
r2_right_knee_ab

#   Left
#     ICC
vicon  <- as.numeric(matrix(data[-1,5]))
curv <- as.numeric(matrix(data[-1,6]))
icc_left_knee_ab <- icc(cbind(vicon, curv), "twoway", "agreement", "single", conf.level = 0.95)
icc_left_knee_ab
#     R2
left_knee_ab.lm <- lm(vicon ~ curv)
r2_left_knee_ab <- summary(left_knee_ab.lm)$r.squared
r2_left_knee_ab

#   KASR
#     ICC
vicon  <- as.numeric(matrix(data[-1,7]))
curv <- as.numeric(matrix(data[-1,8]))
vicon <- vicon[!is.na(vicon)]
curv <- curv[!is.na(curv)]
icc_kasr <- icc(cbind(vicon, curv), "twoway", "agreement", "single", conf.level = 0.95)
icc_kasr
#     R2
kasr.lm <- lm(vicon ~ curv)
r2_kasr <- summary(kasr.lm)$r.squared
r2_kasr