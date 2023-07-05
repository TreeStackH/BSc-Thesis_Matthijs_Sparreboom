source('C:\\Users\\matth\\anaconda3\\envs\\NetworkFusion.R')
source('C:\\Users\\matth\\anaconda3\\envs\\HelperFunctions.R')

library(Rcpp)
library(parallel)
library(Matrix)
library(ggplot2)
library(dplyr)
library(viridis)
setwd("C:\\Users\\matth\\anaconda3\\envs")
# system("R CMD SHLIB projsplx_R.c")
dyn.load("C:\\Users\\matth\\anaconda3\\envs\\projsplx_R.dll")# connect anaconda environment
library(reticulate)
use_virtualenv("C:\\Users\\matth\\anaconda3\\envs\\MDICC_1")

# the path of python
Sys.setenv(RETICULATE_PYTHON="C:\\Users\\matth\\anaconda3\\envs\\MDICC_1\\Scripts\\python.exe")
use_python("C:\\Users\\matth\\anaconda3\\envs\\MDICC_1\\Scripts\\python.exe", required = TRUE)
py_config()
py_available()
source_python("LocalAffinityMatrix.py")
source_python("score.py")
source_python("label.py")# read data
source_python("feature_Selection.py")
source_python("feature_Selection - SelectKBest.py")

#Kbest
division <- 1000
brca_K_Best <- getKBestScores("brca", division)
lihc_K_Best <- getKBestScores("lihc", division)

#Boruta and Lasso
brca_boruta <- basic_MDICC_LB("brca", "boruta")
brca_lasso <- basic_MDICC_LB("brca", "lasso")

lihc_boruta <- basic_MDICC_LB("lihc", "boruta")
lihc_lasso <- basic_MDICC_LB("lihc", "lasso")

#RFECV
brca_RFECV <- basic_MDICC_RFECV("brca")
lihc_RFECV <- basic_MDICC_RFECV("lihc")


#These function are used to plot the individual results. The current implementation with ggplot does not allow multiple graphs.
#To see different results, the plots have to be runl 1 by 1

#Plot Kbest
output_Result(brca_K_Best, lihc_K_Best, 2, division)

#Plot RFECV
plot_Bars_RFECV(brca_RFECV, lihc_RFECV, "ARI", lihcARI)
plot_Bars_RFECV(brca_RFECV, lihc_RFECV, "NMI", lihcNMI)

#Plot Lasso and Boruta
plot_Bars_LB(brca_boruta,  lihc_boruta,  "ARI", "Boruta")
plot_Bars_LB(brca_boruta,  lihc_boruta,  "NMI", "Boruta")
plot_Bars_LB(brca_lasso,  lihc_lasso,  "ARI", "Lasso")
plot_Bars_LB(brca_lasso,  lihc_lasso,  "NMI", "Lasso")

#Values for ARI and NMI of the original MDICC model
lihcARI <-0.85742188
lihcNMI <- 0.75878906

#These lines are used to plot all the values. Due to time constraint, the values are hard-coded.
names <- c("BRCA", "BRCA RFECV", "BRCA boruta", "BRCA lasso", "BRCA Kbest", "LIHC", "LIHC RFECV", "LIHC boruta", "LIHC lasso", "LIHC Kbest")
values_ARI <- c(1, 0.9692383, 0.8413086, 0.4826660, 1, lihcARI, 0.8574219, 0.7822266, -0.08685303, 0.8574219)
values_NMI <- c(1, 0.9272461, 0.7368164, 0.2861328, 1, lihcNMI, 0.7587891, 0.6713867, 0.03948975, 0.7587891)
plot_Bars_ALL(names, values_ARI, "ARI")
plot_Bars_ALL(names, values_NMI, "NMI")


