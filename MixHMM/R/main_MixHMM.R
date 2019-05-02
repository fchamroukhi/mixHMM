###########################################################################
# Clustering and segmentation of time series (including with regime
# changes) by mixture of gaussian Hidden Markov Models (MixFHMMs) and the EM algorithm; i.e functional data
# clustering and segmentation
#
#
#
#
# by Faicel Chamroukhi, 2009
#
## Please cite the following references for this code
#
# @InProceedings{Chamroukhi-IJCNN-2011,
#   author = {F. Chamroukhi and A. Sam\'e  and P. Aknin and G. Govaert},
#   title = {Model-based clustering with Hidden Markov Model regression for time series with regime changes},
#   Booktitle = {Proceedings of the International Joint Conference on Neural Networks (IJCNN), IEEE},
#   Pages = {2814--2821},
#   Adress = {San Jose, California, USA},
#   year = {2011},
#   month = {Jul-Aug},
#   url = {https://chamroukhi.com/papers/Chamroukhi-ijcnn-2011.pdf}
# }
#
# @PhdThesis{Chamroukhi_PhD_2010,
# author = {Chamroukhi, F.},
# title = {Hidden process regression for curve modeling, classification and tracking},
# school = {Universit\'e de Technologie de Compi\`egne},
# month = {13 december},
# year = {2010},
# type = {Ph.D. Thesis},
# url ={https://chamroukhi.com/papers/FChamroukhi-Thesis.pdf}
# }
#
# @article{Chamroukhi-FDA-2018,
# 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
#	  Author = {Faicel Chamroukhi and Hien D. Nguyen},
#	  Note = {DOI: 10.1002/widm.1298.},
#	  Volume = {},
# 	Title = {Model-Based Clustering and Classification of Functional Data},
#	  Year = {2018},
#	  Month = {Dec},
#	  url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
#	}
###########################################################################

rm(list = ls())
source("R/FData.R")
source("R/enums.R")
source("R/ModelMixHMM.R")
source("R/ModelLearner.R")


# Load data
load("data/simulatedTimeSeries.RData")
fData <- FData$new()
fData$setData(X, Y)

# Model specification
K <- 3 # Number of clusters
R <- 3 # Number of regimes (HMM states)
variance_type = variance_types$hetereskedastic

modelMixHMM <- ModelMixHMM(fData, K, R, variance_type)

ordered_states <- TRUE
n_tries <- 1
max_iter <- 1000
init_kmeans <- TRUE
threshold <- 1e-6
verbose <- TRUE

solution <- EM(modelMixHMM, ordered_states, n_tries, max_iter, init_kmeans, threshold, verbose)

solution$plot()
