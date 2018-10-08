#'  Fungal Pathogen Detection Data
#'
#'  A data set containing detections and non-detections of a fungal pathogen's DNA measured at each of 20 ponds located in Coconino National Forest, Arizona, USA.  At each site, water samples were collected to extract the DNA of \emph{Batrachochytrium dendrobatidis} (Bd), a fungal pathogen of various species of amphibians.
#'
#'
#' @format  A data frame consisting of 80 rows and the following 4 columns:
#' \describe{
#' \item{site:}{ location sampled for fungal pathogen }
#' \item{sample:}{ sample number (integer-valued) }
#' \item{pcr1:}{ binary indicator of whether Bd DNA was detected (1) or not detected (0) in 1st PCR replicate of sample }
#' \item{pcr2:}{ binary indicator of whether Bd DNA was detected (1) or not detected (0) in 2nd PCR replicate of sample }
#' }
#'
#' 
#' @docType data
#'
#' 
#' @usage data(fungusDetectionData)
#'
#'
#'
#'
#' @source  Schmidt BR, Kery M, Ursenbacher S, Hyman OJ, and Collins JP (2013) Site occupancy models in the analysis of environmental DNA presence/absence surveys: a case study of an emerging amphibian pathogen.  Methods in Ecology and Evolution 4: 646-653.
#'
#'
#' @references
#' Schmidt BR, Kery M, Ursenbacher S, Hyman OJ, and Collins JP (2013) Site occupancy models in the analysis of environmental DNA presence/absence surveys: a case study of an emerging amphibian pathogen.  Methods in Ecology and Evolution 4: 646-653.
#' 
#' \url{http://dx.doi.org/10.1111/2041-210X.12052}
#'
#' Hyman, OJ, Collins JP (2012) Evaluation of a filtration-based method for detecting  \emph{Batrachochytrium dendrobatidis} in natural bodies of water.  Diseases of Aquatic Organisms 97: 185-195.
#' 
#' \url{ https://doi.org/10.3354/dao02423}
"fungusDetectionData"
