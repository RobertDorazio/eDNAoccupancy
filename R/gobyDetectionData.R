#'  Tidewater Goby Detection Data
#'
#'  A data set containing detections and non-detections of tidewater goby DNA in water samples collected at each of 39 sites located along 400 km of coastline in California and Oregon, USA.  At each site, water samples were collected to extract the DNA of the tidewater goby  (\emph{Eucyclogbobius newberryi}), a federally endangered fish species living in estuaries, lagoons, and sloughs.
#'
#'
#' @format  A data frame consisting of 356 rows and the following 8 columns:
#' \describe{
#' \item{site:}{ location sampled for tidewater goby DNA }
#' \item{sample:}{ sample number (integer-valued) }
#' \item{pcr1:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 1st PCR replicate of sample }
#' \item{pcr2:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 2nd PCR replicate of sample }
#' \item{pcr3:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 3rd PCR replicate of sample }
#' \item{pcr4:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 4th PCR replicate of sample }
#' \item{pcr5:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 5th PCR replicate of sample }
#' \item{pcr6:}{ binary indicator of whether tidewater goby DNA was detected (1) or not detected (0) in 6th PCR replicate of sample }
#' }
#'
#' 
#' @docType data
#'
#' 
#' @usage data(gobyDetectionData)
#'
#'
#'
#'
#' @source  Schmelzle MC, Kinziger AP (2015) Data from: Using occupancy modeling to compare environmental DNA to traditional field methods for regional-scale monitoring of an endangered aquatic species.
#' \describe{
#' \item{Dryad Digital Repository:}{ \url{http://dx.doi.org/10.5061/dryad.6rs23} }
#' \item{File name:}{ "Hierarchical Occupancy Model Input Data.xls" }
#' }
#'
#'
#' @references  Schmelzle MC, Kinziger AP (2016) Using occupancy modelling to compare environmental DNA to traditional field methods for regional-scale monitoring of an endangered aquatic species.  Molecular Ecology Resources 16: 895–908.
#' 
#' \url{http://dx.doi.org/10.1111/1755-0998.12501}
"gobyDetectionData"
