## eDNAoccupancy


eDNAoccupancy is an R package for fitting multiscale occupancy models. Although designed for the analysis of surveys of environmental DNA (eDNA), this package can be used for any occupancy survey that includes three nested levels of sampling:

1. locations within a study area,
2. samples of each location,
3. replicate observations of each sample.

Instructions and examples of using eDNAoccupancy are available in the following publication:

* [Dorazio, RM and Erickson, RA. 2018.  eDNAoccupancy:  An R package for multiscale occupancy modelling of environmental DNA data.  Molecular Ecology Resources 18: 368-380.](https://doi.org/10.1111/1755-0998.12735)



## Instructions for installing eDNAoccupancy

1. Start R

2. Install R packages `mvtnorm`, `pROC`, `mcmcse`, and `devtools` from the CRAN repository.

3. At the R command line, type

        devtools::install_github("RobertDorazio/eDNAoccupancy")

To use the package, at the R command line, type

    library(eDNAoccupancy)

An introductory user's guide is available by typing

    vignette("eDNAoccIntro")



## Contact information

Robert Dorazio (rdorazio@sfsu.edu) developed this R package with assistance from Richard Erickson (rerickson@usgs.gov).



## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the [official USGS copyright policy](https://www2.usgs.gov/visual-id/credit_usgs.html#copyright/).


This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review.  No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.

This software is provided "AS IS".

