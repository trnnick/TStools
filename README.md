TStools 
=======

Collection of R functions for time series analysis

To install use:

> if (!require("devtools")){install.packages("devtools")}

> devtools::install_github("trnnick/TStools")

The package now depends on Rcpp and RcppArmadillo, which will be installed automatically.

However Mac OS users may need to install gfortran libraries in order to use Rcpp. Follow the link for the instructions: http://www.thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

Note
-------
Functions es() and sim.ets() are now maintained in "smooth" package: https://github.com/config-i1/smooth.
The functions will be removed from Tstools by the end of May 2016.
