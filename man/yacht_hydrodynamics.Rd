\name{yacht_hydrodynamics}
\alias{yacht_hydrodynamics}
\docType{data}
\title{Yacht Hydrodynamics Data}
\description{
This dataset contains the hydrodynamic characteristics of sailing yachts, including design parameters and performance metrics.
}
\usage{
yacht_hydrodynamics
}
\format{
A data frame with 308 rows and 7 columns.
\itemize{
\item \code{Residuary Resistance}: Residuary resistance per unit weight of displacement (performance metric).
\item \code{Longitudinal Position of Center of Buoyancy}: Longitudinal position of the center of buoyancy.
\item \code{Prismatic Coefficient}: Prismatic coefficient.
\item \code{Length-Displacement Ratio}: Length-displacement ratio.
\item \code{Beam-Draft Ratio}: Beam-draft ratio.
\item \code{Length-Beam Ratio}: Length-beam ratio.
\item \code{Froude Number}: Froude number.
}
}
\details{
The dataset contains hydrodynamic data for sailing yachts, with the goal of predicting the residuary resistance from various design parameters.
}
\note{
The dataset is commonly used for regression analysis and machine learning tasks to model the relationship between design parameters and performance metrics.
}
\source{
UCI Machine Learning Repository
}
\examples{
# Load the dataset
data(yacht_hydrodynamics)

# Print the first few rows of the dataset
print(head(yacht_hydrodynamics))

}
\keyword{datasets}
