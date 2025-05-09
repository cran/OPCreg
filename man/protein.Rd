\name{protein}
\alias{protein}
\docType{data}
\title{Protein Secondary Structure Data}
\description{
This dataset contains protein sequences and their corresponding secondary structures, including beta-sheets (E), helices (H), and coils (_).
}
\usage{
protein
}
\format{
A data frame with multiple rows and columns representing protein sequences and their secondary structures.
\itemize{
\item \code{Sequence}: Amino acid sequence (using 3-letter codes).
\item \code{Structure}: Secondary structure of the protein (E for beta-sheet, H for helix, _ for coil).
\item \code{Parameters}: Additional parameters for neural networks (to be ignored).
\item \code{Biophysical_Constants}: Biophysical constants (to be ignored).
}
}
\details{
The dataset is used for predicting protein secondary structures from amino acid sequences. The first few numbers in each sequence are parameters for neural networks and should be ignored. The '<' symbol is used as a spacer between proteins and to mark the beginning and end of sequences.
}
\note{
The biophysical constants included in the dataset were found to be unhelpful and are generally ignored in analysis.
}
\source{
Vince G. Sigillito, Applied Physics Laboratory, Johns Hopkins University.
}
\examples{
# Load the dataset
data(protein)

# Print the first few rows of the dataset
print(head(protein))
}
\keyword{datasets}
