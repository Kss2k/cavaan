% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{TPB_1SO}
\alias{TPB_1SO}
\title{TPB_1SO}
\description{
A simulated dataset based on the Theory of Planned Behaviour,
where INT is a higher order construct of ATT, SN, and PBC.
}
\examples{
tpb <- '
  # First order constructs
  ATT =~ att1 + att2 + att3
  SN  =~ sn1 + sn2 + sn3
  PBC =~ pbc1 + pbc2 + pbc3
  BEH =~ b1 + b2

  # Higher order constructs
  INT =~ ATT + PBC + SN
  # Structural model
  BEH ~ PBC + INT
'

}
