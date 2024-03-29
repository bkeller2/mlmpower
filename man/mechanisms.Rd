% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mechanisms.R
\name{mechanisms}
\alias{mechanisms}
\alias{mechanism}
\alias{MCAR}
\alias{MAR}
\title{Helper functions for producing Missing Data Mechanisms}
\usage{
# Generate MCAR data on outcome
MCAR(mis.rate)

# Generate MAR data on outcome due to `cause`
MAR(mis.rate, cause, r2, lower = TRUE)
}
\arguments{
\item{mis.rate}{A proportion for the missing data rate at population level}

\item{cause}{A character for a variable name that is the cause of missingness}

\item{r2}{A proportion of variance explained by the cause in the missing data indicator's latent propensity}

\item{lower}{A logical for the lower or upper tail being more likely to be missing}
}
\description{
Functions to generate data that always follows a specific mechanism
in accordance to a single-level model.
}
\examples{
# Create Model
model <- (
    outcome('Y')
    + within_predictor('X')
    + effect_size(icc = 0.1)
)

# Induce MCAR data on outcome
set.seed(19723)
model |> power_analysis(50, 5, 50, mechanism = MCAR(0.25)) -> powersim_mcar

# Induce MAR data on outcome
set.seed(19723)
model |> power_analysis(
   50, 5, 50,
   mechanism = MAR(0.25, 'X', 0.6)
) -> powersim_mar
}
\seealso{
\link{power_analysis}
}
