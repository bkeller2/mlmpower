% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/actions.R
\name{product}
\alias{product}
\title{Create a Product Term in a Model}
\usage{
product(name1, name2, weight = 1)
}
\arguments{
\item{name1}{a character string that references the first variable's name}

\item{name2}{a character string that references the second variable's name}

\item{weight}{a single numeric value specifying the variable's contribution to the variance explained metric.
Weights are normalized across all variables of the same level.}
}
\value{
A \code{\link{mp_action}} that can be added to a \code{\link{mp_model}}.
}
\description{
Creates a product term between two variables that can be added to a \code{\link{mp_model}}.
}
\details{
Currently the product term is only limited to cross-level
interactions between a level-1 centered within cluster variable (\code{icc = 0})
and level-2 variable.
}
\examples{
# Create Model
model <- (
    outcome('Y')
    + within_predictor('X', icc = 0.0)
    + between_predictor('Z')
)
# Add random slope to the model
model + product('X', 'Z')
}
