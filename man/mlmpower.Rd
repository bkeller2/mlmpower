% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{mlmpower}
\alias{mlmpower}
\alias{mp_model}
\alias{model}
\alias{Modeling}
\alias{modeling}
\alias{`+.mp_base`}
\alias{mp_action}
\alias{mp_base}
\title{\code{mlmpower} Modeling Framework}
\description{
\code{mlmpower} constructs models by adding different features of the model using the plus sign  (\code{+}).

Every model requires an \code{\link{outcome}} and an ICC specified in \code{\link{effect_size}} to be valid.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{model <- outcome('y') + effect_size(icc = 0.1)
}\if{html}{\out{</div>}}

Once a model is constructed, we can add additional features to build the model out more.
For example, we may want to include a level-1 predictor that is centered within cluster.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{model <- model + within_predictor('x', icc = 0.0)
}\if{html}{\out{</div>}}

The additions can be chained together to produce the entire model object.
For example, the previous two code blocks can be combined into one.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{model <- (
    outcome('y')
    + effect_size(icc = 0.1)
    + within_predictor('x', icc = 0.0)
)
}\if{html}{\out{</div>}}

Finally, we can also wrap multiple variables into a list and add that.
This feature can be useful when programmatically generating a model.

\if{html}{\out{<div class="sourceCode r">}}\preformatted{model <- (
    outcome('y')
    + effect_size(icc = 0.1)
    + lapply(1:10, \\(i) within_predictor(paste0('x', i), icc = 0.0))
)
}\if{html}{\out{</div>}}

For more detailed information see the help vignette by running the following:

\if{html}{\out{<div class="sourceCode">}}\preformatted{vignette('mlmpower')
}\if{html}{\out{</div>}}
}
\seealso{
\link{Variables} \code{\link[=effect_size]{effect_size()}}
\code{\link[=correlations]{correlations()}} \code{\link[=random_slope]{random_slope()}}
\code{\link[=product]{product()}}
}
