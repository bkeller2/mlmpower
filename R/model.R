#' @rdname mlmpower
#' @title `mlmpower` Modeling Framework
#' @name mlmpower
#' @aliases mp_model model Modeling modeling `+.mp_base` mp_action mp_base
#' @seealso
#' [mlmpower::Variables] [mlmpower::effect_size()]
#' [mlmpower::correlations()] [mlmpower::random_slope()]
#' [mlmpower::product()]
#' @description
#' `mlmpower` constructs models by adding different features of the model using the plus sign  (`+`).
#'
#' Every model requires an [`mlmpower::outcome`] and an ICC specified in [`mlmpower::effect_size`] to be valid.
#' ```{r}
#' model <- outcome('y') + effect_size(icc = 0.1)
#' ```
#' Once a model is constructed, we can add additional features to build the model out more.
#' For example, we may want to include a level-1 predictor that is centered within cluster.
#' ```{r}
#' model <- model + within_predictor('x', icc = 0.0)
#' ```
#' The additions can be chained together to produce the entire model object.
#' For example, the previous two code blocks can be combined into one.
#' ```{r}
#' model <- (
#'     outcome('y')
#'     + effect_size(icc = 0.1)
#'     + within_predictor('x', icc = 0.0)
#' )
#' ```
#' Finally, we can also wrap multiple variables into a list and add that.
#' This feature can be useful when programmatically generating a model.
#' ```{r}
#' model <- (
#'     outcome('y')
#'     + effect_size(icc = 0.1)
#'     + lapply(1:10, \(i) within_predictor(paste0('x', i), icc = 0.0))
#' )
#' ```
#'
#' For more detailed information see the help vignette by running the following:
#' ```
#' vignette(package = 'mlmpower')
#' ```
NULL

#' Generate Multilevel base class
#' @noRd
make_model <- function(o) {
    # Create environment
    e <- list2env(
        list(
            outcome = o,
            effect_size = default_effect_size(),
            predictors = list(),
            actions = list(),
            corrs = default_correlations()
        )
    )

    # Return object
    structure(e, class = c('mp_model', 'mp_base'))
}

#' Check if `mp_model`
#' @noRd
is.model <- function(x) {
    inherits(x, 'mp_model')
}

#' Check if a Model is Properly Specified
#' @description
#' This function is used to validate if a [`mlmpower::mp_model`] is correct.
#' If the model is incorrect an appopriate error message describing while will be supplied
#' @param x a [`mlmpower::mp_model`]
#' @returns Invisibly returns the original model.
#' @examples
#' # Create Model
#' model <- outcome('Y') + within_predictor('X')
#' # Throws error
#' tryCatch(
#'     is_valid(model),
#'     error = print
#' )
#' # Succeeds
#' is_valid(model + effect_size(icc = 0.1))
#' @export
is_valid <- function(x) {

    # Check that it is a model
    if (!is.model(x)) throw_error(c(
        'The {.cls model} object is not properly constructed.',
        'x' = 'The object is not of class {.cls mp_model}'
    ))

    # Get effect size
    es <- x$effect_size

    # Check if was specified and any null icc in predictors
    if (is.null(es$icc))  {
        sel <- vapply(x$predictors, \(.) is.null(.$icc), logical(1L))
        if (sum(sel) != 0) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'The ICC was never specified for the model.'
        ))
    }

    # Check if random slopes exist with non 0 effect size
    if (es$random_slope > 0) {
        count <- sum(vapply(x$actions, \(.) .$type == "random_slope", logical(1L)))
        if (count == 0) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'An effect size for random slopes was specified with no random slopes.'
        ))
        # Check if correlation is non zero for more than 1 random slope
        else if (count > 1) {
            if (!is_fixed(x$corrs$randeff) | mean(x$corrs$randeff) != 0) {

                # Check if default and swap to zero
                if (isTRUE(attr(x$corrs$randeff, 'default'))) {
                    x$corrs$randeff <- fixed(0)
                } else {
                    throw_error(c(
                        'The {.cls model} object is not properly constructed.',
                        'x' = 'Multiple random slopes exist with non-zero {.arg randeff} correlation.',
                        'i' = 'The {.arg randeff} must be set to 0 with multiple random slopes.'
                    ))
                }
            }
        }
    }

    # Check if product exist with non 0 effect size
    if (es$product > 0) {
        sel <- vapply(x$actions, \(.) .$type == "product", logical(1L))
        if (sum(sel) == 0) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'An effect size for products was specified with no products.'
        ))
    }

    # Check if level-1 exist with non 0 effect size
    if (es$within > 0) {
        sel <- vapply(x$predictors, levels, numeric(1L)) == 1
        if (sum(sel) == 0) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'An effect size for within was specified with no within predictors.'
        ))
    }

    # Check if level-1 exist with non 0 effect size
    if (es$between > 0) {
        sel <- vapply(x$predictors, levels, numeric(1L)) == 2
        if (sum(sel) == 0) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'An effect size for between was specified with no between predictors.'
        ))
    }

    # Check if effect sizes don't sum to 1 or greater
    total <- sum(es$within + es$between + es$product + es$random_slope)
    if (total >= 1) {
        throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'The total ({total}) effect size is greater or equal to 1'
        ))
    }

    # Check if r2 between is too large based on icc
    for (icc in (if (is.null(x$outcome$icc)) es$icc else outcome$icc)) {
        if (es$between >= icc) throw_error(c(
            'The {.cls model} object is not properly constructed.',
            'x' = 'The between effect size is greater or equal to the ICC'
        ))
    }

    # Return object if valid
    invisible(x)
}


#' Internal function for cloning a model
#' @noRd
clone.mp_model <- function(x) {
    # Create environment
    e <- list2env(
        list(
            outcome = x$outcome,
            effect_size = x$effect_size,
            predictors = x$predictors,
            actions = x$actions,
            corrs = x$corrs
        )
    )

    # Return new model
    structure(e, class = c('mp_model', 'mp_base'))
}

#' Subset a [`mlmpower::mp_model`] by Global ICC
#' @description
#' Subsets a [`mlmpower::mp_model`] with multiple ICC values specified in [`mlmpower::effect_size`]
#' into a model with only the single ICC value.
#' @param x a [`mlmpower::mp_model`] object
#' @param icc a single numeric value to subset out of `x`
#' @param ... other arguments not used by this method.
#' @returns A new [`mlmpower::mp_model`] with only the subset ICC
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = cross_sectional)
#' )
#' # Obtain Model with only 0.15 ICC
#' model |> subset(icc = 0.15)
#' @export
subset.mp_model <- function(x, icc, ...) {

    # Validate model first
    is_valid(x)

    if (!is.number(icc)) {
        throw_error( 'icc needs to be a single number')
    }
    else if ( !(icc %in% x$effect_size$icc)) {
        throw_error( 'icc is not in {.cls mp_model}')
    }
    new_x <- clone(x)
    new_x$effect_size <- clone(x$effect_size)
    new_x$effect_size$icc <- icc
    return(new_x)
}

#' Prints a [`mlmpower::mp_model`]
#' @description
#' Prints a [`mlmpower::mp_variable`] in a human readable format.
#' @param x a [`mlmpower::mp_model`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @examples
#' print(
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = cross_sectional)
#' )
#' @export
print.mp_model <- function(x, ...) {

    # Validate model first
    is_valid(x)

    # Print specificaiton
    cli::cli_h2('{.cli mlmpower} model specification')
    cli::cli_h3('{.cli outcome}')
    cli::cli_text('')
    print(x$outcome)
    cli::cli_h3('{.cli predictors}')
    cli::cli_text('')
    print(
        do.call(
            'rbind',
            lapply(x$predictors, variable_to_row)
        )[, c('level', 'icc', 'weight', 'type')],
        right = FALSE
    )

    if (length(x$actions) > 0) {
        cli::cli_h3('{.cli model specifications}')
        cli::cli_text('')
        print(
            do.call(
                'rbind',
                lapply(x$actions, action_to_row)
            ),
            right = FALSE,
            row.names = FALSE
        )
    }

    cli::cli_h3('{.cli effect sizes}')
    cli::cli_text('')
    print(x$effect_size)

    cli::cli_h3('{.cli correlations}')
    cli::cli_text('')
    print(x$corrs)
}


#' Obtain the Parameter Summaries for A [`mlmpower::mp_model`]
#' @description
#' Provide the summarized parameter estimates for a [`mlmpower::mp_model`],
#' including the variance explained break downs.
#'
#' @param object a [`mlmpower::mp_model`]
#' @param ... other arguments not used by this method.
#' @returns
#' A [`mlmpower::mp_parameters`] object that contains the population parameters based on the model.
#' If random correlations are used the average correlation is used to compute the parameters.
#' If multiple ICC's are specified then a named [`base::list`] is
#' returned containing the parameter value for each ICC value.
#' @examples
#' summary(
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = cross_sectional)
#' )
#' @export
summary.mp_model <- function(object, ...) {

    # validate model
    is_valid(object)

    # Get icc
    icc <- object$effect_size$icc

    # Check if fixed and warn if not
    if (!is_fixed_cor(object$corrs))  {
        cli::cli_alert_warning('Parameters are solved based on average correlation.')
    }

    # Handle multiple ICCs
    if (length(icc) > 1) {

        # Create names
        names(icc) <- paste0('icc = ', icc)

        # Run simulation
        return(
            lapply(icc, \(x) {
                object |> subset(icc = x) |> make_avg_parameters() |> clean_parameters()
            })
        )
    }
    # Otherwise run once and return
    object |> make_avg_parameters() |> clean_parameters()
}
