#' Specify random correlations
#' @export
random <- function(lower, upper) {
    force(lower)
    force(upper)
    # lower < upper ifnot error
    if (lower >= upper) {
        cli::cli_abort(
            'lower must be less than upper in {.cls random}'
        )
    }
    `_fixed` <- FALSE
    `_call`  <- call('random', lower, upper)
    `_mean`  <- mean(lower, upper)
    structure(
        function(n) runif(n, lower, upper),
        class = c('mp_corr_func', 'function')
    )
}

#' Specify fixed correlations
#' @export
fixed <- function(value) {
    force(value)
    `_fixed` <- TRUE
    `_call`  <- call('fixed', value)
    `_mean`  <- value
    structure(
        function(n) rep(value, n),
        class = c('mp_corr_func', 'function')
    )
}

#' Internal function to check if a specific corr is fixed
#'
#' @noRd
is_fixed <- function(x) {
    if (!is.function(x)) return(TRUE)
    isTRUE(environment(x)$`_fixed`)
}

#' Convert `mp_corr_func` to a character
#' @noRd
as.character.mp_corr_func <- function(x, ...) {
    if (is_fixed(x)) as.character(environment(x)$`_value`)
    else deparse(environment(x)$`_call`)
}

#' Validate `mp_corr_func`
#' @noRd
is.corr_func <- function(x) {
    inherits(x, "mp_corr_func")
}

#' Obtain mean of `mp_corr_func`
#' @noRd
mean.mp_corr_func <- function(x, ...) {
    environment(x)$`_mean`
}

#' Internal function to check if all correlations are fixed
#'
#' @noRd
is_fixed_cor <- function(x) {
    is_fixed(x$within_cor) &
        is_fixed(x$between_cor) &
        is_fixed(x$randeff_cor)
}

#' Internal function for mp_correlations object
#'
#' @noRd
mp_correlations <- function(within_cor, between_cor, randeff_cor) {
    # Create environment
    e <- list2env(
        list(
            within_cor = within_cor,
            between_cor = between_cor,
            randeff_cor = randeff_cor
        ), parent = emptyenv() )

    # Return object
    structure(e, class = 'mp_correlations')
}

#' Set up correlation defaults
#'
#' @export
correlations <- function(
        within_cor  = random(0.1, 0.3),
        between_cor = random(0.1, 0.3),
        randeff_cor = random(0.1, 0.3)) {

    # Check if single number is specified and set fixed
    if (is.number(within_cor))  within_cor  <- fixed(within_cor)
    if (is.number(between_cor)) between_cor <- fixed(between_cor)
    if (is.number(randeff_cor)) randeff_cor <- fixed(randeff_cor)

    # Validate that they are corr funcs
    if (!is.corr_func(within_cor)) cli::cli_abort(
        '{.cli within_cor} must be a single number or created via `fixed` or `random` functions'
    )
    if (!is.corr_func(between_cor))  cli::cli_abort(
        '{.cli between_cor} must be a single number or created via `fixed` or `random` functions'
    )
    if (!is.corr_func(randeff_cor))  cli::cli_abort(
        '{.cli randeff_cor} must be a single number or created via `fixed` or `random` functions'
    )
    # Construct correlations object
    mp_correlations(within_cor, between_cor, randeff_cor)
}

#' Internal function to create default correlations
#'
#' @noRd
default_correlations <- function() {
    e <- correlations()
    attr(e, 'default') <- TRUE
    return(e)
}

#' Adds correlations to `mp_base` class
#'
#' @noRd
`add.mp_correlations` <- function(x, y) {
    # Add as predictor if model
    if (is.model(x)) {
        # check if predictor already there
        if (is.null(attr(x$corrs, 'default')))  cli::cli_alert_warning(
            'Correlations have already been specified. Overwritting previous ones.'
        )
        x$corrs <- y
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}
