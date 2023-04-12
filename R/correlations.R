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
    `_fixed` <- F
    `_call`  <- call('random', lower, upper)
    function(n) {
        runif(n, lower, upper)
    }
}

#' Specify fixed correlations
#' @export
fixed <- function(value) {
    force(value)
    `_fixed` <- T
    `_call`  <- call('fixed', value)
    function(n) {
        rep(value, n)
    }
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

#' Internal function to check if a specific corr is fixed
#'
#' @noRd
is_fixed_cor <- function(x) {
    if (!is.function(x)) return(TRUE)
    isTRUE(environment(x)$`_fixed`)
}

#' Internal function to check if all correlations are fixed
#'
#' @noRd
is_fixed <- function(x) {
    is_fixed_cor(x$within_cor) &
    is_fixed_cor(x$between_cor) &
    is_fixed_cor(x$randeff_cor)
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
