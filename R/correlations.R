#' Specify random correlations
#' @export
random <- function(lower, upper) {
    force(lower)
    force(upper)
    # lower < upper ifnot error
    if (lower >= upper) {
        throw_error('{.arg lower} must be less than upper in in {.fun random}')
    }
    if (lower < -1 | lower > 1) {
        throw_error('{.arg lower} must be between -1 and 1 in {.fun random}')
    }
    if (upper < -1 | upper > 1) {
        throw_error('{.arg upper} must be between -1 and 1 in {.fun random}')
    }

    `_fixed` <- FALSE
    `_call`  <- call('random', lower, upper)
    `_mean`  <- mean(c(lower, upper))
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

    if (value < -1 | value > 1) {
        throw_error('{.arg value} must be between -1 and 1 in {.fun fixed}')
    }
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
    if (is_fixed(x)) as.character(environment(x)$`value`)
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
    is_fixed(x$within) &
        is_fixed(x$between) &
        is_fixed(x$randeff)
}

#' Check if types are valid for a correlation
#'
#' @noRd
valid_corr_type <- function(x) {
    switch(
        x,
        within = TRUE,
        between = TRUE,
        randeff = TRUE,
        FALSE # Default
    )
}


#' Internal function for mp_correlations object
#'
#' @noRd
new_correlation <- function(type, value) {
    if (missing(type)) {
        throw_error('Must provide effect size type for {.cls mp_correlation}')
    }
    if (!valid_corr_type(type)) {
        throw_error('Invalid type for {.cls mp_correlation}')
    }
    if (!is.corr_func(value)) throw_error(
        '{.arg {type}} must be a single number or created via {.fun fixed} or {.fun random}'
    )
    # Create new effect sizes
    structure(
        list2env(
            list(type = type, value = value),
            parent = emptyenv()
        ),
        class = c('mp_corr', 'mp_base')
    )
}

#' Set up correlation defaults
#'
#' @export
correlations <- function(
        within,
        between,
        randeff) {

    # Create output
    o <- mp_list()

    # Check if single number is specified and set fixed
    if (!missing(within)) {
        if (is.number(within))  within  <- fixed(within)
        o |> add(new_correlation('within', within)) -> o
    }
    if (!missing(between)) {
        if (is.number(between)) between <- fixed(between)
        o |> add(new_correlation('between', between)) -> o
    }
    if (!missing(randeff)) {
        if (is.number(randeff)) randeff <- fixed(randeff)
        o |> add(new_correlation('randeff', randeff)) -> o
    }
    # Return list
    return(o)
}

#' Internal function to create default correlations
#'
#' @noRd
default_correlations <- function() {
    structure(
        list2env(
            list(
                within  = structure(random(0.1, 0.3), default = TRUE),
                between = structure(random(0.1, 0.3), default = TRUE),
                randeff = structure(random(0.1, 0.3), default = TRUE)
            ),
            parent = emptyenv()
        ),
        class = c('mp_correlations', 'mp_base')
    )
}

#' Adds correlations to `mp_base` class
#'
#' @noRd
add.mp_corr <- function(x, y) {
    # Add as correlation if model
    if (is.model(x)) {
        x$corrs[[y$type]] <- y$value
        return(x)
    }

    # Otherwise construct list
    x |> add(mp_list(y))
}


#' Prints `mp_correlations` class
#'
#' @export
print.mp_correlations <- function(x, ...) {
    cli::cli_ul(
        c(
            'WITHIN = {as.character(x$within)}',
            'BETWEEN = {as.character(x$between)}',
            'RANDOM EFFECT = {as.character(x$randeff)}'
        )
    )
    invisible(x)
}
