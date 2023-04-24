#' Check if types are valid for an effect
#' @noRd
valid_effect_type <- function(x) {
    switch(
        x,
        icc = TRUE,
        within = TRUE,
        between = TRUE,
        random_slope = TRUE,
        product = TRUE,
        FALSE # Default
    )
}

#' Generate Base effect size for mlmpower
#' @noRd
new_effsize <- function(type, value) {
    if (missing(type)) {
        throw_error('Must provide effect size type for {.cls mp_effect}')
    }
    if (!valid_effect_type(type)) {
        throw_error('Invalid type for {.cls mp_effect}')
    }
    if (type != 'icc' & !is.number(value)) {
        throw_error(paste0('Effect size for {.cls ', type, '} must be a single number.' ))
    }
    # evaluate if it is a function
    if (is.function(value)) {
        value <- value()
    }
    # Check each number
    if (is.numeric(value)) {
        for (i in value) {
            if (is.number(i)) {
                if (i < 0 | i > 1) {
                    throw_error(paste0('Invalid effect size for {.cls ', type, '}' ))
                }
            } else {
                throw_error(paste0('Effect size for {.cls ', type, '} must be a single number' ))
            }
        }
    } else {
        throw_error(paste0('Invalid effect size for {.cls ', type, '}' ))
    }
    # Create new effect sizes
    structure(
        list2env(
            list(type = type, value = value),
            parent = emptyenv()
        ),
        class = c('mp_effect', 'mp_base')
    )
}

#' Specify the Effect Size for the Model
#' @aliases mp_effect mp_effsize
#' @description
#' Creates a list of effect sizes to be added to a [`mlmpower::mp_model`].
#' @param icc a numeric vector of global ICC values for [`mlmpower::mp_variable`] who are left unspecified.
#' Values must be between 0 and 1.
#' @param within a single numeric value that corresponds to the proportion of variance explained by the within variables.
#' @param between a single numeric value that corresponds to the incremental proportion of variance explained by the between variables.
#' @param random_slope a single numeric value that corresponds to the proportion of variance explained by the random slopes.
#' @param product  a single numeric value that corresponds to the proportion of variance explained by the product terms.
#' @returns A list that corresponds to each R2 value.
#' @examples
#' # Set ICCs
#' (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = c(0.1, 0.2),
#'         within = 0.3
#'     )
#' )
#'
#' # With cross-sectional ICC
#' (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = cross_sectional,
#'         within = 0.3
#'     )
#' )
#'
#' # With longitudinal ICC
#' (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = longitudinal,
#'         within = 0.3
#'     )
#' )
#' @export
effect_size <- function(
        icc,
        within,
        between,
        random_slope,
        product) {

    # Create output
    o <- mp_list()

    # Add icc
    if (!missing(icc)) o |> add(new_effsize('icc', icc)) -> o
    # Add within
    if (!missing(within)) o |> add(new_effsize('within', within)) -> o
    # Add between
    if (!missing(between)) o |> add(new_effsize('between', between)) -> o
    # Add random_slope
    if (!missing(random_slope)) o |> add(new_effsize('random_slope', random_slope)) -> o
    # Add product
    if (!missing(product)) o |> add(new_effsize('product', product)) -> o

    # Return list
    return(o)
}

#' Specifies default effect sizes
#' @noRd
default_effect_size <- function() {
    # Create new effect sizes
    structure(
        list2env(
            list(
                icc = NULL,
                within = 0.0,
                between = 0.0,
                random_slope = 0.0,
                product = 0.0
            ),
            parent = emptyenv()
        ),
        class = c('mp_effsize', 'mp_base')
    )
}

#' Internal function for cloning an effect size
#' @noRd
clone.mp_effsize <- function(x) {
    # Create new effect sizes
    structure(
        list2env(
            list(
                icc = x$icc,
                within = x$within,
                between = x$between,
                random_slope = x$random_slope,
                product = x$product
            ),
            parent = emptyenv()
        ),
        class = c('mp_effsize', 'mp_base')
    )
}

#' Longitudinal icc range
#' @rdname effect_size
#' @description Returns suggested ICC's for cross-sectional studies (0.05, 0.15, and 0.25).
#' @export
cross_sectional <- function() {
    c(0.05, 0.15, 0.25)
}

#' cross-sectional icc range
#' @rdname effect_size
#' @description Returns suggested ICC ranges for longitudinal studies (0.40, 0.50, and 0.60).
#' @export
longitudinal <- function() {
    c(0.40, 0.50, 0.60)
}

#' Adds effect size to `mp_base` class
#' @noRd
add.mp_effect <- function(x, y) {
    # Add as effect size if model
    if (is.model(x)) {
        x$effect_size[[y$type]] <- y$value
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}


#' Prints a [`mlmpower::mp_effect`]
#' @description
#' Prints a [`mlmpower::mp_effect`] in a human readable format.
#' @param x a [`mlmpower::mp_effect`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @noRd
print.mp_effect <- function(x, ...) {
    ulid <- cli::cli_ul()
    cli::cli_li('type   = {x$type}')
    cli::cli_li('value  = {x$value}')
    cli::cli_end(ulid)
    invisible(x)
}


#' Prints a [`mlmpower::mp_effsize`]
#' @description
#' Prints a [`mlmpower::mp_effsize`] in a human readable format.
#' @param x a [`mlmpower::mp_effsize`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @examples
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = c(0.1, 0.2),
#'         within = 0.3
#'     )
#' )
#' # Print effect size only
#' print(model$effect_size)
#' @export
print.mp_effsize <- function(x, ...) {
    cli::cli_ul(
        c(
            'GLOBAL ICC = {x$icc}',
            'WITHIN = {x$within}',
            'BETWEEN = {x$between}',
            'RANDOM SLOPE = {x$random_slope}',
            'PRODUCT = {x$product}'
        )
    )
    invisible(x)
}
