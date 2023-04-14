#' Check if types are valid for an effect
#'
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
#'
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



#' Add Effect sizes
#'
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
#'
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
#' @export
cross_sectional <- function() {
    c(0.05, 0.15, 0.25)
}

#' cross-sectional icc range
#' @export
longitudinal <- function() {
    c(0.40, 0.50, 0.60)
}

#' Adds effect size to `mp_base` class
#'
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


#' Prints `mp_effect` class
#'
#' @export
print.mp_effect <- function(x, ...) {
    ulid <- cli::cli_ul()
    cli::cli_li('type   = {x$type}')
    cli::cli_li('value  = {x$value}')
    cli::cli_end(ulid)
    invisible(x)
}


#' Prints `mp_effsize` class
#'
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
