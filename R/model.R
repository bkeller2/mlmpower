#' Generate Multilevel base class
#'
#' @noRd
new_model <- function(o) {
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

#' Validate model
#' @noRd
is.model <- function(x) {
    inherits(x, 'mp_model')
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

#' Subset by ICC `mp_model`
#'
#' @export
subset.mp_model <- function(x, icc) {
    if (!is.number(icc)) {
        cli::cli_abort( 'icc needs to be a single number')
    }
    else if ( !(icc %in% x$effect_size$icc)) {
        cli::cli_abort( 'icc is not in {.cli mp_model}')
    }
    new_x <- clone(x)
    new_x$effect_size <- clone(x$effect_size)
    new_x$effect_size$icc <- icc
    return(new_x)
}

#' With binding for `mp_model`
#' @export
with.mp_model <- function(data, expr, ...) {
    if (!is.function(expr)) {
        cli::cli_abort('Second argument must be function')
    }
    environment(expr) <- data
    expr(...)
}




# TODO need way to validate model is correctly specified

#' Prints `mp_model` object
#' @export
print.mp_model <- function(x, ...) {
    cli::cli_h2('{.cli mlmpower} model specification')
    cli::cli_h3('{.cli outcome}')
    print(x$outcome)
    cli::cli_h3('{.cli predictors}')
    cli::cli_ol()
    for (i in x$predictors) {
        cli::cli_li('Predictor:')
        print(i)
        cli::cli_text('')
    }
    cli::cli_end()

    cli::cli_h3('{.cli effect sizes}')
    cli::cli_ul()
    cli::cli_li('ICC = {x$effect_size$icc}')
    cli::cli_li('WITHIN = {x$effect_size$within}')
    cli::cli_li('BETWEEN = {x$effect_size$between}')
    cli::cli_li('RANDOM SLOPE = {x$effect_size$random_slope}')
    cli::cli_li('PRODUCT = {x$effect_size$product}')
    cli::cli_end()

    # TODO print actions

    # TODO need way to deal with functions printing
    # cli::cli_h3('{.cli correlations}')
    # cli::cli_ul()
    # cli::cli_li('WITHIN = {x$corrs$within}')
    # cli::cli_li('BETWEEN = {x$corrs$between}')
    # cli::cli_li('RANDOM EFFECT = {x$corrs$randeff_cor}')
    # cli::cli_end()
}


#' Provides parameter summary of a `mp_model` object
#' TODO improvce printing of parameters summary
#' HANDLE random values
#' @export
summary.mp_model <- function(object, ...) {

    # Get icc
    icc <- object$effect_size$icc

    # Handle multiple ICCs
    if (length(icc) > 1) {

        # Create names
        names(icc) <- paste0('icc = ', icc)

        # Run simulation
        results <- lapply(icc, \(x) {
            model |> subset(icc = x) |> new_parameters()
        })

    } else {
        results <- object |> new_parameters()
    }

    results
}
