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

#' Check if `mp_model`
#' @noRd
is.model <- function(x) {
    inherits(x, 'mp_model')
}

#' Validate `mp_model`
#' @noRd
is.valid <- function(x) {
    # TODO check that this is correctly specified
    #  - [ ] Check its a `mp_model` object
    #  - [ ] Check if it has outcome
    #  - [ ] Check if it has effect size
    #  - [ ] Check valid effect size
    #    - (icc, proper predictors based on r2, etc)

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

#' Subset by ICC `mp_model`
#'
#' @export
subset.mp_model <- function(x, icc) {

    # Validate model first
    is.valid(x)

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



#' Prints `mp_model` object
#' @export
print.mp_model <- function(x, ...) {

    # Validate model first
    is.valid(x)

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
            lapply(model$predictors, as.data.frame)
        )[, c('level', 'icc', 'weight', 'type')],
        right = FALSE
    )

    if (length(model$actions) > 0) {
        cli::cli_h3('{.cli model specifications}')
        cli::cli_text('')
        print(
            do.call(
                'rbind',
                lapply(model$actions, as.data.frame)
            ),
            right = FALSE,
            row.names = FALSE
        )
    }

    cli::cli_h3('{.cli effect sizes}')
    cli::cli_text('')
    cli::cli_ul(
        c(
            'GLOBAL ICC = {x$effect_size$icc}',
            'WITHIN = {x$effect_size$within}',
            'BETWEEN = {x$effect_size$between}',
            'RANDOM SLOPE = {x$effect_size$random_slope}',
            'PRODUCT = {x$effect_size$product}'
        )
    )

    cli::cli_h3('{.cli correlations}')
    cli::cli_text('')
    cli::cli_ul(
        c(
            'WITHIN = {as.character(x$corrs$within_cor)}',
            'BETWEEN = {as.character(x$corrs$between_cor)}',
            'RANDOM EFFECT = {as.character(x$corrs$randeff_cor)}'
        )
    )
    cli::cli_end()
}


#' Provides parameter summary of a `mp_model` object
#' @export
summary.mp_model <- function(object, ...) {

    # validate model
    is.valid(object)

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
        results <- lapply(icc, \(x) {
            model |> subset(icc = x) |> new_parameters_mean()
        })

    } else {
        results <- object |> new_parameters_mean()
    }

    results # Return results
}
