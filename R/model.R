#' Generate Multilevel base class
#'
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

#' Check if a `mp_model` is properly constructed
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

#' Subset by ICC `mp_model`
#'
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

#' Internal function to call function for a `mp_model`
#' @noRd
with_model <- function(model, expr, ...) {
    if (!is.function(expr)) {
        throw_error('Second argument must be function')
    }
    environment(expr) <- model
    expr(...)
}



#' Prints `mp_model` object
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
            lapply(x$predictors, as.data.frame)
        )[, c('level', 'icc', 'weight', 'type')],
        right = FALSE
    )

    if (length(x$actions) > 0) {
        cli::cli_h3('{.cli model specifications}')
        cli::cli_text('')
        print(
            do.call(
                'rbind',
                lapply(x$actions, as.data.frame)
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


#' Provides parameter summary of a `mp_model` object
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
