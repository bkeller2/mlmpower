
#' Check if types are valid for an action
#'
#' @noRd
valid_action_type <- function(x) {
    switch(
        x,
        product = TRUE,
        random_slope = TRUE,
        FALSE # Default
    )
}

#' Generate Base action for mlmpower
#'
#' @noRd
new_action <- function(type, name, weight) {
    if (missing(type)) {
        cli::cli_abort('Must provide effect size type for {.cls mp_action}')
    }
    if (!valid_action_type(type)) {
        cli::cli_abort('Invalid type for {.cls mp_action}')
    }
    if (!is.number(weight)) {
        cli::cli_abort(paste0('Weight for {.cls ', type, '} must be a single number.' ))
    }

    # Create new effect sizes
    structure(
        list2env(
            list(type = type, name = name, weight = weight),
            parent = emptyenv()
        ),
        class = c('mp_action', 'mp_base')
    )
}



#' Add a random slope
#'
#' @export
random_slope <- function(name, weight = 1) {
    if (!is.character(name)) {
        cli::cli_abort('Name for {.cls random slope} must be a characters')
    }

    # Handle one random slope
    if (length(name) == 1 & length(weight) == 1) {
        return(new_action('random_slope', name, weight))
    }

    # Handle multiple random slopes
    if( length(name) != length(weight)) {
        if (length(weight) == 1) {
            weight <- rep(weight, length(name))
        } else {
            cli::cli_abort('To specify multiple random slopes you need multiple weights.')
        }
    }

    # Create output
    o <- mp_list()
    for (i in Map(random_slope, name, weight)) {
        o |> add(i) -> o
    }
    return(o)
}


#' Add an interaction effect
#'
#' @export
product <- function(name1, name2, weight = 1) {
    if (!is.character(name1) | !is.character(name2)) {
        cli::cli_abort('Name for {.cls product} must be a characters')
    }
    # Check names
    if (length(name1) != 1 | length(name2) != 1) {
        cli::cli_abort('Name for {.cls product} cannot be a vector of names')
    }
    # Check weights
    if (!is.number(weight)) {
        cli::cli_abort('Weight for {.cls product} must be a single numeric value')
    }

    # Create product
    new_action('product', c(name1, name2), weight)
}


#' Adds effect size to `mp_base` class
#'
#' @noRd
`add.mp_action` <- function(x, y) {
    # Add as actions if model
    if (is.model(x)) {

        # Make sure not an outcome
        if (x$outcome$name %in% y$name) {
            cli::cli_abort(
                paste0('Cannot specify {.cls ', y$type, '} with an outcome ({y$name})')
            )
        }
        for (i in y$name) {
            if (is.null(x$predictors[[i]])) {
                cli::cli_abort(
                    paste0('Variable {i} has not been specified in {.cls ', y$type, '}.')
                )
            }
        }
        # Check products are cross-level
        if (y$type == 'product') {
            v1_level <- levels(x$predictors[[y$name[1]]])
            v2_level <- levels(x$predictors[[y$name[2]]])
            if (v1_level == v2_level) {
                cli::cli_abort(c(
                    'Only cross-level interactions are supported.',
                    "x" = 'Variables {y$name} are at the same level.'
                ))
            }
            # Make first name level-1
            new_name <- vector('character', 2L)
            new_name[v1_level] <- y$name[1]
            new_name[v2_level] <- y$name[2]
            y$name <- new_name
            # Verify level-1 variable has no icc
            x1_icc <- x$predictors[[new_name[1]]]$icc
            if (is.null(x1_icc)) {
                cli::cli_abort(c(
                    'Cross-level products require the level-1 variable to have icc = 0',
                    'x' = 'Variable {new_name[1]} has an unspecified ICC set.'
                ))
            }
            else if (x1_icc != 0 ) {
                cli::cli_abort(c(
                    'Cross-level products require the level-1 variable to have icc = 0',
                    'x' = 'Variable {new_name[1]} has an ICC = {x1_icc}'
                ))
            }
        }
        x$actions[[length(x$actions) + 1]] <- y
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}


#' Converts a `mp_action` class into a single row `data.frame`
#' @export
as.data.frame.mp_action <- function(x) {
    data.frame(
        `type` = switch(
            x$type,
            product      = 'Product',
            random_slope = 'Random Slope',
            NA
        ),
        `variable(s)` = paste0(x$name, collapse = ':'),
        `weight` = x$weight
    )
}
