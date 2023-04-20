#' Check if types are valid for an action
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
#' @noRd
new_action <- function(type, name, weight) {
    if (missing(type)) {
        throw_error('Must provide effect size type for {.cls mp_action}')
    }
    if (!valid_action_type(type)) {
        throw_error('Invalid type for {.cls mp_action}')
    }
    if (!is.number(weight)) {
        throw_error(paste0('Weight for {.cls ', type, '} must be a single number.' ))
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

#' Check if it is an `mp_action`
#' @noRd
is.action <- function(x) {
    inherits(x, 'mp_action')
}

#' Create a Random Slope in a Model
#' @description
#' Creates a random slope that can be added to a [`mlmpower::mp_model`].
#' @param name a character string that references a variable's name
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @returns A [`mlmpower::mp_action`] that can be added to a [`mlmpower::mp_model`].
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = 0.1,
#'         within = 0.1,
#'         random_slope = 0.03
#'     )
#' )
#' # Add random slope to the model
#' model + random_slope('X')
#' @export
random_slope <- function(name, weight = 1) {
    if (!is.character(name)) {
        throw_error('Name for {.cls random slope} must be a characters')
    }

    # Handle one random slope
    if (length(name) == 1 & length(weight) == 1) {
        return(new_action('random_slope', name, weight))
    }

    # Handle multiple random slopes
    if(length(name) != length(weight)) {
        if (length(weight) == 1) {
            weight <- rep(weight, length(name))
        } else {
            throw_error('To specify multiple random slopes you need multiple weights.')
        }
    }
    for (i in weight) {
        if (i < 0 ) throw_error('Weight for {.cls random slope} cannot be negative.')
    }

    # Create output
    o <- mp_list()
    for (i in Map(random_slope, name, weight)) {
        o |> add(i) -> o
    }
    return(o)
}


#' Create a Product Term in a Model
#' @description
#' Creates a product term between two variables that can be added to a [`mlmpower::mp_model`].
#' @param name1 a character string that references the first variable's name
#' @param name2 a character string that references the second variable's name
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @details
#' Currently the product term is only limited to cross-level
#' interactions between a level-1 centered within cluster variable (`icc = 0`)
#' and level-2 variable.
#'
#' @returns A [`mlmpower::mp_action`] that can be added to a [`mlmpower::mp_model`].
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X', icc = 0.0)
#'     + between_predictor('Z')
#' )
#' # Add random slope to the model
#' model + product('X', 'Z')
#' @export
product <- function(name1, name2, weight = 1) {
    if (!is.character(name1) | !is.character(name2)) {
        throw_error('Name for {.cls product} must be a characters')
    }
    # Check names
    if (length(name1) != 1 | length(name2) != 1) {
        throw_error('Name for {.cls product} cannot be a vector of names')
    }
    # Check weights
    if (!is.number(weight)) {
        throw_error('Weight for {.cls product} must be a single numeric value')
    }
    else if (weight < 0) {
        throw_error('Weight for {.cls product} cannot be negative.')
    }

    # Create product
    new_action('product', c(name1, name2), weight)
}


#' Adds effect size to `mp_base` class
#' @noRd
`add.mp_action` <- function(x, y) {
    # Add as actions if model
    if (is.model(x)) {

        # Make sure not an outcome
        if (x$outcome$name %in% y$name) {
            throw_error(
                paste0('Cannot specify {.cls ', y$type, '} with an outcome ({y$name})')
            )
        }
        for (i in y$name) {
            if (is.null(x$predictors[[i]])) {
                throw_error(
                    paste0('Variable {i} has not been specified in {.cls ', y$type, '}.')
                )
            }
        }
        # Check products are cross-level
        if (y$type == 'product') {
            v1_level <- levels(x$predictors[[y$name[1]]])
            v2_level <- levels(x$predictors[[y$name[2]]])
            if (v1_level == v2_level) {
                throw_error(c(
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
                throw_error(c(
                    'Cross-level products require the level-1 variable to have icc = 0',
                    'x' = 'Variable {new_name[1]} has an unspecified ICC set.'
                ))
            }
            else if (x1_icc != 0 ) {
                throw_error(c(
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
#' @noRd
action_to_row <- function(x) {
    data.frame(
        `type` = switch(
            x$type,
            product      = 'Product',
            random_slope = 'Random Slope',
            NA
        ),
        `variable(s)` = paste0(x$name, collapse = ':'),
        `weight` = x$weight,
        check.names = FALSE
    )
}
