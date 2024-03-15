#' @rdname mp_variable
#' @title Functions for Creating Variables
#' @name Variables
#' @aliases variables variable mp_variable
#' @description
#' These functions are the building blocks used to create the multilevel model
#' and are used to specify the names, properties, and variable types.
#' @returns Returns a `mp_variable` object based on the variable's type.
#' @details
#' Note that specifying an `icc = 0` in `within_predictor()`
#' will result in a centered within cluster (CWC) predictor.
#'
#' See vignettes for more details.
#' ```
#' vignette(package = 'mlmpower')
#' ```
NULL

#' Check if types are valid for a variable
#' @noRd
valid_variable_type <- function(x) {
    switch(
        x,
        outcome = TRUE,
        predictor = TRUE,
        binary = TRUE,
        timevar = TRUE,
        FALSE # Default
    )
}

#' Generate Base variable for mlmpower
#' @noRd
make_variable <- function(type, name, weight, mean, sd, icc) {
    if (missing(name)) {
        throw_error('Must provide name for {.cls variable}.')
    }
    if (missing(type)) {
        throw_error('Must provide variable type for {.cls variable}.')
    }
    if (!valid_variable_type(type[1])) {
        throw_error('Invalid type for {.cls variable}.')
    }
    if (name == '_id') {
        throw_error(c(
            'Invalid variable name ({name}).',
            'x' = 'Variable names cannot be `_id`'
        ))
    }
    if (name == 'all') {
        throw_error(c(
            'Invalid variable name ({name}).',
            'x' = 'Variable names cannot be `all`'
        ))
    }
    if (grepl('\\(|\\)', name, perl = T)) {
        throw_error(c(
            'Invalid variable name ({name}).',
            'x' = 'Variable names cannot contain parentheses.'
        ))
    }
    if (grepl(' ', name, fixed = T)) {
        throw_error(c(
            'Invalid variable name ({name}).',
            'x' = 'Variable names cannot contain spaces.'
        ))
    }
    if (type[1] != 'outcome') {
        if (!is.number(weight)) throw_error(
            'The weight must be a single number for {.cls variable} ({name}).'
        )
        else if (weight < 0) throw_error(
            'The weight must be a positive number for {.cls variable} ({name}).'
        )
    }
    if (!is.number(mean)) throw_error(
        'The mean must be a single number for {.cls variable} ({name}).'
    )
    if (!is.number(sd)) throw_error(
        'The standard deviation must be a single number for {.cls variable} ({name}).'
    )
    if (!is.null(icc)) {
        if (!is.na(icc)) {
            if (!is.number(icc)) {
                throw_error(c(
                    'Invalid ICC for {.cls variable} ({name}).',
                    'x' = 'It must be a single number.'
                ))
            }
            else if (icc < 0 | icc > 1) {
                throw_error(c(
                    'Invalid ICC for {.cls variable} ({name}).',
                    'x' = 'It must be between 0 and 1.'
                ))
            }
        }
    }
    if (type[1] != 'outcome') {
        if (weight < 0) throw_error(c(
            'Invalid Weight for {.cls variable} ({name}).',
            'x' = 'Weights cannot be negative.'
        ))
    }

    # Return variable
    structure(
        list(
            name = name,
            weight = weight,
            mean = mean,
            sd = sd,
            icc = icc
        ), class = c( paste0('mp_', type), 'mp_variable', 'mp_base')
    )
}

#' Validate classes
#' @noRd
is.variable <- function(x) {
    inherits(x, 'mp_variable')
}

#' Check if a variable has already been specified by name
#' @noRd
has_variable <- function(x, y) {
    !is.null(x[[y$name]])
}

# Create outcome
#' @rdname mp_variable
#' @aliases outcome
#' @param name a character string for the specific variable's name
#' @param mean a single numeric value that specifies the variable's mean
#' @param sd a single numeric value that specifies the variable's standard deviation
#' @param icc a single numeric value between 0 and 1 that specifies the variable's intraclass correlation.
#' If `NULL` then the global ICC specified in [mlmpower::effect_size()] is used instead.
#' @export
outcome <- function(name, mean = 10, sd = 5, icc = NULL) {
    make_variable('outcome', name, NA, mean, sd, icc)
}

#' Validate classes
#' @noRd
is.outcome <- function(x) {
    inherits(x, 'mp_outcome')
}

# Create within predictor
#' @rdname mp_variable
#' @aliases within_predictor
#' @param name a character string for the specific variable's name
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @param mean a single numeric value that specifies the variable's mean
#' @param sd a single numeric value that specifies the variable's standard deviation
#' @param icc a single numeric value between 0 and 1 that specifies the variable's intraclass correlation.
#' If `NULL` then the global ICC specified in [`mlmpower::effect_size()`] is used instead.
#' @export
within_predictor <- function(name, weight = 1, mean = 0, sd = 1, icc = NULL) {
    if (is.null(icc)) { } # do nothing
    else if (!is.number(icc)) {
        throw_error('ICC must be a single number {.cls within_predictor}')
    }
    else if (icc < 0 | icc > 1) {
        throw_error('ICC must be between 0.0 and 1.0 in {.cls within_predictor}')
    }
    make_variable('predictor', name, weight, mean, sd, icc)
}

# Create within time predictor
#' @rdname mp_variable
#' @aliases within_time_predictor
#' @param name a character string for the specific variable's name
#' @param values a numeric vector specifying the time scores that will be repeated within each cluster.
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @export
within_time_predictor <- function(name, values, weight = 1) {

    if (!is.numeric(values)) throw_error(
        'Values must be a numeric vector in {.cls within_time_predictor}.'
    )

    v <- make_variable(
        c('timevar', 'predictor'),
        name, weight,
        mean(values),
        # Compute stddev manually
        sqrt(mean((values - mean(values))^2)),
        icc = 0.0
    )
    v$values <- values
    return(v)
}

# Create Between predictor
#' @rdname mp_variable
#' @aliases between_predictor
#' @param name a character string for the specific variable's name
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @param mean a single numeric value that specifies the variable's mean
#' @param sd a single numeric value that specifies the variable's standard deviation
#' @export
between_predictor <- function(name, weight = 1, mean = 0, sd = 1) {
    make_variable('predictor', name, weight, mean, sd, NA)
}

# Create Between binary predictor
#' @rdname mp_variable
#' @aliases between_binary_predictor
#' @param name a character string for the specific variable's name
#' @param proportion a single numeric value between 0 and 1 that specifies the proportion of 1's at the population.
#' @param weight a single numeric value specifying the variable's contribution to the variance explained metric.
#' Weights are normalized across all variables of the same level.
#' @export
between_binary_predictor <- function(name, proportion = 0.5, weight = 1) {
    if (!is.number(proportion)) {
        throw_error('Proportion must be a single number {.cls between_binary_predictor}')
    }
    else if (proportion < 0 | proportion > 1) {
        throw_error('Proportion must be between 0.0 and 1.0 in {.cls between_binary_predictor}')
    }
    make_variable(
        c('binary', 'predictor'),
        name, weight,
        proportion,
        sqrt(proportion*(1.0 - proportion)),
        icc = NA
    )
}


#' Obtain Level of Observation for a Variable
#' @description
#' Returns which level a variable is observed at in the multilevel model.
#' @param x a [`mlmpower::mp_variable`].
#' @returns Returns a single integer of the level of observation
#' @examples
#' # Returns 1
#' levels(
#'     within_predictor(
#'         'X',
#'         weight = 1,
#'         mean = 5,
#'         sd = 10,
#'         icc = 0.1
#'     )
#' )
#' @export
levels.mp_variable <- function(x) {
    if (is.null(x$icc)) return(1)
    else if (is.na(x$icc)) return(2)
    else return(1)
}

#' Adds outcome to `mp_base` class
#' @noRd
add.mp_outcome <- function(x, y) {
    # Add as outcome if model
    if (is.model(x)) {
        if (!is.null(x$outcome)) throw_error(
            'Outcome already specified in {.cls mp_model}. '
        )
        x$outcome <- y
        return(x)
    }
    # Otherwise create model and add x
    y |> make_model() |> add(x)
}

#' Adds timevar to `mp_base` class
#' @noRd
add.mp_timevar <- function(x, y) {
    # Add as predictor if model
    if (is.model(x)) {
        # check if predictor already there
        if (x$predictors |> has_variable(y)) throw_error(
            'A predictor variable has already specified by the name "{y$name}".'
        )
        timevars <- vapply(
            x$predictors,
            \(.)  'mp_timevar' %in% class(.),  # Select timevar
            logical(1L)
        )
        if (TRUE %in% timevars)  throw_error( c(
            'A time variable has already been specified.',
            'x' = 'Only one time variable is allowed in the model.',
            'i' = 'Problem variable is {y$name}.'
        ))
        x$predictors[[y$name]] <- y
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}

#' Adds predictor to `mp_base` class
#' @noRd
add.mp_predictor <- function(x, y) {
    # Add as predictor if model
    if (is.model(x)) {
        # check if predictor already there
        if (x$predictors |> has_variable(y)) throw_error(
            'A predictor variable has already specified by the name "{y$name}".'
        )
        x$predictors[[y$name]] <- y
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}

#' Converts a `mp_variable` class into a single row `data.frame`
#' @noRd
variable_to_row <- function(x) {
    with(x,
         data.frame(
             name = name,
             weight = weight,
             mean = mean,
             sd = sd,
             icc = if (is.null(icc)) '' else if (is.na(icc)) 'NA' else icc,
             type = if ('mp_timevar' %in% class(x)) 'timevar'
             else if ('mp_binary' %in% class(x)) 'binary'
             else 'continuous',
             level = levels(x)
         )
    )
}

#' Prints a [`mlmpower::mp_variable`]
#' @description
#' Prints a [`mlmpower::mp_variable`] in a human readable format.
#' @param x a [`mlmpower::mp_variable`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @examples
#' print(
#'     within_predictor(
#'         'X',
#'         weight = 1,
#'         mean = 5,
#'         sd = 10,
#'         icc = 0.1
#'     )
#' )
#' @export
print.mp_variable <- function(x, ...) {
    cli::cli_ul(
        c(
            'name   = {x$name}',
            'weight = {x$weight}',
            'mean   = {x$mean}',
            'sd     = {x$sd}',
            if (is.null(x$icc))
                'icc    = NULL'
            else
                'icc    = {x$icc}'
        )
    )
    invisible(x)
}


