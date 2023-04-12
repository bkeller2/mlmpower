#' Check if types are valid for a variable
#'
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
#'
#' @noRd
new_variable <- function(type, name, weight, mean, sd, icc) {
    if (missing(name)) {
        cli::cli_abort('Must provide name for {.cls variable}')
    }
    if (missing(type)) {
        cli::cli_abort('Must provide variable type for {.cls variable}')
    }
    if (!valid_variable_type(type[1])) {
        cli::cli_abort('Invalid type for {.cls variable}')
    }
    if (name == '_id') {
        cli::cli_abort(c(
            'Invalid variable name ({name})',
            'x' = 'Variable names cannot be `_id`'
        ))
    }
    if (grepl('\\(|\\)', name, perl = T)) {
        cli::cli_abort(c(
            'Invalid variable name ({name})',
            'x' = 'Variable names cannot contain parentheses.'
        ))
    }
    if (grepl(' ', name, fixed = T)) {
        cli::cli_abort(c(
            'Invalid variable name ({name})',
            'x' = 'Variable names cannot contain spaces.'
        ))
    }
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
#'
#' @noRd
is.variable <- function(x) {
    inherits(x, 'mp_variable')
}

#' Check if a variable has already been specified by name
#'
#' @noRd
has_variable <- function(x, y) {
    !is.null(x[[y$name]])
}

#' Create outcome
#'
#' @export
outcome <- function(name, mean = 0, sd = 1, icc = NULL) {
    new_variable('outcome', name, NA, mean, sd, icc)
}

#' Validate classes
#'
#' @noRd
is.outcome <- function(x) {
    inherits(x, 'mp_outcome')
}

#' Create within predictor
#'
#' @export
within_predictor <- function(name, weight = 1, mean = 0, sd = 1, icc = NULL) {

    # TODO validate inputs

    if (is.null(icc)) { } # do nothing
    else if (!is.number(icc)) {
        cli::cli_abort('ICC must be a single number {.cls within_predictor}')
    }
    else if (icc < 0 | icc > 1) {
        cli::cli_abort('ICC must be between 0.0 and 1.0 in {.cls within_predictor}')
    }
    new_variable('predictor', name, weight, mean, sd, icc)
}

#' Create within predictor
#'
#' @export
within_time_predictor <- function(name, weight = 1, values) {

    # TODO validate values
    # TODO require 0 in values otherwise error
    v <- new_variable(
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

#' Create Between predictor
#'
#' @export
between_predictor <- function(name, weight = 1, mean = 0, sd = 1) {
    # TODO validate inputs
    new_variable('predictor', name, weight, mean, sd, NA)
}

#' Create Between binary predictor
#'
#' @export
between_binary_predictor <- function(name, weight = 1, proportion = 0.5) {
    # TODO validate inputs
    new_variable(
        c('binary', 'predictor'),
        name, weight,
        proportion,
        sqrt(proportion*(1.0 - proportion)),
        icc = NA
    )
}


#' Obtain levels for a variable
#'
#' @export
levels.mp_variable <- function(x) {
    if (is.null(x$icc)) return(1)
    else if (is.na(x$icc)) return(2)
    else return(1)
}

#' Adds outcome to `mp_base` class
#'
#' @noRd
`add.mp_outcome` <- function(x, y) {
    # Add as outcome if model
    if (is.model(x)) {
        if (!is.null(x$outcome)) cli::cli_abort(
            'Outcome already specified in {.cls mp_model}. '
        )
        x$outcome <- y
        return(x)
    }
    # Otherwise create model and add x
    y |> new_model() |> add(x)
}

#' Adds timevar to `mp_base` class
#'
#' @noRd
`add.mp_timevar` <- function(x, y) {
    # Add as predictor if model
    if (is.model(x)) {
        # check if predictor already there
        if (x$predictors |> has_variable(y)) cli::cli_abort(
            'A predictor variable has already specified by the name "{y$name}".'
        )
        timevars <- vapply(
            x$predictors,
            \(.)  'mp_timevar' %in% class(.),  # Select timevar
            logical(1L)
        )
        if (TRUE %in% timevars)  cli::cli_abort( c(
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
#'
#' @noRd
`add.mp_predictor` <- function(x, y) {
    # Add as predictor if model
    if (is.model(x)) {
        # check if predictor already there
        if (x$predictors |> has_variable(y)) cli::cli_abort(
            'A predictor variable has already specified by the name "{y$name}".'
        )
        x$predictors[[y$name]] <- y
        return(x)
    }
    # Otherwise construct list
    x |> add(mp_list(y))
}

#' Prints `mp_variable` class
#'
#' @export
`print.mp_variable` <- function(x, ...) {
    ulid <- cli::cli_ul()
    cli::cli_li('name   = {x$name}')
    cli::cli_li('weight = {x$weight}')
    cli::cli_li('mean   = {x$mean}')
    cli::cli_li('sd     = {x$sd}')
    if (is.null(x$icc)) {
        cli::cli_li('icc    = NULL')
    } else {
        cli::cli_li('icc    = {x$icc}')
    }
    cli::cli_end(ulid)
    invisible(x)
}



