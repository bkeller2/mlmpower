#' `mp_list` class used for storing multiple mp objects
#'
#' @noRd
mp_list <- function(...) {
    structure(list(...), class = c('mp_list', 'mp_base'))
}

#' Check if it is a  list
#' @noRd
is.mp_list <- function(x) {
    inherits(x, 'mp_list')
}

#' combine `mp_list`
#'
#' @noRd
as.mp_list <- function(object) {
    structure(object, class = c('mp_list', 'mp_base'))
}

#' Adds `mp_list` to `mp_base` class
#'
#' @noRd
add.mp_list <- function(x, y) {

    # Check if x is list and append
    if (is.mp_list(x)) {
        return(as.mp_list(c(y, x)))
    }

    # If x is an outcome create model
    if (is.outcome(x)) {
        return(x |> make_model() |> add(y))
    }

    # Check if x is model
    # Add to model and return
    if (is.model(x)) {
        # Run actions second to prevent errors
        actions <- vapply(y, is.action, logical(1))
        for (i in y[!actions]) x |> add(i) -> x
        for (i in y[actions])  x |> add(i) -> x
        return(x)
    }

    # Otherwise
    # Validate
    if (!is.base(x)) {
        throw_error('Must provide base class in {.cls mp_list}')
    }

    # Otherwise append and return list
    y[[length(y) + 1]] <- x

    # Return list
    return(y)
}


#' Wrapper for base `list`
#'
#' @noRd
add.list <- function(x, y) {
    add.mp_list(x, y)
}
