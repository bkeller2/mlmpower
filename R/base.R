
#' Internal command to add to model
#' @noRd
add <- function(x, y) {
    UseMethod('add', y)
}

#' Internal command to clone
#' @noRd
clone <- function(x) {
    UseMethod('clone', x)
}

#' Validate is base
#' @noRd
is.base <- function(x) {
    inherits(x, 'mp_base')
}

#' Add `mp_base` objects together
#' @noRd
#' @export
`+.mp_base` <- function(e1, e2) {
    # Check if one is missing
    if (missing(e1)) return(e2)
    if (missing(e2)) return(e1)

    # Don't allow two models
    if (is.model(e1) & is.model(e2)) {
        throw_error(c('x' = 'Cannot add two models together.'))
    }

    # Check if first is a model class
    if (is.model(e1)) return(e1 |> clone() |> add(e2))

    # Check if second is model class
    if (is.model(e2)) return(e2 |> clone() |> add(e1))

    # Otherwise perform add and return
    e1 |> add(e2)
}
