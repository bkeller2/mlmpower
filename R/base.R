
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
        throw_error('Cannot add two model sets together.')
    }
    # Check if second is models class or outcome
    if (is.model(e2)) return(e2 |> add(e1))

    # Otherwise check if one is a list
    e1 |> add(e2)
}
