#' @name mechanisms
#' @title Helper functions for producing Missing Data Mechanisms
#' @aliases mechanism MCAR MAR
#' @seealso
#' [mlmpower::power_analysis]
#' @description
#' Functions to generate data that always follows a specific mechanism
#' in accordance to a single-level model.
#' @param mis.rate A proportion for the missing data rate at population level
#' @param cause A character for a variable name that is the cause of missingness
#' @param r2 A proportion of variance explained by the cause in the missing data indicator's latent propensity
#' @param lower A logical for the lower or upper tail being more likely to be missing
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
NULL

#' @rdname mechanisms
#' @usage
#' # Generate MCAR data on outcome
#' MCAR(mis.rate)
#' @examples
#'
#' # Induce MCAR data on outcome
#' set.seed(19723)
#' model |> power_analysis(50, 5, 50, mechanism = MCAR(0.25)) -> powersim_mcar
#' @export
MCAR <- function(mis.rate) {
    if (!is.numeric(mis.rate) || length(mis.rate) != 1 || mis.rate < 0 || mis.rate > 1) {
        throw_error(
            c(x = 'The {.arg mis.rate} must be a numeric between 0 and 1.')
        )
    }
    force(mis.rate)
    function(data) {
        data[runif(NROW(data)) < mis.rate, 2] <- NA
        data
    }
}
#' @rdname mechanisms
#' @usage
#' # Generate MAR data on outcome due to `cause`
#' MAR(mis.rate, cause, r2, lower = TRUE)
#' @examples
#'
#' # Induce MAR data on outcome
#' set.seed(19723)
#' model |> power_analysis(
#'    50, 5, 50,
#'    mechanism = MAR(0.25, 'X', 0.6)
#' ) -> powersim_mar
#' @export
MAR <- function(mis.rate, cause, r2, lower = TRUE) {
    if (!is.numeric(mis.rate) || length(mis.rate) != 1 || mis.rate < 0 || mis.rate > 1) {
        throw_error(
            c(x = 'The {.arg mis.rate} must be a numeric between 0 and 1.')
        )
    }
    if (!is.numeric(r2) || length(r2) != 1 || r2 < 0 || r2 > 1) {
        throw_error(
            c(x = 'The {.arg r2} must be a numeric between 0 and 1.')
        )
    }
    if (!is.character(cause) || length(cause) != 1) {
        throw_error(
            c(x = 'The {.arg cause} must be a character vector of length 1.')
        )
    }
    if (!is.logical(lower) || length(lower) != 1) {
        throw_error(
            c(x = 'The {.arg lower} must be a logical vector of length 1.')
        )
    }
    force(cause)
    force(mis.rate)
    function(data) {
        if (!(cause %in% names(data))) throw_error(
           c(x = 'The variable "{cause}" is not in the data set.')
        )
        p <- data |> attr('parameters')
        if (length(which(names(p$mean_X) %in% cause)) != 0) {
            ind <- which(names(p$mean_X) %in% cause)
            mu <- p$mean_X[ind]
            s <- sqrt((p$phi_w[ind,ind] + p$phi_b[ind,ind]))
        }
        else if (length(which(names(p$mean_Z) %in% cause)) != 0) {
            ind <- which(names(p$mean_Z) %in% cause)
            mu <- p$mean_Z[ind]
            s <- sqrt(p$phi_b[ind + length(p$mean_X), ind + length(p$mean_X)])
        }
        else throw_error(
            c(x = 'The variable "{cause}" is not a level-1 or level-2 predictor')
        )
        val <- (data[,cause] - mu) / s
        phi <- sqrt(r2) / sqrt(1 - r2)
        M.star <- val*phi + rnorm(NROW(data))
        M <- M.star < qnorm(mis.rate, sd = sqrt(phi^2 + 1), lower.tail = lower)
        if (!lower) M <- !M
        data[M, 2] <- NA
        data
    }
}


