#' validates input single number
#' @noRd
is.number <- function(x) {
    length(x) == 1 & is.numeric(x)
}

#' Internal helper for multivariate normal draw with means
#' @noRd
rmvnorm <- function(n, mu, Sigma) {
    p <- NCOL(Sigma)
    if (p == 0) return(matrix(0, n, 0))
    (matrix(mu, nrow = n, ncol = p, byrow = T)
        + matrix(rnorm(p * n), ncol = p)
        %*% with(
            eigen(Sigma, symmetric = TRUE), {
                t(vectors %*% (t(vectors) * sqrt(pmax(values, 0))))
            }
        )
    )
}

#' Internal helper for multivariate normal draw without means
#' @noRd
rmvnorm_nomean <- function(n, Sigma) {
    p <- NCOL(Sigma)
    if (p == 0) return(matrix(0, n, 0))
    matrix(rnorm(p * n), ncol = p) %*% with(
        eigen(Sigma, symmetric = TRUE), {
            t(vectors %*% (t(vectors) * sqrt(pmax(values, 0))))
        }
    )
}

#' Internal helper to cwc for lme4
#' This assumes equal cluster-size
#' @noRd
make_cwc <- function(id) {
    force(id)
    function(x) {
        x - aggregate(x, by = list(id), mean)[id, 2]
    }
}

#' Internal helper to centering functions to attach to a formula
#' @noRd
centering_env <- function(id) {
    force(id)
    e <- list2env(list(
        cwc = make_cwc(id),
        cgm = \(x) x - mean(x)
    ))
}


#' Internal helper to extract lme4 results in a nice way
#' @noRd
extract_results <- function(x, ddf = 'Satterthwaite') {
    # Extract regression coefficients
    reg_coef <- coefficients(summary(x, ddf = ddf))
    # Extract variance/covariance
    vc <- lme4::VarCorr(x)
    vars <- diag(vc$`_id`)

    # Handle covars
    covs <- with(vc, `_id`[lower.tri(`_id`)])
    names(covs) <- unlist(
        lapply(seq_along(vars), \(i) {
            vapply(seq_along(vars)[-seq_len(i)], \(j) {
                paste0( 'cov[', names(vars)[i],', ', names(vars)[j], ']')
            }, character(1L))
        })
    )

    # Return parameters
    list(
        fixed = reg_coef[, 1],
        random = c(vars, covs),
        within_variance = sigma(x)^2
    )
}


#' Internal helper for computing jackknife standard errors
#' @noRd
jack_se <- function(data, func) {
    jack <- vapply(seq_len(NCOL(data)), \(i) func(data[,-i, drop = F]), numeric(NROW(data)))
    v <- ((NCOL(jack) -1) / NCOL(jack)) * rowSums((jack - rowMeans(jack))^2)
    sqrt(v)
}

#' Internal helper for margin of error with jackknife
#' @noRd
mc_error <- function(data, func, alpha = 0.05) {
    jack_se(data, func) * qnorm(1 - alpha/2)
}


#' Internal function for suppressing output
#' @noRd
quiet <- function(x) {
    invisible(suppressMessages(suppressWarnings(force(x))))
}

#' Internal function for parsing parameter names
#' @noRd
parse_param_names <- function(x) {
    vapply(strsplit(x, '.', fixed = T), \(x) {
        if (x[1] == 'fixed') {
            paste0('Fixed:  ', paste0(x[-1], collapse = '.'))
        }
        else {
            paste0('Random: Slopes Test')
        }
    }, character(1L))
}

#' Internal function for pretty output printing of power tables
#' @noRd
power_table <- function(x) {
    d <- data.frame(
        Power = paste(
            formatC(x[,1], digits = 2, format = 'f' ),
            paste('\u00B1', formatC(x[,2], digits = 2, format = 'f'))
        )
    )
    rownames(d) <- parse_param_names(row.names(x))
    print(d, right = F)
    invisible(d)
}


#' Internal function to keep online mean
#' @noRd
online_mean <- function(e1, e2) {

    # Check if one is missing
    if (missing(e1)) return(e2)
    if (missing(e2)) return(e1)

    # Check they are same classes
    if (!identical(class(e1), class(e2))) {
        throw_error('The two must have identical classes')
    }
    e1_count <- attr(e1, '_count')
    e2_count <- attr(e2, '_count')

    # Return Reverse
    if (is.null(e1_count) && !is.null(e2_count)) {
        return(online_mean(e2, e1))
    }
    # Construct new one
    else if (is.null(e1_count) && is.null(e2_count)) {
        e1 <- (e1 + e2) / 2
        attr(e1, '_count') <- 2
    }
    # Handle both averages
    else if (!is.null(e1_count) && !is.null(e2_count)) {
        d <- e1 - e2
        total <- e1_count + e2_count
        e1 <- (e1_count * e1  + e2_count * e2) / (total);
        attr(e1, '_count') <- total
    }
    else {
        e1_count <- e1_count +  1
        e1 <- e1 + (e2 - e1) / (e1_count)
        attr(e1, '_count') <- e1_count
    }
    # otherwise handle adding one element
    return(e1)
}

#' Internal function to take diagonal safely even if its one element
#' @noRd
diagonal <- function(x) {
    if (is.matrix(x)) diag(x)
    else diag(x, NROW(x), NROW(x))
}

#' Internal map function that uses all combinations of inputs
#' @noRd
cmap <- function(f, ...) {
    l <- expand.grid(..., stringsAsFactors=FALSE)
    do.call(mapply, c(FUN = f, l, SIMPLIFY = FALSE))
}

#' Internal error function
#' Wrapper for `cli_abort` to not specify the call
#' @noRd
throw_error <- function(message, ..., .envir = parent.frame(), .frame = .envir) {
    cli::cli_abort(message, ..., .envir = .envir, .frame = .frame, call = NULL)
}

