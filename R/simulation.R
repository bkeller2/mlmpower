#' Generates Data Sets Based on a [`mlmpower::mp_model`]
#' @aliases mp_data
#' @description
#' Generates data sets based on a [`mlmpower::mp_model`].
#' These data sets will be returned as a [`data.frame`] and include the
#' solved parameters as an attribute to the [`data.frame`].
#' @param model a [`mlmpower::mp_model`].
#' @param n_within a single positive integer of the desired within cluster observations.
#' @param n_between a single positive integer of the desired between cluster observations.
#' @param ndata a single positive integer of the number of desired data sets.
#' @details
#' Note that there must only be one global ICC in [`mlmpower::mp_model`].
#' @returns
#' For `ndata = 1` a single [`data.frame`] is returned.
#' The first variable is the cluster identifier labeled `_id`.
#' This object is also of class `mp_data` which means that it was generated
#' based on a specified model.
#' If multiple data sets are requested then they will be contained in a [`list`].
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
#' # Set seed
#' set.seed(198723)
#' # Create data set
#' model |> generate(5, 50) -> mydata
#' @export
generate <- function(model, n_within, n_between, ndata = 1) {

    # Validate Inputs
    model |> is_valid()

    if (!is.number(n_within)) throw_error(
        "{.arg n_within} must be a single integer >= 1."
    )
    if (!is.number(n_between)) throw_error(
        "{.arg n_between} must be a single integer >= 1."
    )
    if (!is.number(ndata)) throw_error(
        "{.arg ndata} must be a single integer >= 1."
    )
    if (n_within < 1) throw_error(
        "{.arg n_within} must be a single integer >= 1."
    )
    if (n_between < 1) throw_error(
        "{.arg n_between} must be a single integer >= 1."
    )
    if (ndata < 1) throw_error(
        "{.arg ndata} must be a single integer >= 1."
    )
    if (length(model$effect_size$icc) != 1) throw_error(
        "{.arg model} must have only one global ICC specified."
    )

    # Obtain levels
    lvls <- vapply(model$predictors, levels, numeric(1L))

    # Obtain timevar level-1 indicators
    timevar_l1 <- unlist(lapply(
        model$predictors[lvls == 1],  # Subset level-1
        \(.)  'mp_timevar' %in% class(.)  # Select timevar
    ))
    if (sum(timevar_l1) > 1) throw_error("Only one timevar should be specified")

    # Change within sample size if timevar specified
    if (TRUE %in% timevar_l1) {
        len <- length(model$predictors[lvls == 1][[which(timevar_l1)]]$values)
        if (missing(n_within)) {
            n_within <- len
        } else if (len != n_within) {
            cli::cli_alert_info(
                "Setting n_within = {len} to match time variable's length"
            )
            n_within <- len
        }
    }

    # Return list if ndata isn't 1
    if (ndata > 1) {
        return(
            lapply(seq_len(ndata), \(.) {
                model |> generate(n_within, n_between, ndata = 1)
            })
        )
    }

    # Obtain binary level-2 indicators
    binary_l2 <- sapply(
        # Subset level-2
        model$predictors[lvls == 2],  # Subset level-2
        \(.)  'mp_binary' %in% class(.)   # Select binary
    )

    # Create parameters
    model |> make_parameters() -> p

    # useful precomputed values
    N <- n_within * n_between # Total sample size
    l1 <- length(p$mean_X)
    l2 <- length(p$mean_Z)

    # Generate ID variable
    `_id` <- seq_len(n_between) %x% matrix(1, n_within)


    # Generate X predictors

    # Generate Level-1 Time variables if exists
    if (TRUE %in% timevar_l1) {
        # Get values
        vals <- model$predictors[lvls == 1][[which(timevar_l1)]]$values

        # Create X_w matrix
        X_w <- matrix(0, N, l1)
        X_w[,  timevar_l1] <-  matrix(vals, nrow = N, ncol = 1)

        # Generate X's conditional on timevars
        if (l1 > 1) {
            cond_var <- (
                p$phi_w[ !timevar_l1, !timevar_l1, drop = F]
                - p$phi_w[!timevar_l1, timevar_l1, drop = F]
                %*% solve(
                    p$phi_w[timevar_l1, timevar_l1, drop = F],
                    p$phi_w[timevar_l1, !timevar_l1, drop = F]
                )
            )
            cond_mean <- t(
                p$phi_w[!timevar_l1, timevar_l1, drop = F]
                %*% solve(
                    p$phi_w[timevar_l1, timevar_l1, drop = F],
                    t(X_w[,  timevar_l1] - p$mean_X[timevar_l1, drop = F])
                )
            )
            X_w[, !timevar_l1] <- cond_mean + rmvnorm_nomean(N, cond_var)
        }

        # Draw between based on conditionals
        s <- c(!timevar_l1, rep(TRUE, l2)) # level-2 selector for timevar
        X_b <- matrix(0, N, l1 + l2)
        X_b[, s] <- rmvnorm(
            n = N,
            mu = c(p$mean_X[!timevar_l1], rep(0, l2)),
            Sigma = p$phi_b[s,s, drop = F]
        )[`_id`, , drop = F]

    } else {
        X_w <- rmvnorm_nomean(N, p$phi_w)
        X_b <- rmvnorm(N, c(p$mean_X, rep(0, l2)), p$phi_b)[`_id`, , drop = F]
    }

    # Generate Binary Between Predictors
    if (TRUE %in% binary_l2) {
        bin_select <- seq_len(l2)[binary_l2] + l1
        X_binary <- X_b[ , bin_select, drop = F]
        # Obtain thresholds
        thresh <- qnorm(p$mean_Z[binary_l2], sd = sqrt(diag(p$phi_b)[bin_select]))
        # Generate data
        for (i in seq_along(thresh)) {
            X_b[ , bin_select[i]] <- ifelse(
                X_b[ , bin_select[i]] < thresh[i],
                0.0, 1.0
            ) - p$mean_Z[seq_len(l2)[i]] # Subtract mean for centering
        }
    }

    # Create all interactions
    inter <- do.call('cbind', apply(
        X_w, 2,
        \(.). * X_b[, seq_len(l2) + l1, drop = F],
        simplify = F
    ))

    # Create predictors matrix
    X <- cbind(1, X_w, inter, X_b)

    # Generate level-1 residuals
    e_ij <- rnorm(N, 0, sqrt(p$var_e))

    # Generate level-2 residuals
    u_j <- rmvnorm_nomean(n_between, p$tau)[`_id`, , drop = F]

    # Generate Outcome
    Y <- X %*% p$gammas + rowSums(X[,seq_len(l1 + 1), drop = F] * u_j) + e_ij

    # Return data set
    d <- data.frame(
        `_id`,
        Y,
        X_w + X_b[ , seq_len(l1), drop = F],
        (X_b[ , seq_len(l2) + l1, drop = F]
         # Add means for Z
         + matrix(
             p$mean_Z,
             nrow = N,
             ncol = l2,
             byrow = T
         )
        )
    )
    names(d) <- c(
        '_id',
        names(p$mean_Y),
        names(which(vapply(model$predictors, levels, numeric(1L)) == 1)),
        names(which(vapply(model$predictors, levels, numeric(1L)) == 2))
    )
    # Add timevar_l1 to parameters
    p$timevar_l1 <- timevar_l1
    attr(d, 'parameters') <- p

    # Return data
    structure(
        d,
        class = c('mp_data', 'data.frame')
    )
}

#' Check if it is a  list
#' @noRd
is.mp_data <- function(x) {
    inherits(x, 'mp_data')
}


#' Analyzes a single [`mlmpower::mp_data`] using [`lme4::lmer`]
#' @description
#' Analyzes a single [`mlmpower::mp_data`] based on the data generating model.
#' @param data a [`mlmpower::mp_data`].
#' @param alpha the significance level to determine if an effect is statistically significant.
#' If `NULL` then no nested model testing is conducted.
#' @param no_lrt do not perform additional likelihood ratio tests.
#' Setting to `TRUE` will speed up the analysis because the model is only fit once.
#' @param ... other arguments passed to [`lme4::lmer()`].
#' @returns A [`list`] that with the following named elements:
#' - `estimates`: The estimates from fitting the model.
#' - `sig_test`: The logical if the estimates were statistically siginificant based on `alpha`.
#' - `parameters`: The [`mlmpower::mp_parameters`] extracted from `data`.
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
#' # Set seed
#' set.seed(198723)
#' # Create data set and analyze
#' model |> generate(5, 50) |> analyze() -> results
#' @export
analyze <- function(data, alpha = 0.05, no_lrt = FALSE, ...) {

    # Check Inputs
    if (!is.mp_data(data)) throw_error(
        "{.arg data} must be of a {.cli mp_data} object."
    )
    # Allow null
    if (!is.null(alpha)) {
        if (!is.number(alpha)) throw_error(
            "{.arg alpha} must be a single number between 0 and 1."
        )
        if (alpha <= 0 | alpha >= 1) throw_error(
            "{.arg alpha} must be a single number between 0 and 1."
        )
    }
    if (!is.logical(no_lrt) & length(no_lrt) == 1) throw_error(
        "{.arg no_lrt} must be a single logical value"
    )

    # Make centering env
    e <- centering_env(data$`_id`)

    # Get formulas for model
    attr(data, 'parameters') |> to_formula(e) -> f
    attr(data, 'parameters') |> to_formula(e, nested = T) -> f_nest

    # Fit Model with lme4 and return results
    full_reml <- quiet(lmerTest::lmer(f, data, ...))

    # Obtain LRT
    if (!is.null(alpha) & !no_lrt & !identical(f, f_nest)) {
        full_ml <- quiet(lme4::lmer(f, data, REML = FALSE, ...))
        nested <- quiet(lme4::lmer(f_nest, data, REML = FALSE, ...))
        lrt <- varTestnlme::varCompTest(full_ml, nested, pval.comp = "bounds", output = FALSE)
        rand_result <- c(omnibus_test = lrt$p.value[[4]] < alpha)
    } else {
        rand_result <- c(omnibus_test = NA)
    }

    # Return output
    list(
        estimates = full_reml |> extract_results(),
        sig_test = if (is.null(alpha)) list() else list(
            fixed  = coefficients(summary(full_reml))[,'Pr(>|t|)'] < alpha,
            random = rand_result
        ),
        parameters =  attr(data, 'parameters')
    )
}


#' Internal power analysis function based on `mp_model`
#' @noRd
`_power_analysis` <- function(
        model,
        replications,
        n_within,
        n_between,
        ...) {

    # Double check valid inputs
    if (n_within < 1) throw_error(
        "{.arg n_within} must be a single integer >= 1."
    )
    if (n_between < 1) throw_error(
        "{.arg n_between} must be a single integer >= 1."
    )

    # Run reps and include progress bar
    results <- lapply(
        cli::cli_progress_along(
            seq_len(replications),
            clear = FALSE,
            format = paste0(
                "  Simulating {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
                "| Elapsed: {cli::pb_elapsed}",
                "| ETA: {cli::pb_eta}"
            )
        ), \(.) {
            model |> generate(n_within, n_between) |> analyze(...)
        }
    )

    # Extract results
    e <- sapply(results, \(.) unlist(.$estimates), simplify = 'array')
    s <- sapply(results, \(.) unlist(.$sig_test), simplify = 'array')
    p <- Reduce(average, lapply(results, \(.) .$parameters))

    # Drop omnibus test if missing
    s <- na.omit(s)

    # Output table
    structure(
        list(
            sim = list(
                model = model,
                replications = replications,
                n_within = n_within,
                n_between = n_between
            ),
            power = cbind(
                value = rowMeans(s),
                mc_moe = mc_error(s, rowMeans)
            ),
            estimates = cbind(
                mean       = rowMeans(e),
                sampling_sd = apply(e, 1, sd),
                mc_moe = mc_error(e, rowMeans)
            ),
            mean_parameters = p |> clean_parameters()
        ),
        class = c('mp_power', 'mp_base')
    )

}


#' Conduct a Power Analysis Based on [`mlmpower::mp_model`]
#' @aliases mp_power
#' @description
#' This function will construct a multilevel power analysis via a Monte Carlo Simulation
#' based on a constructed [`mlmpower::mp_model`].
#' @param model a [`mlmpower::mp_model`].
#' @param replications a single positive integer of the number of replications per condition.
#' @param n_within an integer vector of the desired within cluster observations.
#' @param n_between an integer vector of the desired between cluster observations.
#' @param ... other arguments passed to [`mlmpower::analyze()`].
#' @details
#' Specifying multiple `n_within` and `n_between` will produce a full factorial simulation design.
#' @returns A `mp_power` object that contains the results.
#' See [`mlmpower::print.mp_power`] for more information.
#' The object has the following slots:
#' - `sim`: The information about the simulation
#' - `power`: The power power results per condition.
#' - `estimates`: The simulation summaries of the parameter estimates per condition.
#' - `mean_parameters`: The average population parameter per condition.
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
#' # Set seed
#' set.seed(19723)
#' # Create data set and analyze
#' # Note: Generally Use more than 50 replications
#' model |> power_analysis(50, 5, 50)
#' @export
power_analysis <- function(
        model,
        replications,
        n_within,
        n_between,
        ...) {

    # Validate inputs

    # Validate model
    is_valid(model)

    if (!is.numeric(n_within)) throw_error(
        '{.arg n_within} must be a numeric vector'
    )
    if (!is.numeric(n_between)) throw_error(
        '{.arg n_between} must be a numeric vector'
    )
    if (replications < 2) throw_error(
        '{.arg replications} must be greater than 1.'
    )

    # Get icc
    icc <- model$effect_size$icc

    # Check for timevars and adjust n_within accordingly
    timevars <- which(vapply(
        model$predictors,
        \(.)  'mp_timevar' %in% class(.),  # Select timevar
        logical(1L)
    ))

    if (missing(n_within)) {
        if (length(timevars) == 0) throw_error(
            '{.arg n_within} is missing with no time variable.'
        )
        n_within <- length(model$predictors[[timevars]]$values)
    } else if (length(timevars) == 1) {
        len <- length(model$predictors[[timevars]]$values)
        # Loop over n_within
        alert <- TRUE # Only alert once
        for (i in seq_along(n_within)) {
            if (len < n_within[i]) {
                if (alert) {
                    cli::cli_alert_info(
                        "Setting n_within = {len} to match time variable's length"
                    )
                    alert <- FALSE
                }
                n_within[i] <- len
            }
        }
    }
    else if (length(timevars) > 1) {
        throw_error('Only one time variable is currently allowed.')
    }

    # Run simulation based on all combos
    results <- cmap(
        \(icc, n_within, n_between) {
            cli::cli_alert(c(
                "ICC = {formatC(icc, digits = 2, format = 'f' )}, ",
                "n_within = {n_within}, ",
                "n_between = {n_between}"
            ))
            model |> subset(icc = icc) |> `_power_analysis`(
                replications,
                n_within,
                n_between,
                ...
            )
        },
        icc = icc,
        n_within = n_within,
        n_between = n_between
    )

    # Return object if only length of 1
    if (length(results) == 1) return(results[[1]])

    # Otherwise reconstruct into single object

    # Get named list
    nam <-  apply(
        expand.grid(icc = icc, n_within = n_within, n_between = n_between),
        MARGIN = 1,
        FUN = \(.)  paste(names(.), '=' ,., collapse = ', ')
    )

    # Return results
    return(structure(list(
        sim = list(
            model = model,
            replications = replications,
            n_within = n_within,
            n_between = n_between
        ),
        power           = lapply(setNames(results, nam), \(.) .$power),
        estimates       = lapply(setNames(results, nam), \(.) .$estimates),
        mean_parameters = lapply(setNames(results, nam), \(.) .$mean_parameters)
    ), class = c('mp_power', 'mp_base')))

}

#' Prints a [`mlmpower::mp_power`]
#' @description
#' Prints a [`mlmpower::mp_power`] in a human readable format.
#' @param x a [`mlmpower::mp_power`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
#' # Set seed
#' set.seed(19723)
#' # Create data set and analyze
#' # Note: Generally Use more than 50 replications
#' model |> power_analysis(50, 5, 50) -> powersim
#' # Print results
#' print(powersim)
#' @export
print.mp_power <- function(x, ...) {
    cli::cli_h2('Power Analysis Specification')
    cli::cli_text('Replications = {x$sim$replications}')
    cli::cli_text('N_Within = {x$sim$n_within}')
    cli::cli_text('N_Between = {x$sim$n_between}')
    cli::cli_h3('lme4 Model Syntax')
    cli::cli_text('')
    # Handle only one icc case
    if (is.null(names(x$power))) {
        cli::cli_text(deparse(to_formula(x$mean_parameters)))
    } else {
        cli::cli_text(deparse(to_formula(x$mean_parameters[[1]])))
    }
    cli::cli_text('')
    cli::cli_alert_info('cwc() = centering within cluster')
    cli::cli_alert_info('cgm() = centering with grand mean')
    cli::cli_h3('Effect Sizes')
    cli::cli_text('')
    cli::cli_ul()
    cli::cli_li('WITHIN = {x$sim$model$effect_size$within}')
    cli::cli_li('BETWEEN = {x$sim$model$effect_size$between}')
    cli::cli_li('RANDOM SLOPE = {x$sim$model$effect_size$random_slope}')
    cli::cli_li('PRODUCT = {x$sim$model$effect_size$product}')
    cli::cli_end()

    cli::cli_h2('Power Analysis Results')
    # Handle only one icc case
    if (is.null(names(x$power))) {
        cli::cli_h3('global icc = {x$sim$model$effect_size$icc}')
        power_table(x$power)
    } else {
        for (i in names(x$power)) {
            cli::cli_h3('global {.cli {i}}')
            power_table(x$power[[i]])
        }
    }
    cli::cli_text('')
    cli::cli_alert_info('Margin of error (MOE) computed as \u00B1 1.96 * standard error')
    invisible(x)
}

#' Summarizes a [`mlmpower::mp_power`]
#' @description
#' Summarizes a [`mlmpower::mp_power`] in a human readable format.
#' This is a simple wrapper for [`mlmpower::print.mp_power`].
#' @param object a [`mlmpower::mp_power`].
#' @param ... other arguments not used by this method.
#' @returns Invisibly returns the original variable.
#' @examples
#' # Create Model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(icc = 0.1)
#' )
#' # Set seed
#' set.seed(19723)
#' # Create data set and analyze
#' # Note: Generally Use more than 50 replications
#' model |> power_analysis(50, 5, 50) -> powersim
#' # Summarizes results
#' summary(powersim)
#' @export
summary.mp_power <- function(object, ...) {
    print(object)
}
