#' Generate a data set based on `mp_model`
#'
#' @export
simulate.mp_model <- function(object, n_within, n_between, nsim = 1) {

    # TODO validate inputs

    # Obtain levels
    lvls <- sapply(object$predictors, levels)

    # Obtain timevar level-1 indicators
    timevar_l1 <- sapply(
        object$predictors[lvls == 1],  # Subset level-1
        \(.)  'mp_timevar' %in% class(.)  # Select timevar
    )
    if (sum(timevar_l1) > 1) cli::cli_abort("Only one timevar should be specified")

    # Change within sample size if timevar specified
    if (TRUE %in% timevar_l1) {
        len <- length(object$predictors[lvls == 1][[which(timevar_l1)]]$values)
        if (missing(n_within)) {
            n_within <- len
        } else if (len != n_within) {
            cli::cli_alert_info(
                "Setting n_within = {len} to match time variable's length"
            )
            n_within <- len
        }
    }

    # Return list if nsim isn't 1
    if (nsim > 1) {
        return(
            lapply(seq_len(nsim), \(.) {
                object |> simulate(n_within, n_between)
            })
        )
    }

    # Obtain binary level-2 indicators
    binary_l2 <- sapply(
        # Subset level-2
        object$predictors[lvls == 2],  # Subset level-2
        \(.)  'mp_binary' %in% class(.)   # Select binary
    )

    # Create parameters
    p <- new_parameters(object)

    # useful precomputed values
    N <- n_within * n_between # Total sample size
    l1 <- length(p$mean_X)
    l2 <- length(p$mean_W)

    # Generate ID variable
    `_id` <- seq_len(n_between) %x% matrix(1, n_within)


    # Generate X predictors

    # TODO deal with timevars
    # Generate Level-1 Time variables if exists
    if (TRUE %in% timevar_l1) {
        # Get values
        vals <- object$predictors[lvls == 1][[which(timevar_l1)]]$values

        # Create X_w matrix
        X_w <- matrix(0, N, l1)
        X_w[,  timevar_l1] <-  matrix(vals, nrow = N, ncol = 1)

        # Generate X's conditional on timevars
        if (l1 > 1) {
            cond_var <- (
                p$phi_XX_w[ !timevar_l1, !timevar_l1, drop = F]
                - p$phi_XX_w[!timevar_l1, timevar_l1, drop = F]
                %*% solve(
                    p$phi_XX_w[timevar_l1, timevar_l1, drop = F],
                    p$phi_XX_w[timevar_l1, !timevar_l1, drop = F]
                )
            )
            cond_mean <- t(
                p$phi_XX_w[!timevar_l1, timevar_l1, drop = F]
                %*% solve(
                    p$phi_XX_w[timevar_l1, timevar_l1, drop = F],
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
        X_w <- rmvnorm_nomean(N, p$phi_XX_w)
        X_b <- rmvnorm(N, c(p$mean_X, rep(0, l2)), p$phi_b)[`_id`, , drop = F]
    }

    # Generate Binary Between Predictors
    if (TRUE %in% binary_l2) {
        bin_select <- seq_len(l2)[binary_l2] + l1
        X_binary <- X_b[ , bin_select, drop = F]
        # Obtain thresholds
        thresh <- qnorm(p$mean_W[binary_l2], sd = sqrt(diag(p$phi_b)[bin_select]))
        # Generate data
        for (i in seq_along(thresh)) {
            X_b[ , bin_select[i]] <- ifelse(
                X_b[ , bin_select[i]] < thresh[i],
                0.0, 1.0
            ) - p$mean_W[binary_l2] # Subtract mean for centering
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
         # Add means for W
         + matrix(
             ifelse(binary_l2, 0, p$mean_W),
             nrow = N,
             ncol = l2,
             byrow = T
         )
        )
    )
    names(d) <- c(
        '_id',
        names(p$mean_Y),
        names(which(sapply(object$predictors, levels) == 1)),
        names(which(sapply(object$predictors, levels) == 2))
    )
    # Add timevar_l1 to parameters
    p$timevar_l1 <- timevar_l1
    attr(d, 'parameters') <- p
    return(d)
}



#' Internal function to Analyze one replication for `mp_model`
#'
#' @noRd
analyze <- function(model, n_within, n_between, alpha = 0.05) {

    # Obtain data
    model |> simulate(n_within, n_between) -> d

    # Make centering env
    e <- centering_env(d$`_id`)

    # Get formula for model
    attr(d, 'parameters') |> to_formula(e) -> f

    # Fit Model with lme4 and return results
    full_reml <- quiet(lmerTest::lmer(f, d))
    full_ml   <- quiet(lme4::lmer(f, d, REML = FALSE))
    nested    <- quiet(lme4::lmer(attr(d, 'parameters') |> to_formula(e, nested = T), d, REML = FALSE))

    # Obtain LRT
    lrt <- varTestnlme::varCompTest(full_ml, nested, pval.comp = "bounds", output = FALSE)

    # Return output
    list(
        estimates = full_reml |> extract_results(),
        sig_test = list(
            fixed  = coefficients(summary(full_reml) )[,'Pr(>|t|)'] < alpha,
            random = c(omnibus_test = lrt$p.value[[4]] < alpha)
        ),
        parameters =  attr(d, 'parameters')
    )
}



#' Conduct a power analysis based on `mp_model`
#'
#' @export
power_analysis <- function(
        model,
        replications,
        n_within,
        n_between) {


    # Get icc
    icc <- model$effect_size$icc

    # Check for missing n_within
    if (missing(n_within)) {
        timevars <- which(sapply(
            model$predictors,
            \(.)  'mp_timevar' %in% class(.)  # Select timevar
        ))
        if (length(timevars) == 0) cli::cli_abort(
            '{.cli n_within} is missing with no time variable.'
        )
        len <- length(model$predictors[[timevars]]$values)
        if (missing(n_within)) {
            n_within <- len
        } else if (len != n_within) {
            cli::cli_alert_info(
                "Setting n_within = {len} to match time variable's length"
            )
            n_within <- len
        }
    }

    # TODO validate inputs
    # TODO dont' allow reps to be 2

    # Handle multiple ICCs
    if (length(icc) > 1) {

        # Create names
        names(icc) <- paste0('icc = ', icc)

        # Run simulation
        results <- lapply(icc, \(x) {
            model |> subset(icc = x) |> power_analysis(
                replications,
                n_within,
                n_between
            )
        })

        # Return results
        return(structure(list(
            sim = list(
                model = model,
                replications = replications,
                n_within = n_within,
                n_between = n_between
            ),
            power = lapply(results, \(.) .$power),
            estimates = lapply(results, \(.) .$estimates),
            mean_parameters = lapply(results, \(.) .$mean_parameters)
        ), class = c('mp_power', 'mp_base')))
    }
    # Otherwise handle single replication

    # Run reps and include progress bar
    results <- lapply(
        cli::cli_progress_along(
            seq_len(replications),
            clear = FALSE,
            format = paste0(
                "   Simulating ICC = ", icc ,
                ": {cli::pb_bar} {cli::pb_current}/{cli::pb_total} ",
                "| Elapsed: {cli::pb_elapsed}",
                "| ETA: {cli::pb_eta}"
            )
        ), \(.) {
            model |> analyze(n_within, n_between)
        }
    )

    # Extract results
    e <- sapply(results, \(.) unlist(.$estimates))
    s <- sapply(results, \(.) unlist(.$sig_test))
    p <- Reduce('+', lapply(results, \(.) .$parameters))

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
            mean_parameters = p
        ),
        class = c('mp_power', 'mp_base')
    )
}

#' Print a `mp_power`
#'
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
    cli::cli_alert_info('Margin of Error (MOE) Computed as ±1.96 * standard error')
    invisible(x)
}

#' Summary for `mp_power`
#'
#' @export
summary.mp_power <- function(object, ...) {
    object
}
