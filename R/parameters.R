#' @rdname Parameters
#' @title `mlmpower` Modeling Framework
#' @name Parameters
#' @aliases mp_parameters parameters
#' @description Test
# TODO finish this documentation
NULL

#' Internal function to convert `mp_model` object to an environment
#' @noRd
to_env <- function(model) {

    # Compute with model
    l <- with(model, {

        # Obtain level-1 and level-2 vars
        l1 <- which(vapply(predictors, levels, numeric(1L)) == 1)
        l2 <- which(vapply(predictors, levels, numeric(1L)) == 2)

        # Selectors for products
        sel_p <- vapply(actions, \(.) .$type == "product", logical(1L))
        # Selectors for random slopes
        sel_r <- vapply(actions, \(.) .$type == "random_slope", logical(1L))

        # Make product weight matrix
        # NOTE Assumes first name is always level-1
        pmat <- matrix(0, length(l2), length(l1))
        for (a in actions[sel_p]) {
            j <- which(names(l1) %in% a$name[1])
            i <- which(names(l2) %in% a$name[2])
            pmat[i, j] <- a$weight
        }
        # Make random slope weights
        rvec <- vector("numeric", length(l1))
        for (a in actions[sel_r]) {
            i <- which(names(l1) %in% a$name)
            rvec[i] <- a$weight
        }

        # Create list
        within(list(), {
            # Get means
            mean_Y <- outcome$mean
            names(mean_Y) <- outcome$name
            mean_X <- vapply(predictors[l1], \(.) .$mean, numeric(1L))
            mean_Z <- vapply(predictors[l2], \(.) .$mean, numeric(1L))

            # Get model implied means of X based on centering
            model_mean_X <- ifelse(
                vapply(predictors[l1], \(.) 'mp_timevar' %in% class(.), logical(1L)),
                mean_X,
                0
            )
            # Get variances
            var_Y <- outcome$sd^2
            var_X <- vapply(predictors[l1], \(.) .$sd^2, numeric(1L))
            var_Z <- vapply(predictors[l2], \(.) .$sd^2, numeric(1L))
            # Get correlations
            corr_X <- corrs$within
            corr_Z <- corrs$between
            corr_X_Z <- corrs$between
            corr_raneffects <- corrs$randeff
            # Get icc (Assumes only one)
            icc_Y <- if (is.null(outcome$icc)) effect_size$icc else outcome$icc
            icc_X <- vapply(predictors[l1], \(.) {
                if (is.null(.$icc)) effect_size$icc else .$icc
            }, numeric(1L))
            # Get R2s
            R2_X_w <- effect_size$within
            R2_increment_b <- effect_size$between
            R2_XXproduct_w <- NULL # Place holders (not used currently)
            R2_ZZproduct_b <- NULL # Place holders (not used currently)
            R2_XZproduct_w <- effect_size$product
            R2_ranslopes_w <- effect_size$random_slope

            # Get weights
            weights_X_w <- vapply(predictors[l1], \(.) .$weight, numeric(1L))
            weights_XXproduct_w <- NULL # Place holders (not used currently)
            weights_ZZproduct_b <- NULL # Place holders (not used currently)
            weights_XZproduct_w <- c(pmat)
            weights_ranslopes_w <- rvec
            weights_increment_b <- c(
                vapply(icc_X, \(.) if (. == 0) 0 else NA, numeric(1L)),
                vapply(predictors[l2], \(.) .$weight, numeric(1L))
            )

            # Set binary Z for level2
            binary_Z <- vapply(predictors[l2], \(.) {
                if ("mp_binary" %in% class(.)) .$mean else 0.0
            }, numeric(1L))
        })
    })

    # Return list as env
    return(list2env(l))
}

#' Internal function to solve parameters based on
#' a converted `mp_model` object.
#' @noRd
make_parameters <- function(model) {
    # Convert model to old specification
    model |> to_env() -> env

    # Create list using env
    l <- with(env, {
        # variable counts
        num_X <- length(mean_X)
        num_Z <- length(var_Z)

        # set binary variable variance
        var_Z[binary_Z > 0 & binary_Z < 1] <- binary_Z[binary_Z > 0 & binary_Z < 1] * (1 - binary_Z[binary_Z > 0 & binary_Z < 1])

        # compute within and between-cluster variances
        icc_X[weights_increment_b[seq_len(num_X)] == 0] <- 0
        var_X_b <- var_X * icc_X
        var_X_w <- var_X - var_X_b
        var_Y_b <- var_Y * icc_Y
        var_Y_w <- var_Y - var_Y_b

        # construct the within-cluster correlation and covariance matrix of the level-1 Xs
        cor_XX_w <- matrix(1, nrow = num_X, ncol = num_X)
        cor_XX_w[lower.tri(cor_XX_w)] <- cor_XX_w[upper.tri(cor_XX_w)] <- corr_X(num_X * (num_X - 1) / 2)
        phi_XX_w <- diagonal(sqrt(var_X_w)) %*% cor_XX_w %*% diagonal(sqrt(var_X_w))

        # construct the between-cluster covariance matrix of X's level-2 group means
        phi_XX_b <- diagonal(sqrt(var_X_b)) %*% cor_XX_w %*% diagonal(sqrt(var_X_b))

        # construct the between-cluster correlation and covariance matrix of the level-2 Zs
        cor_ZZ_b <- matrix(1, nrow = num_Z, ncol = num_Z)
        cor_ZZ_b[lower.tri(cor_ZZ_b)] <- cor_ZZ_b[upper.tri(cor_ZZ_b)] <- corr_Z(num_Z * (num_Z - 1) / 2)
        phi_ZZ_b <- diagonal(sqrt(var_Z)) %*% cor_ZZ_b %*% diagonal(sqrt(var_Z))

        # for binary level-2 variables, convert inputted correlations to point-biserial correlations then to covariances
        if (is.null(binary_Z)) {
            binary_Z_ind <- rep(FALSE, length(var_Z))
        } else {
            binary_Z_ind <- binary_Z > 0 & binary_Z < 1
        }
        for (r in seq_len(num_Z)) {
            for (i in seq_len(num_Z)) {
                if (binary_Z_ind[r] == T & i != r) {
                    cor_pointbi <- cor_ZZ_b[r, i] / sqrt(binary_Z[r] * (1 - binary_Z[r])) * dnorm(qnorm(binary_Z[r]))
                    cor_ZZ_b[r, i] <- cor_ZZ_b[i, r] <- cor_pointbi
                    phi_ZZ_b[r, i] <- phi_ZZ_b[i, r] <- cor_ZZ_b[r, i] * sqrt(phi_ZZ_b[r, r] * phi_ZZ_b[i, i])
                }
            }
        }

        # construct the between-cluster covariance matrix of the X cluster means and level-2 Zs
        cor_XZ_b <- matrix(corr_X_Z(num_X * num_Z), nrow = num_Z, ncol = num_X)
        phi_XZ_b <- diagonal(sqrt(var_Z)) %*% cor_XZ_b %*% diagonal(sqrt(var_X_b))

        # construct the between-cluster covariance matrix
        phi_b <- rbind(cbind(phi_XX_b, t(phi_XZ_b)), cbind(phi_XZ_b, t(phi_ZZ_b)))

        # construct the within-cluster covariance matrix
        phi_XZwithXZ_w <- phi_XX_w %x% phi_ZZ_b
        phi_XwithXZ_w <- matrix(0, nrow = NROW(phi_XZwithXZ_w), ncol = NCOL(phi_XX_w))
        phi_w <- rbind(cbind(phi_XX_w, t(phi_XwithXZ_w)), cbind(phi_XwithXZ_w, t(phi_XZwithXZ_w)))

        # solve for the within-cluster regression coefficients
        weights_scaled <- 1 / sqrt(diagonal(phi_XX_w)) * weights_X_w
        gamma_X_w <- weights_scaled * c(sqrt((var_Y * R2_X_w) / t(weights_scaled) %*% phi_XX_w %*% weights_scaled))

        # solve for the cross-level product coefficients
        weights_scaled <- 1 / sqrt(diagonal(phi_XZwithXZ_w)) * weights_XZproduct_w
        if (sum(weights_XZproduct_w) == 0) {
            gamma_XZ_w <- weights_XZproduct_w
        } else {
            gamma_XZ_w <- weights_scaled * c(sqrt((var_Y * R2_XZproduct_w) / t(weights_scaled) %*% phi_XZwithXZ_w %*% weights_scaled))
        }
        gamma_w <- c(gamma_X_w, gamma_XZ_w)

        # compute the within-cluster residual variance
        var_e_w <- (1 - icc_Y - R2_X_w - R2_XZproduct_w - R2_ranslopes_w) * var_Y

        # solve for the between-cluster coefficients
        select_weighted <- !is.na(weights_increment_b)
        select_nonweighted <- is.na(weights_increment_b)
        if (sum(select_weighted) != NROW(phi_b)) {
            phi_nonweighted_b <- phi_b[select_nonweighted, select_nonweighted, drop = F]
            phi_weighted_b <- phi_b[select_weighted, select_weighted, drop = F]
            phi_covs_b <- phi_b[select_nonweighted, select_weighted, drop = F]
            resvar_Z_b <- phi_weighted_b - t(solve(phi_nonweighted_b) %*% phi_covs_b) %*% phi_covs_b

            # Predefine and only compute for non zero phi_b
            weights_scaled <- numeric(sum(select_weighted))
            sel_weight <- diagonal(resvar_Z_b) != 0 # Select out only non-zero diagonals
            weights_scaled[sel_weight] <- 1 / sqrt(diagonal(resvar_Z_b)[sel_weight]) * weights_increment_b[select_weighted][sel_weight]
            gamma_weighted_b <- weights_scaled * c(sqrt((var_Y * R2_increment_b) / t(weights_scaled) %*% resvar_Z_b %*% weights_scaled))
            gamma_b <- c(gamma_w[seq_len(num_X)], rep(0, num_Z))
            gamma_b[select_weighted] <- gamma_weighted_b
            gamma_nonweighted_b <- gamma_b[select_nonweighted]
        } else {
            # Predefine and only compute for non zero phi_b
            weights_scaled <- numeric(NROW(phi_b))
            sel_weight <- diagonal(phi_b) != 0 # Select out only non-zero diagonals
            weights_scaled[sel_weight] <- 1 / sqrt(diagonal(phi_b)[sel_weight]) * weights_increment_b[sel_weight]
            gamma_b <- weights_scaled * c(sqrt((var_Y * R2_increment_b) / t(weights_scaled) %*% phi_b %*% weights_scaled))
        }
        # Replace NaN with 0
        gamma_b[is.na(gamma_b)] <- 0

        # compute the random effect correlation matrix
        cor_raneffects <- diag(nrow = num_X + 1, ncol = num_X + 1)
        cor_sel <- which(weights_ranslopes_w != 0)
        cor_raneffects[1, cor_sel + 1] <- cor_raneffects[cor_sel + 1, 1] <- corr_raneffects(length(cor_sel))

        # solve for the random slope variances
        cor_ranslopes <- cor_raneffects[-1, -1, drop = F]
        tau_trace <- (var_Y * R2_ranslopes_w) / sum(diagonal(cor_ranslopes %*% phi_XX_w %*% diagonal(weights_ranslopes_w / diagonal(phi_XX_w))))
        var_ranslopes <- if (is.finite(tau_trace)) weights_ranslopes_w * tau_trace / diagonal(phi_XX_w) else rep(0.0, NROW(weights_ranslopes_w))

        # vector of correlations between random intercept and slopes
        cor_is <- cor_raneffects[-1, 1, drop = F]

        # covariance matrix of just the random slopes
        tau_ranslopes <- diagonal(sqrt(c(var_ranslopes))) %*% cor_raneffects[-1, -1, drop = F] %*% diagonal(sqrt(c(var_ranslopes)))

        # compute random intercept variance
        # explained level-2 variation
        b <- t(gamma_b) %*% phi_b %*% gamma_b
        # random intercept variation due to non-zero level-1 means
        a <- t(model_mean_X) %*% (cor_is * sqrt(var_ranslopes))
        s <- t(model_mean_X) %*% tau_ranslopes %*% model_mean_X

        # Precompute portion of tau00
        comp1 <- -4 * b * a^2 + a^4 - 4 * a^2 * s + 4 * a^2 * var_Y_b
        # Check if comp1 is positive
        if (is.na(comp1) | comp1 < 0) {
            throw_error(c(
                'The random intercept variance is negative.',
                'x' = 'This is caused because the effect size specified are impossible.',
                'i' = 'The between effect size is most likely too large compared to the ICC.',
                '>' = 'ICC = {icc_Y} and Between R2 = {R2_increment_b}'
            ))
        }
        tau00 <- 0.5 * (-2 * b + a^2 - 2 * s + 2 * var_Y_b) - 0.5 * sqrt(comp1)
        # Check if tau00 is positive
        if (is.na(tau00) | tau00 < 0) {
            throw_error(c(
                'The random intercept variance is negative.',
                'x' = 'This is caused because the effect size specified are impossible.',
                'i' = 'The between effect size is most likely too large compared to the ICC.',
                '>' = 'ICC = {icc_Y} and Between R2 = {R2_increment_b}'
            ))
        }

        # compute intercept-slope covariance and construct tau matrix
        tau <- diagonal(sqrt(c(tau00, var_ranslopes))) %*% cor_raneffects %*% diagonal(sqrt(c(tau00, var_ranslopes)))
        cor_raneffects[tau == 0] <- 0

        # compute fixed intercept and construct coefficient matrix
        # the mean of the product from Bohrnstedt & Goldberger Equation 3 simplifies because cov(X_w,Z) = 0
        means <- c(model_mean_X, rep(0, num_X + num_Z + num_X*num_Z)) # Set all Z means to 0
        gamma00 <- mean_Y - c(gamma_w, gamma_b) %*% means
        gammas <- c(gamma00, gamma_w, gamma_b)

        # variable names
        if (length(names(mean_X)) != 0) {
            vars_Xw <- paste0(names(mean_X), "_w")
            vars_Xb <- paste0(names(mean_X), "_b")
        } else {
            vars_Xw <- vars_Xb <- names(mean_X)
        }
        vars_Z <- names(mean_Z)
        vars_Z[binary_Z_ind == T] <- paste0(vars_Z[binary_Z_ind == T], "_binary")
        vars_XZ <- unlist(lapply(vars_Xw, \(x) {
            vapply(vars_Z, \(w) {
                paste0(x, "*", w)
            }, character(1L))
        }))

        # collect parameters and construct names
        params_coeff <- matrix(c(gammas), ncol = 1)
        params_res <- matrix(c(var_e_w), ncol = 1)
        row.names(tau) <- colnames(tau) <- c("Icept", vars_Xw)
        row.names(params_coeff) <- c("Icept", vars_Xw, vars_XZ, vars_Xb, vars_Z)
        row.names(params_res) <- c("Res. Var.")
        colnames(params_coeff) <- colnames(params_res) <- "Value"
        row.names(phi_XX_w) <- colnames(phi_XX_w) <- vars_Xw
        row.names(phi_b) <- colnames(phi_b) <- c(vars_Xb, vars_Z)
        row.names(phi_XZwithXZ_w) <- colnames(phi_XZwithXZ_w) <- vars_XZ

        # R-square summary
        check_var_Y <- (
            t(gamma_w) %*% phi_w %*% gamma_w + t(gamma_b) %*%
                phi_b %*% gamma_b + sum(diagonal(tau[-1, -1, drop = F] %*% phi_XX_w))
            + t(model_mean_X) %*% tau_ranslopes %*% model_mean_X + t(model_mean_X) %*%
                tau[-1, 1, drop = F] + tau00 + var_e_w
        )
        # check_var_Y_w <- (
        #     t(gamma_w) %*% phi_w %*% gamma_w
        #     + sum(diag(tau[-1, -1, drop = F] %*% phi_XX_w)) + var_e_w
        # )
        # check_var_Y_b <- (
        #     t(gamma_b) %*% phi_b %*% gamma_b + t(model_mean_X) %*%
        #         tau_ranslopes %*% model_mean_X + t(model_mean_X) %*% tau[-1,1]
        #     + tau00
        # )
        R2check_X_w <- t(gamma_X_w) %*% phi_XX_w %*% gamma_X_w / check_var_Y
        R2check_XZ_w <- t(gamma_XZ_w) %*% phi_XZwithXZ_w %*% gamma_XZ_w / check_var_Y
        R2check_ranslopes_w <- sum(diagonal(tau_ranslopes %*% phi_XX_w)) / check_var_Y
        R2check_var_e <- ((1 - icc_Y) * var_Y - t(gamma_w) %*% phi_w %*% gamma_w - sum(diagonal(tau_ranslopes %*% phi_XX_w))) / check_var_Y
        R2check_XZ_b <- t(gamma_b) %*% phi_b %*% gamma_b / check_var_Y
        if (sum(select_weighted) != NROW(phi_b)) {
            R2check_increment_b <- t(gamma_weighted_b) %*% resvar_Z_b %*% gamma_weighted_b / check_var_Y
        } else {
            R2check_increment_b <- R2check_XZ_b
        }
        R2check_totalminusincrement_b <- R2check_XZ_b - R2check_increment_b
        R2check_ranicept <- tau00 / check_var_Y

        # Collect r2 summaries
        r2  <- matrix(
            c(
                R2check_X_w,
                R2check_XZ_w,
                R2check_ranslopes_w,
                R2check_var_e,
                R2check_XZ_b,
                R2check_increment_b,
                R2check_totalminusincrement_b,
                R2check_ranicept
            ),
            dimnames = list(c(
                "Variance Within-Cluster Fixed Effects",
                "Variance Cross-Level Interaction Effects",
                "Variance Random Slopes",
                "Within-Cluster Error Variance",
                "Variance Between-Cluster Fixed Effects",
                "Incremental Variance Level-2 Predictors",
                "Between Variance Level-1 Covariates",
                "Variance Random Intercepts"
            ), 'Proportion')
        )

        # Subset phi_XZwithXZ_w based on actual requested products
        phi_XZwithXZ_w <- phi_XZwithXZ_w[gamma_XZ_w != 0, gamma_XZ_w != 0, drop = F]

        # Return final output
        list(
            mean_Y = mean_Y,
            gammas = params_coeff,
            tau = tau,
            var_e = params_res,
            mean_X = mean_X,
            mean_Z = mean_Z,
            phi_w = phi_XX_w,
            phi_p = phi_XZwithXZ_w,
            phi_b = phi_b,
            r2 = r2
        )
    })

    # Return mp_parameters class
    structure(list2env(l), class = c("mp_parameters", "mp_base"))
}

#' Internal function to solve parameters based on
#' a converted `mp_model` object across multiple reps
#' @noRd
make_avg_parameters <- function(model) {
    # Check if fixed and return normally
    if (is_fixed_cor(model$corrs)) return(model |> make_parameters())

    # Otherwise find average correlation
    # Create temp model with average corrs
    new_model <- (
        clone(model)
        + correlations(
            within  = fixed(mean(model$corrs$within)),
            between = fixed(mean(model$corrs$between)),
            randeff = fixed(mean(model$corrs$randeff))
        )
    )

    # output parameters
    new_model |> make_parameters()
}

#' Internal function to clean parameters print outs
#' @noRd
clean_parameters <- function(x) {

    # Keep copy of original
    attr(x, '_gammas') <- x$gammas
    attr(x, '_phi_b')  <- x$phi_b
    attr(x, '_tau')    <- x$tau

    # Drop 0 gammas
    x$gammas <- x$gammas[x$gammas != 0, , drop = F]
    x$phi_b <- x$phi_b[diag(x$phi_b) != 0, diag(x$phi_b) != 0, drop = F]
    x$tau <- x$tau[diag(x$tau) != 0, diag(x$tau) != 0, drop = F]

    # Return
    x
}

#' Validate parameters
#' @noRd
is.parameters <- function(x) {
    inherits(x, "mp_parameters")
}

#' Convert `mp_parameters` object to a formula for `lme4`
#' @noRd
to_formula <- function(x, e = globalenv(), nested = FALSE) {

    # Get regression coefficients
    gammas <- if (is.null(attr(x, '_gammas'))) x$gammas else attr(x, '_gammas')
    tau    <- if (is.null(attr(x, '_tau')))    x$tau    else attr(x, '_tau')

    # Obtain number of l1 and l2
    n_l1 <- NROW(x$mean_X)
    n_l2 <- NROW(x$mean_Z)

    # Determine which predictors are CGM
    # NOTE Assumes always same regression coefficient means cgm
    cgm_sel <- gammas[seq_len(n_l1) + 1] == gammas[seq_len(n_l1) + (1 + n_l1 + n_l1 * n_l2)]

    # Variable names
    var_l1 <- ifelse(
        cgm_sel,
        paste0("cgm(", names(x$mean_X), ")"),
        paste0("cwc(", names(x$mean_X), ")")
    )

    # Check for timevar
    if (!is.null(x$timevar_l1)) {
        var_l1[x$timevar_l1] <- names(x$mean_X[x$timevar_l1])
    }

    var_l2 <- if (length(names(x$mean_Z)) == 0) {
        character()
    } else {
        paste0("cgm(", names(x$mean_Z), ")")
    }

    # Create products
    var_prod <- unlist(lapply(var_l1, \(x) {
        vapply(var_l2, \(w) {
            paste0(x, ":", w)
        }, character(1L))
    }))

    # Set between to 0 so its dropped
    gammas[seq_len(n_l1) + (1 + n_l1 + n_l1 * n_l2)][cgm_sel] <- 0

    # Create Formula
    model_fixed <- paste0(c(
        # Intercept
        "1",
        # Level-1 predictors
        var_l1,
        # Interactions
        var_prod[gammas[seq_len(n_l1 * n_l2) + 1 + n_l1] != 0],
        # level-2
        var_l2
    ), collapse = " + ")

    if (nested) {
        model_random <- "(1 | `_id`)"
    } else {
        raneff_model <- c("1", var_l1[diagonal(tau)[-1] != 0])
        raneff_model <- paste0(raneff_model, collapse = " + ")
        model_random <- paste0("(", raneff_model, " | ", "`_id`)")
    }
    lme4_model <- paste0(names(x$mean_Y), " ~ ", model_fixed, " + ", model_random)
    return(as.formula(lme4_model, e))
}

#' Convert [`mlmpower::mp_parameters`] to a [`list`]
#' @description
#' A wrapper to coerce [`mlmpower::mp_parameters`] to a [`list`].
#' @param x the [`mlmpower::mp_parameters`] to be coered.
#' @param ... additional arguments passed to [`as.list`]
#' @returns a [`list`]
#' @examples
#' # Specify model
#' model <- (
#'     outcome('Y')
#'     + within_predictor('X')
#'     + effect_size(
#'         icc = c(0.1, 0.2),
#'         within = 0.3
#'     )
#' )
#' # Obtain parameters and convert to a list
#' model |> summary() |> as.list() -> param_list
#' @export
as.list.mp_parameters <- function(x, ...) {
    as.list.environment(x, ...)
}

#' Prints a [`mlmpower::mp_parameters`]
#' @description
#' Prints a [`mlmpower::mp_parameters`] in a human readable format.
#' @param x a [`mlmpower::mp_parameters`].
#' @param ... arguments passed to [`print()`].
#' @returns Invisibly returns the original variable.
#' @examples
#' print(
#'     summary(
#'         outcome('Y')
#'         + within_predictor('X')
#'         + effect_size(icc = cross_sectional)
#'     )
#' )
#' @export
print.mp_parameters <- function(x, ...) {
    l <- lapply(as.list(x), function(x) {
        attr(x, '_count') <- NULL
        x
    })
    print(l[names(l) != 'timevar_l1'], ...)
    invisible(x)
}

#' Internal function to keep online average for `mp_parameters`
#' @noRd
average <- function(e1, e2) {
    # Check if one is missing
    if (missing(e1)) return(e2)
    if (missing(e2)) return(e1)

    # Don't allow two models
    if (!is.parameters(e1) | !is.parameters(e2)) {
        throw_error('Can only average two parameters together')
    }

    # Iterate over everything in e1 and e2
    for (i in ls(envir = e1)) {
        if (i == 'timevar_l1') next  # skip timevar
        assign(
            i,
            online_mean(get(i, e1), get(i, e2)),
            e1
        )
    }
    return(e1)
}

