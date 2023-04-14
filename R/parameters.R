#' Internal function to convert `mp_model` object to old set up
#'
#' @noRd
`_convert` <- function() {
    # obtain level-1 and level-2 vars
    l1 <- which(vapply(predictors, levels, numeric(1L)) == 1)
    l2 <- which(vapply(predictors, levels, numeric(1L)) == 2)

    # Selectors for products
    sel_p <- vapply(actions, \(.) .$type == "product", logical(1L))
    # selectors for random slopes
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
    l <- within(list(), {
        # Get means
        mean_Y <- outcome$mean
        names(mean_Y) <- outcome$name
        mean_X <- vapply(predictors[l1], \(.) .$mean, numeric(1L))
        mean_W <- vapply(predictors[l2], \(.) .$mean, numeric(1L))

        # Get model implied means of X based on centering
        model_mean_X <- ifelse(
            vapply(predictors[l1], \(.) 'mp_timevar' %in% class(.), logical(1L)),
            mean_X,
            0
        )
        # Get variances
        var_Y <- outcome$sd^2
        var_X <- vapply(predictors[l1], \(.) .$sd^2, numeric(1L))
        var_W <- vapply(predictors[l2], \(.) .$sd^2, numeric(1L))
        # Get correlations
        corr_X <- corrs$within_cor
        corr_W <- corrs$between_cor
        corr_X_W <- corrs$between_cor
        corr_raneffects <- corrs$randeff_cor
        # Get icc (Assumes only one)
        icc_Y <- if (is.null(outcome$icc)) effect_size$icc else outcome$icc
        icc_X <- vapply(predictors[l1], \(.) {
            if (is.null(.$icc)) effect_size$icc else .$icc
        }, numeric(1L))
        # Get R2s
        R2_X_w <- effect_size$within
        R2_increment_b <- effect_size$between
        R2_XXproduct_w <- NULL # Place holders (not used currently)
        R2_WWproduct_b <- NULL # Place holders (not used currently)
        R2_XWproduct_w <- effect_size$product
        R2_ranslopes_w <- effect_size$random_slope

        # Get weights
        weights_X_w <- vapply(predictors[l1], \(.) .$weight, numeric(1L))
        weights_XXproduct_w <- NULL # Place holders (not used currently)
        weights_WWproduct_b <- NULL # Place holders (not used currently)
        weights_XWproduct_w <- c(pmat)
        weights_ranslopes_w <- rvec
        weights_increment_b <- c(
            vapply(icc_X, \(.) if (. == 0) 0 else NA, numeric(1L)),
            vapply(predictors[l2], \(.) .$weight, numeric(1L))
        )

        # Set binary W for level2
        binary_W <- vapply(predictors[l2], \(.) {
            if ("mp_binary" %in% class(.)) .$mean else 0.0
        }, numeric(1L))
    })

    # Return list as env
    return(list2env(l))
}

#' Internal function to solve parameters based on
#' a converted `mp_model` object
#'
#' @noRd
make_parameters <- function(model) {
    # Convert model to old specification
    model |> with(`_convert`) -> env

    # Create list using env
    l <- with(env, {
        # variable counts
        num_X <- length(mean_X)
        num_W <- length(var_W)

        # set binary variable variance
        var_W[binary_W > 0 & binary_W < 1] <- binary_W[binary_W > 0 & binary_W < 1] * (1 - binary_W[binary_W > 0 & binary_W < 1])

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

        # construct the between-cluster correlation and covariance matrix of the level-2 Ws
        cor_WW_b <- matrix(1, nrow = num_W, ncol = num_W)
        cor_WW_b[lower.tri(cor_WW_b)] <- cor_WW_b[upper.tri(cor_WW_b)] <- corr_W(num_W * (num_W - 1) / 2)
        phi_WW_b <- diagonal(sqrt(var_W)) %*% cor_WW_b %*% diagonal(sqrt(var_W))

        # for binary level-2 variables, convert inputted correlations to point-biserial correlations then to covariances
        if (is.null(binary_W)) {
            binary_W_ind <- rep(FALSE, length(var_W))
        } else {
            binary_W_ind <- binary_W > 0 & binary_W < 1
        }
        for (r in seq_len(num_W)) {
            for (i in seq_len(num_W)) {
                if (binary_W_ind[r] == T & i != r) {
                    cor_pointbi <- cor_WW_b[r, i] / sqrt(binary_W[r] * (1 - binary_W[r])) * dnorm(qnorm(binary_W[r]))
                    cor_WW_b[r, i] <- cor_WW_b[i, r] <- cor_pointbi
                    phi_WW_b[r, i] <- phi_WW_b[i, r] <- cor_WW_b[r, i] * sqrt(phi_WW_b[r, r] * phi_WW_b[i, i])
                }
            }
        }

        # construct the between-cluster covariance matrix of the X cluster means and level-2 Ws
        cor_XW_b <- matrix(corr_X_W(num_X * num_W), nrow = num_W, ncol = num_X)
        phi_XW_b <- diagonal(sqrt(var_W)) %*% cor_XW_b %*% diagonal(sqrt(var_X_b))

        # construct the between-cluster covariance matrix
        phi_b <- rbind(cbind(phi_XX_b, t(phi_XW_b)), cbind(phi_XW_b, t(phi_WW_b)))

        # construct the within-cluster covariance matrix
        phi_XWwithXW_w <- phi_XX_w %x% phi_WW_b
        phi_XwithXW_w <- matrix(0, nrow = NROW(phi_XWwithXW_w), ncol = NCOL(phi_XX_w))
        phi_w <- rbind(cbind(phi_XX_w, t(phi_XwithXW_w)), cbind(phi_XwithXW_w, t(phi_XWwithXW_w)))

        # solve for the within-cluster regression coefficients
        weights_scaled <- 1 / sqrt(diagonal(phi_XX_w)) * weights_X_w
        gamma_X_w <- weights_scaled * c(sqrt((var_Y * R2_X_w) / t(weights_scaled) %*% phi_XX_w %*% weights_scaled))

        # solve for the cross-level product coefficients
        weights_scaled <- 1 / sqrt(diagonal(phi_XWwithXW_w)) * weights_XWproduct_w
        if (sum(weights_XWproduct_w) == 0) {
            gamma_XW_w <- weights_XWproduct_w
        } else {
            gamma_XW_w <- weights_scaled * c(sqrt((var_Y * R2_XWproduct_w) / t(weights_scaled) %*% phi_XWwithXW_w %*% weights_scaled))
        }
        gamma_w <- c(gamma_X_w, gamma_XW_w)

        # compute the within-cluster residual variance
        var_e_w <- (1 - icc_Y - R2_X_w - R2_XWproduct_w - R2_ranslopes_w) * var_Y

        # solve for the between-cluster coefficients
        select_weighted <- !is.na(weights_increment_b)
        select_nonweighted <- is.na(weights_increment_b)
        if (sum(select_weighted) != NROW(phi_b)) {
            phi_nonweighted_b <- phi_b[select_nonweighted, select_nonweighted, drop = F]
            phi_weighted_b <- phi_b[select_weighted, select_weighted, drop = F]
            phi_covs_b <- phi_b[select_nonweighted, select_weighted, drop = F]
            resvar_W_b <- phi_weighted_b - t(solve(phi_nonweighted_b) %*% phi_covs_b) %*% phi_covs_b

            # Predefine and only compute for non zero phi_b
            weights_scaled <- numeric(sum(select_weighted))
            sel_weight <- diagonal(resvar_W_b) != 0 # Select out only non-zero diagonals
            weights_scaled[sel_weight] <- 1 / sqrt(diagonal(resvar_W_b)[sel_weight]) * weights_increment_b[select_weighted][sel_weight]
            gamma_weighted_b <- weights_scaled * c(sqrt((var_Y * R2_increment_b) / t(weights_scaled) %*% resvar_W_b %*% weights_scaled))
            gamma_b <- c(gamma_w[seq_len(num_X)], rep(0, num_W))
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
        cor_raneffects <- matrix(1, nrow = num_X + 1, ncol = num_X + 1)
        corvec <- corr_raneffects((num_X + 1) * ((num_X + 1) - 1) / 2)
        cor_raneffects[upper.tri(cor_raneffects)] <- cor_raneffects[lower.tri(cor_raneffects)] <- corvec

        # solve for the random slope variances
        cor_ranslopes <- cor_raneffects[-1, -1, drop = F]
        tau_trace <- (var_Y * R2_ranslopes_w) / sum(diagonal(cor_ranslopes %*% phi_XX_w %*% diagonal(weights_ranslopes_w)))
        var_ranslopes <- if (is.finite(tau_trace)) weights_ranslopes_w * tau_trace else rep(0.0, NROW(weights_ranslopes_w))

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
        tau00 <- 0.5 * (-2 * b + a^2 - 2 * s + 2 * var_Y_b) + 0.5 * sqrt(-4 * b * a^2 + a^4 - 4 * a^2 * s + 4 * a^2 * var_Y_b)

        # compute intercept-slope covariance and construct tau matrix
        tau <- diagonal(sqrt(c(tau00, var_ranslopes))) %*% cor_raneffects %*% diagonal(sqrt(c(tau00, var_ranslopes)))
        cor_raneffects[tau == 0] <- 0

        # compute fixed intercept and construct coefficient matrix
        # the mean of the product from Bohrnstedt & Goldberger Equation 3 simplifies because cov(X_w,W) = 0
        means <- c(model_mean_X, rep(0, num_X + num_W + num_X*num_W)) # Set all W means to 0
        gamma00 <- mean_Y - c(gamma_w, gamma_b) %*% means
        gammas <- c(gamma00, gamma_w, gamma_b)

        # TODO figure out way to check for invalid R2
        # R-square summary
        check_var_Y <- t(gamma_w) %*% phi_w %*% gamma_w + t(gamma_b) %*% phi_b %*% gamma_b + sum(diagonal(tau[-1, -1] %*% phi_XX_w)) + tau00 + var_e_w
        R2check_X_w <- t(gamma_X_w) %*% phi_XX_w %*% gamma_X_w / check_var_Y
        R2check_XW_w <- t(gamma_XW_w) %*% phi_XWwithXW_w %*% gamma_XW_w / check_var_Y
        R2check_ranslopes_w <- sum(diagonal(tau_ranslopes %*% phi_XX_w)) / check_var_Y
        R2check_var_e <- ((1 - icc_Y) * var_Y - t(gamma_w) %*% phi_w %*% gamma_w - sum(diagonal(tau_ranslopes %*% phi_XX_w))) / check_var_Y
        R2check_XW_b <- t(gamma_b) %*% phi_b %*% gamma_b / check_var_Y
        if (sum(select_weighted) != NROW(phi_b)) {
            R2check_increment_b <- t(gamma_weighted_b) %*% resvar_W_b %*% gamma_weighted_b / check_var_Y
        } else {
            R2check_increment_b <- R2check_XW_b
        }
        R2check_totalminusincrement_b <- R2check_XW_b - R2check_increment_b
        R2check_ranicept <- tau00 / check_var_Y

        # variable names
        if (length(names(mean_X)) != 0) {
            vars_Xw <- paste0(names(mean_X), "_w")
            vars_Xb <- paste0(names(mean_X), "_b")
            vars_Xt <- paste0(names(mean_X), "_t")
        } else {
            vars_Xw <- vars_Xb <- vars_Xt <- names(mean_X)
        }
        vars_W <- names(mean_W)
        vars_W[binary_W_ind == T] <- paste0(vars_W[binary_W_ind == T], "_binary")
        vars_XW <- unlist(lapply(vars_Xw, \(x) {
            vapply(vars_W, \(w) {
                paste0(x, "*", w)
            }, character(1L))
        }))

        # collect parameters and construct names
        params_coeff <- matrix(c(gammas), ncol = 1)
        params_res <- matrix(c(var_e_w), ncol = 1)
        row.names(tau) <- colnames(tau) <- c("Icept", vars_Xw)
        row.names(params_coeff) <- c("Icept", vars_Xw, vars_XW, vars_Xb, vars_W)
        row.names(params_res) <- c("Res. Var.")
        colnames(params_coeff) <- colnames(params_res) <- "Value"

        # Return final output
        list(
            mean_Y = mean_Y,
            gammas = params_coeff,
            tau = tau,
            var_e = params_res,
            mean_X = mean_X,
            mean_W = mean_W,
            phi_XX_w = phi_XX_w,
            phi_b = phi_b
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

    # Create temp model with default corrs
    new_model <- clone(model)

    # Find average
    new_model$corrs <- correlations(
        within_cor  = fixed(mean(model$corrs$within_cor)),
        between_cor = fixed(mean(model$corrs$between_cor)),
        randeff_cor = fixed(mean(model$corrs$randeff_cor))
    )

    # output parameters
    new_model |> make_parameters()
}


#' Validate parameters
#' @noRd
is.parameters <- function(x) {
    inherits(x, "mp_parameters")
}

#' Convert `mp_parameters` object to a formula for `lme4`
#'
#' @noRd
to_formula <- function(x, e = globalenv(), nested = FALSE) {

    # Get regression coefficients
    gammas <- x$gammas

    # Obtain number of l1 and l2
    n_l1 <- NROW(x$mean_X)
    n_l2 <- NROW(x$mean_W)

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

    var_l2 <- if (length(names(x$mean_W)) == 0) {
        character()
    } else {
        paste0("cgm(", names(x$mean_W), ")")
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
        raneff_model <- c("1", var_l1[diagonal(x$tau)[-1] != 0])
        raneff_model <- paste0(raneff_model, collapse = " + ")
        model_random <- paste0("(", raneff_model, " | ", "`_id`)")
    }
    lme4_model <- paste0(names(x$mean_Y), " ~ ", model_fixed, " + ", model_random)
    return(as.formula(lme4_model, e))
}

#' Convert `mp_parameters` to list
#' @export
as.list.mp_parameters <- function(x, ...) {
    as.list.environment(x, ...)
}

#' print `mp_parameters` to list
#' @export
print.mp_parameters <- function(x, ...) {
    l <- lapply(as.list(x), function(x) {
        attr(x, '_count') <- NULL
        x
    })
    print(l[names(l) != 'timevar_l1'], ...)
}

#' Internal function to keep online average for `mp_parameters`
#'
#' @export
`+.mp_parameters` <- function(e1, e2) {
    # Check if one is missing
    if (missing(e1)) return(e2)
    if (missing(e2)) return(e1)

    # Don't allow two models
    if (!is.parameters(e1) | !is.parameters(e2)) {
        cli::cli_abort('Can only average two parameters together')
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

