ci.reliability <- function(data = NULL, S = NULL, N = NULL, type = "omega", analysis.type = "default", 
    interval.type = "normal-theory", B = 1000, conf.level = 0.95) {
    input <- match.call()
    if (!(is.vector(B) && (length(B) == 1) && is.numeric(B))) 
        stop("Please put a number in the number of bootstrap argument!")
    if (!(is.vector(conf.level) && (length(conf.level) == 1) && is.numeric(conf.level))) 
        stop("Please put a number in the confidence level!")
    alpha <- 1 - conf.level
    if (!is.null(data)) {
        S <- var(na.omit(data))
        N <- dim(data)[1]
    }
    if (!isSymmetric(S, tol = 1e-05)) 
        stop("Input a symmetric covariance or correlation matrix 'S.'")
    if (is.null(N)) 
        stop("Since only 'S' is entered, 'N' is also needed.")
    N <- as.numeric(N)
    q <- nrow(S)
    type1 <- c("True Score", "True Score Equivalent", "True-Score Equivalent", "Equivalent", 
        "Tau Equivalent", "Chronbach", "Cronbach", "cronbach", "alpha", "true score", 
        "tau-equivalent", "Tau-Equivalent", "True-Score", "true-score")
    type2 <- c("Congeneric", "congeneric", "omega", "Omega")
    analysis.type1 <- c("analytic", "theoretical", "formula")
    analysis.type2 <- c("factor", "cfa", "factor analytic", "factor analysis")
    interval.type0 <- c("none", "na")
    interval.type1 <- c("parallel", "sb")
    interval.type2 <- c("normal-theory", "wald", "normal")
    interval.type3 <- c("adf", "wls")
    interval.type4 <- c("fisher", "naivefisher")
    interval.type5 <- c("bonett")
    interval.type6 <- c("ll", "likelihood")
    interval.type7 <- c("perc", "percentile", "percentile ci")
    interval.type8 <- c("bca", "boot", "bootstrap", "bias-corrected and acceleration")
    interval.type9 <- c("feldt", "feldt65", "feldt1965", "f", "fdist")
    interval.type10 <- c("hakstianwhalen", "hw", "hakstian", "hw76", "hw1976", "cuberoot")
    interval.type11 <- c("hakstianbarchard", "hb", "hb00", "hb2000", "randomitem")
    interval.type12 <- c("siotani", "shf", "shf85", "siotani85", "siotani1985", "f2")
    interval.type13 <- c("intraclass correlation", "icc", "modifiedfisher", "improvedfisher")
    interval.type14 <- c("bsi", "standardboot", "sdboot")
    
    type <- tolower(type)
    if (!type %in% c(type1, type2)) 
        stop("Please provide a correct type of reliability: 'alpha' or 'omega'.")
    analysis.type <- tolower(analysis.type)
    if (analysis.type %in% "default") {
        ifelse(type %in% type1, analysis.type <- analysis.type1[1], analysis.type <- analysis.type2[1])
    }
    if (!analysis.type %in% c(analysis.type1, analysis.type2)) 
        stop("Please provide how to analyze a reliability coefficient: 'analytic' or 'cfa'.")
    
    interval.type <- tolower(interval.type)
    if (!interval.type %in% c(interval.type0, interval.type1, interval.type2, interval.type3, 
        interval.type4, interval.type5, interval.type6, interval.type7, interval.type8, 
        interval.type9, interval.type10, interval.type11, interval.type12, interval.type13)) 
        stop("Please provide how to analyze the confidence interval of reliability. Specify as 'none' if do not need one.")
    
    boot.out <- NULL
    
    # Find the point estimate
    relia <- NULL
    if (analysis.type %in% analysis.type2) {
        result <- NULL
        if (type %in% type1) {
            if (interval.type %in% interval.type1) {
                relia <- .ciCFA(S, N, type = 1, interval = FALSE)$relia
            } else {
                relia <- .ciCFA(S, N, type = 2, interval = FALSE)$relia
            }
        } else {
            relia <- .ciCFA(S, N, type = 3, interval = FALSE)$relia
        }
    } else if (analysis.type %in% analysis.type1) {
        if (type %in% type2) 
            stop("Coefficient omega cannot be analyzed analytically.")
        sigma.jj <- sum(diag(S))
        sigma2.Y <- sum(S)
        relia <- (q/(q - 1)) * (1 - sigma.jj/sigma2.Y)
    }
    # model = 'True-Score Equivalent', type = 'Factor Analytic', interval = TRUE,
    # Bootstrap = FALSE, BootstrapCI = 'BCa'
    
    # model.type1 <- c('Parallel', 'SB', 'Spearman Brown', 'Spearman-Brown', 'sb',
    # 'parallel') model.type2 <- c('True Score', 'True Score Equivalent',
    # 'True-Score Equivalent', 'Equivalent', 'Tau Equivalent', 'Chronbach',
    # 'Cronbach', 'cronbach', 'alpha', 'true score', 'tau-equivalent',
    # 'Tau-Equivalent', 'True-Score', 'true-score') model.type3 <- c('Congeneric',
    # 'congeneric', 'omega', 'Omega') if (sum(model == model.type1, model ==
    # model.type2, model == model.type3) != 1) stop('Assign one and only one of the
    # three types of models to 'model'.')
    
    
    # Find Interval Estimate
    ci.lower <- NULL
    ci.upper <- NULL
    se <- NULL
    if (interval.type %in% interval.type0) {
        ci.lower <- ci.upper <- se <- -1
    } else if (interval.type %in% interval.type1) {
        if (type %in% type2) 
            stop("Paralel model-based confidence interval cannot be applied to coefficient omega.")
        if (analysis.type %in% analysis.type1) {
            Crit.Value <- qnorm(1 - (1 - conf.level)/2)
            variance <- (2 * (1 - relia)^2 * q)/((N - 1) * (q - 1))
            se <- as.vector(sqrt((2 * (1 - relia)^2 * q)/((N - 1) * (q - 1))))
            ci.upper <- relia + (Crit.Value * se)
            if (ci.upper > 1) 
                ci.upper = 1
            ci.lower <- relia - (Crit.Value * se)
            if (ci.lower < 0) 
                ci.lower = 0
        } else if (analysis.type %in% analysis.type2) {
            result <- .ciCFA(S, N, conf.level = conf.level, type = 1, interval = TRUE)
            ci.lower <- result$ci.lower
            ci.upper <- result$ci.upper
            se <- result$se
        }
    } else if (interval.type %in% interval.type2) {
        if (type %in% type1) {
            if (analysis.type %in% analysis.type1) {
                cor.mat <- cov2cor(S)
                j <- cbind(rep(1, times = q))
                step.1 <- (q^2/(q - 1)^2)
                gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
                gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
                gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
                gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
                gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% cor.mat) %*% 
                  j)
                gamma.2 <- gamma.2.1 - gamma.2.2
                gamma.final <- gamma.1 * gamma.2
                variance <- (step.1 * gamma.final)/(N - 1)
                se <- as.vector(sqrt(variance))
                Crit.Value <- qnorm(1 - (1 - conf.level)/2)
                ci.upper <- relia + (Crit.Value * se)
                if (ci.upper > 1) 
                  ci.upper = 1
                ci.lower <- relia - (Crit.Value * se)
                if (ci.lower < 0) 
                  ci.lower = 0
            } else if (analysis.type %in% analysis.type2) {
                result <- .ciCFA(S, N, conf.level = conf.level, type = 2, interval = TRUE)
                ci.lower <- result$ci.lower
                ci.upper <- result$ci.upper
                se <- result$se
            }
        } else if (type %in% type2) {
            if (analysis.type %in% analysis.type1) {
                stop("Coefficient omega cannot be analyzed analytically.")
            } else if (analysis.type %in% analysis.type2) {
                result <- .ciCFA(S, N, conf.level = conf.level, type = 3, interval = TRUE)
                ci.lower <- result$ci.lower
                ci.upper <- result$ci.upper
                se <- result$se
            }
        }
    } else if (interval.type %in% interval.type3) {
        if (is.null(data)) 
            stop("Data is required for running ADF confidence interval")
        if (type %in% type1) {
            se <- .seReliabilityAdf(data)
        } else {
            result <- .wlsOmegaInterval(data)
            relia <- result$relia
            se <- result$se
        }
        Crit.Value <- qnorm(1 - (1 - conf.level)/2)
        ci.upper <- relia + (Crit.Value * se)
        if (ci.upper > 1) 
            ci.upper <- 1
        ci.lower <- relia - (Crit.Value * se)
        if (ci.lower < 0) 
            ci.lower <- 0
    } else if (interval.type %in% interval.type4) {
        result <- .fisherCIrelia(relia, N, alpha)
        se <- result$se
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
    } else if (interval.type %in% interval.type5) {
        result <- .iccCIrelia(relia, N, q, alpha, bonett = TRUE)
        se <- result$se
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
    } else if (interval.type %in% interval.type6) {
        result <- .likelihoodCIrelia(S, N)
        if (type %in% type1) {
            ci.lower <- result[1, 1]
            ci.upper <- result[1, 2]
        } else {
            ci.lower <- result[2, 1]
            ci.upper <- result[2, 2]
        }
    } else if (interval.type %in% interval.type9) {
        result <- .fMethod(relia, N - 1, (N - 1) * (q - 1), 1 - conf.level)
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
    } else if (interval.type %in% c(interval.type10, interval.type11)) {
        Correction <- 1
        if (interval.type %in% interval.type11) {
            averageCov <- mean(S[lower.tri(S, diag = FALSE)])
            Snew <- matrix(averageCov, q, q)
            diag(Snew) <- mean(diag(S))
            lstar <- det(S)/det(Snew)
            l <- (-N * log(lstar))/((q^2 + q - 4)/2)
            cstar <- 1.452 - (0.464 * l) + (0.046 * l^2)
            Correction <- min(cstar, 1)
        }
        dfn <- (N - 1) * Correction
        dfd <- (N - 1) * (q - 1) * Correction
        a <- (1 - (2/(9 * dfn)))/(1 - (2/(9 * dfd)))
        sigma.num <- (2 * q/(9 * dfd)) * ((1 - relia)^(2/3))
        sigma.denom <- (1 - (2/(9 * dfn)))^2
        se <- sqrt(sigma.num/sigma.denom)
        Crit.Value <- qnorm(1 - (1 - conf.level)/2)
        transform.relia <- (1 - relia)^(1/3)
        transform.lower <- transform.relia + (Crit.Value * se)
        transform.upper <- transform.relia - (Crit.Value * se)
        ci.lower <- 1 - ((a^3) * (transform.lower^3))
        ci.upper <- 1 - ((a^3) * (transform.upper^3))
        try(if (ci.lower < 0) 
            ci.lower = 0)
        try(if (ci.upper > 1) 
            ci.upper = 1)
    } else if (interval.type %in% interval.type12) {
        result <- .fMethod(relia, N, N * (q - 1), 1 - conf.level)
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
    } else if (interval.type %in% interval.type13) {
        result <- .iccCIrelia(relia, N, q, alpha, bonett = FALSE)
        se <- result$se
        ci.lower <- result$ci.lower
        ci.upper <- result$ci.upper
    } else {
        if (!require(boot)) 
            stop("This function depends on the 'boot' package. Please install the 'boot' package \n    as you installed the 'MBESS' package")
        if (is.null(data)) 
            stop("Data is required for running bootstrap CI")
        boot.out <- NULL
        if (analysis.type %in% analysis.type2) {
            reliaType <- 3
            if (type %in% type1) {
                ifelse(interval.type %in% interval.type1, reliaType <- 1, reliaType <- 2)
            }
            .bs.omega <- function(data, reliaType, i) {
                S <- var(data[i, ], y = NULL, na.rm = TRUE)
                N <- dim(data)[1]
                omega <- .ciCFA(S, N, type = reliaType, interval = FALSE)$relia
                return(omega)
            }
            boot.out <- boot(data = data, statistic = .bs.omega, R = B, stype = "i", 
                reliaType = reliaType)
        } else if (analysis.type %in% analysis.type1) {
            .bs.alpha <- function(data, i) {
                cov.data <- var(data[i, ], y = NULL, na.rm = TRUE)
                N <- dim(data)[1]
                sigma.jj <- sum(diag(cov.data))
                sigma2.Y <- sum(cov.data)
                p <- ncol(cov.data)
                alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
                return(alpha)
            }
            boot.out <- boot(data = data, statistic = .bs.alpha, R = B, stype = "i")
        }
        ci.output <- NULL
        if (interval.type %in% interval.type7) {
            ci.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
            ci.lower <- ci.output[4]
            ci.upper <- ci.output[5]
        } else if (interval.type %in% interval.type8) {
            if (B < nrow(data)) 
                warnings("The number of bootstrap samples is less than the number of cases.\nIf 'bca' cannot be calculated, please make sure that the number of cases is greater than\nthe number of bootstrap samples.")
            # From https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html
            ci.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
            ci.lower <- ci.output[4]
            ci.upper <- ci.output[5]
        } else if (interval.type %in% interval.type14) {
            se <- apply(boot.out$t, 2, sd, na.rm = TRUE)
            Crit.Value <- qnorm(1 - (1 - conf.level)/2)
            ci.upper <- relia + (Crit.Value * se)
            if (ci.upper > 1) 
                ci.upper <- 1
            ci.lower <- relia - (Crit.Value * se)
            if (ci.lower < 0) 
                ci.lower <- 0
        }
    }
    if (!is.numeric(ci.lower)) 
        ci.lower <- NA
    if (!is.numeric(ci.upper)) 
        ci.upper <- NA
    out <- list(est = relia, se = se, ci.lower = ci.lower, ci.upper = ci.upper, conf.level = conf.level, 
        est.type = type, analysis.type = analysis.type, interval.type = interval.type, 
        call = input, boot = boot.out)
    return(out)
}

.ciCFA <- function(S, N, conf.level = 0.95, type = 3, interval = TRUE, ...) {
    usedSE <- "standard"
    if (!interval) 
        usedSE <- "none"
    q <- nrow(S)
    if (type == 1) 
        result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE, 
            package = "lavaan", se = usedSE, ...)
    if (type == 2) 
        result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE, 
            package = "lavaan", se = usedSE, ...)
    if (type == 3) 
        result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE, 
            package = "lavaan", se = usedSE, ...)
    omega <- NULL
    se.omega <- NULL
    ci.lower <- NULL
    ci.upper <- NULL
    if (result$converged) {
        l <- length(result$Factor.Loadings)
        p <- length(result$Indicator.var)
        k <- nrow(result$Parameter.cov)
        if (l == 1) {
            u <- result$Factor.Loadings * q
            names(u) <- NULL
        } else {
            u <- sum(result$Factor.Loadings)
        }
        if (p == 1) {
            v <- result$Indicator.var * q
            names(v) <- NULL
        } else {
            v <- sum(result$Indicator.var)
        }
        omega <- u^2/(u^2 + v)
    }
    if (!is.null(omega) && (omega > 1 | omega < 0)) 
        result$converged <- FALSE
    if (result$converged) {
        if (interval) {
            if (l == 1) {
                var.u <- result$Parameter.cov[1, 1] * q * q
            } else {
                var.u <- sum(result$Parameter.cov[1:q, 1:q])
            }
            if (p == 1) {
                var.v <- result$Parameter.cov[k, k] * q * q
            } else {
                var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
            }
            if (l != 1 & p != 1) 
                cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
            if (l == 1 & p != 1) 
                cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
            if (l == 1 & p == 1) 
                cov.u.v <- result$Parameter.cov[2, 1] * q * q
            D1 <- 2 * u * v/(u^2 + v)^2
            D2 <- (-1) * u^2/(u^2 + v)^2
            se.omega <- sqrt(D1^2 * var.u + D2^2 * var.v + 2 * D1 * D2 * cov.u.v)
            z <- qnorm(1 - (1 - conf.level)/2)
            ci.lower <- omega - z * se.omega
            if (ci.lower < 0) 
                ci.lower = 0
            ci.upper <- omega + z * se.omega
            if (ci.upper > 1) 
                ci.upper = 1
        }
        return(list(relia = omega, se = se.omega, ci.lower = ci.lower, ci.upper = ci.upper))
    } else {
        return(list(relia = NA, se = NA, ci.lower = NA, ci.upper = NA))
    }
}

.wlsOmegaInterval <- function(data) {
    library(lavaan)
    p <- ncol(data)
    colnames(data) <- paste("y", 1:p, sep = "")
    loadingName <- paste("a", 1:p, sep = "")
    errorName <- paste("b", 1:p, sep = "")
    model <- "f1 =~ NA*y1 + "
    loadingLine <- paste(paste(loadingName, "*", colnames(data), sep = ""), collapse = " + ")
    factorLine <- "f1 ~~ 1*f1\n"
    errorLine <- paste(paste(colnames(data), " ~~ ", errorName, "*", colnames(data), 
        sep = ""), collapse = "\n")
    sumLoading <- paste("loading :=", paste(loadingName, collapse = " + "), "\n")
    sumError <- paste("error :=", paste(errorName, collapse = " + "), "\n")
    relia <- "relia := (loading^2) / ((loading^2) + error) \n"
    model <- paste(model, loadingLine, "\n", factorLine, errorLine, "\n", sumLoading, 
        sumError, relia)
    fit <- lavaan(model, data = data, estimator = "WLS")
    allEst <- fit@Fit@est
    est <- allEst[length(allEst)]
    allSE <- fit@Fit@se
    se <- allSE[length(allSE)]
    return(list(relia = est, se = se))
}

.seReliabilityAdf <- function(data) {
    data <- na.omit(data)
    N <- dim(data)[1]
    ni <- dim(data)[2]
    S <- cov(data)
    vecs.S <- NULL
    for (i in 1:ni) {
        vecs.S <- c(vecs.S, S[i:ni, i])
    }
    item.mean <- apply(data, 2, mean)
    off.diag.1 <- ((2 * ni)/(ni - 1))
    off.diag.2 <- (sum(diag(S))/(sum(S)^2))
    off.diag <- off.diag.1 * off.diag.2
    delta <- matrix(rep(off.diag, (ni * ni)), nrow = ni)
    on.diag.1 <- ((-ni)/(ni - 1))
    on.diag.2 <- ((sum(S) - sum(diag(S)))/(sum(S)^2))
    on.diag <- on.diag.1 * on.diag.2
    diag(delta) <- on.diag
    est.delta <- NULL
    for (i in 1:ni) {
        est.delta <- c(est.delta, delta[i:ni, i])
    }
    add <- 0
    for (i in 1:N) {
        dev.y <- scale(data, center = TRUE, scale = FALSE)
        temp <- dev.y %*% t(dev.y)
        si <- NULL
        for (j in 1:ni) {
            si <- c(si, temp[j:ni, j])
        }
        d.si <- si - vecs.S
        product <- t(est.delta) %*% d.si %*% t(d.si) %*% est.delta
        add <- add + product
    }
    result <- add/(N * (N - 1))
    return(as.vector(sqrt(result)))
}

.fisherCIrelia <- function(relia, n, alpha = 0.05) {
    Z <- 0.5 * log((1 + relia)/(1 - relia))
    SE <- sqrt(1/(n - 3))
    Crit <- qnorm(1 - (alpha/2))
    MarginError <- SE * Crit
    ZLower <- Z - MarginError
    ZUpper <- Z + MarginError
    Lower <- (exp(2 * ZLower) - 1)/(exp(2 * ZLower) + 1)
    Upper <- (exp(2 * ZUpper) - 1)/(exp(2 * ZUpper) + 1)
    try(if (Lower < 0) 
        Lower = 0)
    try(if (Upper > 1) 
        Upper = 1)
    return(list(z = Z, se = SE, ci.lower = Lower, ci.upper = Upper))
}

.iccCIrelia <- function(relia, n, k, alpha = 0.05, bonett = TRUE) {
    PointEstimate <- log(1 - relia)
    Crit <- qnorm(1 - (alpha/2))
    SE <- NULL
    if (bonett) {
        SE <- sqrt(2 * k/((k - 1) * (n - 2)))
    } else {
        SE <- sqrt(2 * k/((k - 1) * n))
    }
    MarginError <- SE * Crit
    Lower <- 1 - exp(PointEstimate + MarginError)
    Upper <- 1 - exp(PointEstimate - MarginError)
    try(if (Lower < 0) 
        Lower = 0)
    try(if (Upper > 1) 
        Upper = 1)
    return(list(z = PointEstimate, se = SE, ci.lower = Lower, ci.upper = Upper))
}

.fMethod <- function(relia, df1, df2, alpha = 0.05) {
    fa <- qf(alpha/2, df1, df2)
    fb <- qf(1 - alpha/2, df1, df2)
    ci.lower <- 1 - ((1 - relia) * fb)
    ci.upper <- 1 - ((1 - relia) * fa)
    return(list(ci.lower = ci.lower, ci.upper = ci.upper))
}


.likelihoodCIrelia <- function(CM, N) {
    if (!require(OpenMx)) {
        try(source("http://openmx.psyc.virginia.edu/getOpenMx.R"), silent = TRUE)
        tryCatch(library(OpenMx), error = function(e) {
            stop("The OpenMx program cannot be loaded. Please search 'how to install OpenMx' in the Internet.")
        })
    }
    p <- dim(CM)[1]
    vecs.CM <- NULL
    for (i in 1:p) {
        vecs.CM <- c(vecs.CM, CM[1:p, i])
    }
    varnames <- colnames(CM)
    if (is.null(varnames)) {
        for (i in 1:p) {
            temp <- paste("x", i, sep = "")
            varnames <- c(varnames, temp)
        }
    }
    # finding Alpha
    A <- NULL
    U <- NULL
    L <- NULL
    matrixC <- mxMatrix(type = "Symm", nrow = p, ncol = p, free = TRUE, values = vecs.CM, 
        name = "C")
    matrixP <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, values = p, 
        name = "p")
    ALPHA <- mxAlgebra(expression = (p/(p - 1)) * (1 - (tr(C)/sum(C))), name = "Alpha")
    ALPHAObj <- mxMLObjective(covariance = "C", dimnames = varnames)
    Data <- mxData(observed = CM, type = "cov", numObs = N)
    ALPHAci <- mxCI("Alpha")
    AlphaModel <- mxModel("FindAlpha", matrixC, matrixP, ALPHA, ALPHAObj, Data, ALPHAci)
    AlphaFit <- mxRun(AlphaModel, intervals = T)
    # finding Omega
    matrixA <- mxMatrix(type = "Full", nrow = p, ncol = 1, free = TRUE, values = 0.5, 
        name = "A")
    matrixL <- mxMatrix(type = "Symm", nrow = 1, ncol = 1, free = FALSE, values = 1, 
        name = "L")
    matrixU <- mxMatrix(type = "Diag", nrow = p, ncol = p, free = TRUE, values = 0.2, 
        lbound = 1e-04, name = "U")
    A <- matrixA
    L <- matrixL
    U <- matrixU    
    OMEGA <- mxAlgebra(expression = (sum(A) * sum(A))/((sum(A) * sum(A)) + sum(U)), 
        name = "Omega")
    algebraR <- mxAlgebra(expression = A %*% L %*% t(A) + U, name = "R")
    OMEGAObj <- mxMLObjective(covariance = "R", dimnames = varnames)
    OMEGAci <- mxCI("Omega")
    OmegaModel <- mxModel("FindOmega", matrixA, matrixL, matrixU, OMEGA, algebraR, 
        OMEGAObj, Data, OMEGAci)
    OmegaFit <- mxRun(OmegaModel, intervals = T)
    # Return the matrix of result, which CI of alpha is in the top row and CI of
    # omega is in the bottom row.
    result <- matrix(rep(NA, 4), ncol = 2)
    rownames(result) <- c("Alpha", "Omega")
    colnames(result) <- c("Lower", "Upper")
    try(result[1, 1:2] <- AlphaFit@output$confidenceIntervals[1:2])
    try(result[2, 1:2] <- OmegaFit@output$confidenceIntervals[1:2])
    return(result)
} 
