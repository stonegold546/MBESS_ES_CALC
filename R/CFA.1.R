CFA.1 <- function(S, N, equal.loading = FALSE, equal.error = FALSE, package = "sem", se = "standard", ...) {
    if (!isSymmetric(S, tol = 1e-05)) 
        stop("Input a symmetric covariance or correlation matrix 'S'")
    q <- nrow(S)
    if (package == "sem") {
        if (!require(sem)) 
            stop("This function depends on the 'sem' package. Please install the 'sem' package first")
        
        try(detach("package:lavaan"), silent = TRUE)
        x <- matrix(NA, nrow = q, ncol = 1)
        x <- paste("x", row(x), sep = "")
        
        if (equal.loading) {
            lamda <- matrix(rep("lamda", q), nrow = q, ncol = 1)
        } else {
            lamda <- matrix(NA, nrow = q, ncol = 1)
            lamda <- paste("lamda", row(lamda), sep = "")
        }
        
        if (equal.error) {
            psi.sq <- matrix(rep("psi.sq", q), nrow = q, ncol = 1)
        } else {
            psi.sq <- matrix(NA, nrow = q, ncol = 1)
            psi.sq <- paste("psi.sq", row(psi.sq), sep = "")
        }
        
        model.1 <- cbind(paste("ksi", "->", x), lamda)
        model.2 <- cbind(paste(x, "<->", x), psi.sq)
        model <- rbind(model.1, model.2)
        start <- matrix(NA, nrow = nrow(model), ncol = 1)
        model <- cbind(model, start)
        model <- rbind(model, c(paste("ksi", "<->", "ksi"), NA, 1))
        class(model) <- "mod"
        
        rownames(S) <- x
        colnames(S) <- x
        model.fitted <- sem(model, S, N)
        converged <- model.fitted$convergence
        if (converged) {
            if (equal.loading) 
                k <- 1 else k <- q
            Factor.Loadings <- model.fitted$coeff[1:k]
            if (equal.error) 
                m <- 1 else m <- q
            Indicator.var <- model.fitted$coeff[k + 1:m]
            Parameter.cov <- model.fitted$vcov
        } else {
            Factor.Loadings <- NA
            Indicator.var <- NA
            Parameter.cov <- NULL
        }
        result <- list(Model = model, Factor.Loadings = Factor.Loadings, Indicator.var = Indicator.var, 
            Parameter.cov = Parameter.cov, converged = converged)
        return(result)
    } else if (package == "lavaan") {
        if (!require(lavaan)) 
            stop("This function depends on the 'lavaan' package. Please install the 'lavaan' package first")
        
        try(detach("package:sem"), silent = TRUE)
        
        
        colnames(S) <- rownames(S) <- paste("y", 1:q, sep = "")
        
        if (equal.loading) {
            loadingName <- rep("a1", q)
        } else {
            loadingName <- paste("a", 1:q, sep = "")
        }
        
        if (equal.error) {
            errorName <- rep("b1", q)
        } else {
            errorName <- paste("b", 1:q, sep = "")
        }
        model <- "f1 =~ NA*y1 + "
        loadingLine <- paste(paste(loadingName, "*", colnames(S), sep = ""), collapse = " + ")
        factorLine <- "f1 ~~ 1*f1\n"
        errorLine <- paste(paste(colnames(S), " ~~ ", errorName, "*", colnames(S), 
            sep = ""), collapse = "\n")
        model <- paste(model, loadingLine, "\n", factorLine, errorLine, "\n")
        try(fit <- lavaan(model, sample.cov = S, sample.nobs = N, se = se, ...), 
            silent = TRUE)
        converged <- fit@Fit@converged
        loading <- unique(as.vector(fit@Model@GLIST$lambda))
        error <- unique(diag(fit@Model@GLIST$theta))
        if (!all(error > 0)) 
            converged <- FALSE
        if (converged) {
            if (se == "none") {
                paramCov <- NULL
            } else {
                paramCov <- vcov(fit)
            }
        } else {
            loading <- NA
            error <- NA
            if (se == "none") {
                paramCov <- NULL
            } else {
                paramCov <- NA
            }
        }
        result <- list(Model = model, Factor.Loadings = loading, Indicator.var = error, 
            Parameter.cov = paramCov, converged = converged)
        return(result)
        
    } else {
        stop("The package should be either 'sem' or 'lavaan' package.")
    }
} 
