
#### ci.reliability
#### Updated 4/24/10

ci.reliability <- function (S = NULL, data = NULL, N = NULL, model = "True-Score Equivalent", 
    type = "Factor Analytic", conf.level = 0.95, interval = TRUE, Bootstrap = FALSE, B = 10000, BootstrapCI = "BCa") 
{
	if(!require(sem)) 
		stop("This function depends on the 'sem' package. Please install the 'sem' package \n    as you installed the 'MBESS' package")
    model.type1 <- c("Parallel", "SB", "Spearman Brown", "Spearman-Brown", 
        "sb", "parallel")
    model.type2 <- c("True Score", "True Score Equivalent", "True-Score Equivalent", 
        "Equivalent", "Tau Equivalent", "Chronbach", "Cronbach", 
        "cronbach", "alpha", "true score", "tau-equivalent", 
        "Tau-Equivalent", "True-Score", "true-score")
    model.type3 <- c("Congeneric", "congeneric", "omega", "Omega")
    if (sum(model == model.type1, model == model.type2, model == 
        model.type3) != 1) 
        stop("Assign one and only one of the three types of models to 'model'.")
    alpha <- 1 - conf.level
	if(Bootstrap == FALSE) {
		if (type == "Factor Analytic") {
			if (!is.null(data)) {
				S <- var(na.omit(data))
				N <- dim(data)[1]
			}
			if (!isSymmetric(S, tol = 1e-05)) 
				stop("Input a symmetric covariance or correlation matrix 'S.'")
			if (is.null(N)) 
				stop("Since only 'S' is entered, 'N' is also needed.")
			q <- nrow(S)
			if (sum(model == model.type1) == 1) 
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
			if (sum(model == model.type2) == 1) 
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
			if (sum(model == model.type3) == 1) 
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
			l <- length(result$Factor.Loadings)
			p <- length(result$Indicator.var)
			k <- nrow(result$Parameter.cov)
			if (l == 1) {
				u <- result$Factor.Loadings * q
				names(u) <- NULL
				var.u <- result$Parameter.cov[1, 1] * q * q
			}
			else {
				u <- sum(result$Factor.Loadings)
				var.u <- sum(result$Parameter.cov[1:q, 1:q])
			}
			if (p == 1) {
				v <- result$Indicator.var * q
				names(v) <- NULL
				var.v <- result$Parameter.cov[k, k] * q * q
			}
			else {
				v <- sum(result$Indicator.var)
				var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
			}
			if (l != 1 & p != 1) 
				cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
			if (l == 1 & p != 1) 
				cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
			if (l == 1 & p == 1) 
				cov.u.v <- result$Parameter.cov[2, 1] * q * q
			omega <- u^2/(u^2 + v)
			if(interval == "T" | interval == "t") {
				D1 <- 2 * u * v/(u^2 + v)^2
				D2 <- (-1) * u^2/(u^2 + v)^2
				se.omega <- sqrt(D1^2 * var.u + D2^2 * var.v + 2 * D1 * D2 * cov.u.v)
				z <- qnorm(1 - alpha/2)
				ci.lower <- omega - z * se.omega
				if (ci.lower < 0) {
					ci.lower = 0
				}
				ci.upper <- omega + z * se.omega
				if (ci.upper > 1) {
            		ci.upper = 1
				}
			} else {
				ci.lower <- ci.upper <- se.omega <- -1
			}
			ci <- list(CI.lower = ci.lower, CI.upper = ci.upper, 
            Estimated.reliability = omega, SE.reliability = se.omega, 
            Conf.Level = conf.level)
			return(ci)
		}
		if (type == "Normal Theory") {
			if (sum(model == model.type3) == 1) 
				stop("The Congeneric model can not be used with `Normal Theory.'")
			if (!is.null(data)) { 
				data <- na.omit(data)
				if (dim(data)[1] == dim(data)[2]) {
					if (sum(round(data, 4) == t(round(data, 4))) == (dim(data)[1] * dim(data)[2])) 
					cov.data <- data
				}
				if (dim(data)[1] != dim(data)[2]) {
					cov.data <- var(data, y = NULL, na.rm = TRUE)
					N <- dim(data)[1]
				}
			}
			if (!is.null(S)) 
				cov.data <- S
			if (is.null(N)) 
				stop("Since only 'S' is entered, 'N' is also needed.")
			if (sum(model == model.type1) == 1) {
				sigma.jj <- sum(diag(cov.data))
				sigma2.Y <- sum(cov.data)
				p <- ncol(cov.data)
				alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
				if(interval) {
					Crit.Value <- qnorm(1 - (1 - conf.level)/2)
					variance <- (2 * (1 - alpha)^2 * p)/((N - 1) * (p - 1))
					SE <- sqrt((2 * (1 - alpha)^2 * p)/((N - 1) * (p - 1)))
					UB <- alpha + (Crit.Value * SE)
					if (UB > 1) {
						UB = 1
					}
					LB <- alpha - (Crit.Value * SE)
					if (LB < 0) {
						LB = 0
					}
				} else {
					LB <- UB <- SE <- -1
				}
				CI <- list(CI.lower = LB, CI.upper = UB, Estimated.reliability = alpha, 
					SE.reliability = SE, Conf.Level = conf.level)
			}
			if (sum(model == model.type2) == 1) {
				sigma.jj <- sum(diag(cov.data))
				sigma2.Y <- sum(cov.data)
				p <- ncol(cov.data)
				alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
				if(interval == "T" | interval == "t") {
					cor.mat <- cov2cor(cov.data)
					p <- ncol(cor.mat)
					j <- cbind(rep(1, times = p))
					step.1 <- (p^2/(p - 1)^2)
					gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
					gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
					gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
					gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
					gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% cor.mat) %*% j)
					gamma.2 <- gamma.2.1 - gamma.2.2
					gamma.final <- gamma.1 * gamma.2
					n <- N - 1
					variance <- (step.1 * gamma.final)/n
					SE <- sqrt(variance)
					Crit.Value <- qnorm(1 - (1 - conf.level)/2)
					UB <- alpha + (Crit.Value * SE)
					if (UB > 1) {
						UB = 1
					}
					LB <- alpha - (Crit.Value * SE)
					if (LB < 0) {
						LB = 0
					}
				} else {
					LB <- UB <- SE <- -1
				}
				CI <- list(CI.lower = LB, CI.upper = UB, Estimated.reliability = alpha, 
					SE.reliability = SE, Conf.Level = conf.level)
			}
			return(CI)
		}
	} else {
		if(!require(boot))
			stop("This function depends on the 'boot' package. Please install the 'boot' package \n    as you installed the 'MBESS' package")
		BootstrapCI.type1 <- c("perc", "Perc", "percentile", "Percentile", "percentile ci", "Percentile ci", "Percentile CI", "percentile CI")
		BootstrapCI.type2 <- c("bca", "Bca", "BCa", "BCA", "Bias-corrected and acceleration", "Bias-Corrected and Acceleration", 
		"bias-corrected and acceleration") 
		if (sum(BootstrapCI == BootstrapCI.type1, BootstrapCI == BootstrapCI.type2) != 1) 
			stop("Assign one and only one of the two types of percentile methodmodels to 'model'.")
		if (is.null(data))
			stop("Data is required for running bootstrap CI")
		if (type == "Factor Analytic") {
			model <- 0
			S <- var(data, y = NULL, na.rm = TRUE)
			N <- dim(data)[1]
			q <- nrow(S)
			if (sum(model == model.type1) == 1) {
				model <- 1
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
			} else if(sum(model == model.type2) == 1)  {
				model <- 2
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
			} else {
				model <- 3
				result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
			}
			l <- length(result$Factor.Loadings)
			p <- length(result$Indicator.var)
			k <- nrow(result$Parameter.cov)
			if (l == 1) {
				u <- result$Factor.Loadings * q
				names(u) <- NULL
				var.u <- result$Parameter.cov[1, 1] * q * q
			} else {
				u <- sum(result$Factor.Loadings)
				var.u <- sum(result$Parameter.cov[1:q, 1:q])
			}
			if (p == 1) {
				v <- result$Indicator.var * q
				names(v) <- NULL
				var.v <- result$Parameter.cov[k, k] * q * q
			} else {
				v <- sum(result$Indicator.var)
				var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
			}
			if (l != 1 & p != 1) 
				cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
			if (l == 1 & p != 1) 
				cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
			if (l == 1 & p == 1) 
				cov.u.v <- result$Parameter.cov[2, 1] * q * q
			omega <- u^2/(u^2 + v)
			.bs.omega <- function(data, model, i) {
				S <- var(data[i, ], y = NULL, na.rm = TRUE)
				N <- dim(data)[1]
				q <- nrow(S)
				if (model == 1) {
					result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
				} else if (model == 2) {
					result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
				} else {
					result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
				}
				l <- length(result$Factor.Loadings)
				p <- length(result$Indicator.var)
				k <- nrow(result$Parameter.cov)
				if (l == 1) {
					u <- result$Factor.Loadings * q
					names(u) <- NULL
					var.u <- result$Parameter.cov[1, 1] * q * q
				} else {
					u <- sum(result$Factor.Loadings)
					var.u <- sum(result$Parameter.cov[1:q, 1:q])
				}
				if (p == 1) {
					v <- result$Indicator.var * q
					names(v) <- NULL
					var.v <- result$Parameter.cov[k, k] * q * q
				} else {
					v <- sum(result$Indicator.var)
					var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
				}
				if (l != 1 & p != 1) 
					cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
				if (l == 1 & p != 1) 
					cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
				if (l == 1 & p == 1) 
					cov.u.v <- result$Parameter.cov[2, 1] * q * q
				omega <- u^2/(u^2 + v)    
				return(omega)
			}
			boot.out <- boot(data = data, statistic = .bs.omega, R = B, stype = "i", model = model)
			se.omega <- sd(boot.out$t)
			CI.output <- NULL
			CI.lower <- NULL
			CI.upper <- NULL
			if (sum(BootstrapCI == BootstrapCI.type1) == 1) {
				CI.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
				CI.lower <- CI.output[4]
				CI.upper <- CI.output[5]
			} else { 
			CI.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
			CI.lower <- CI.output[4]
			CI.upper <- CI.output[5]
			}
			ci <- list(CI.lower = CI.lower, CI.upper = CI.upper, 
				Estimated.reliability = omega, SE.reliability = se.omega, 
				Conf.Level = conf.level)	 	
			return(ci)
		}
		if (type == "Normal Theory") {
			if (sum(model == model.type3) == 1) 
				stop("The Congeneric model can not be used with `Normal Theory.'")
			cov.data <- var(data, y = NULL, na.rm = TRUE)
			N <- dim(data)[1]
			sigma.jj <- sum(diag(cov.data))
			sigma2.Y <- sum(cov.data)
			p <- ncol(cov.data)
			alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
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
			se.alpha <- sd(boot.out$t)
			CI.output <- NULL
			CI.lower <- NULL
			CI.upper <- NULL
			if (sum(BootstrapCI == BootstrapCI.type1) == 1) {
				CI.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
				CI.lower <- CI.output[4]
				CI.upper <- CI.output[5]
			} else { 
				CI.output <- boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
				CI.lower <- CI.output[4]
				CI.upper <- CI.output[5]
			}
			ci <- list(CI.lower = CI.lower, CI.upper = CI.upper, 
            Estimated.reliability = alpha, SE.reliability = se.alpha, 
            Conf.Level = conf.level)	 	
			return(ci)
		}
	}
}
