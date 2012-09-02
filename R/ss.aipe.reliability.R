ss.aipe.reliability<- function (model = NULL, type = NULL, width = NULL, S = NULL, conf.level = 0.95, assurance = NULL, data = NULL, i = NULL, cor.est = NULL, lambda = NULL, psi.square = NULL, initial.iter = 500, final.iter = 5000, start.ss = NULL) 
{
    if (sum(search() == "package:sem") != 1) 
        stop("This function depends on the package 'sem', please load 'sem'.")
    if (sum(search() == "package:MASS") != 1) 
        stop("This function depends on the package 'MASS', please load 'MASS'.")
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
    if (sum(model == model.type1) == 1) 
        Model.to.Use <- "Parallel"
    if (sum(model == model.type2) == 1) 
        Model.to.Use <- "True Score"
    if (sum(model == model.type3) == 1) 
        Model.to.Use <- "Congeneric"
    type.1 <- c("Normal Theory", "Normal theory", "normal theory", 
        "nt", "NT")
    type.2 <- c("Factor Analytic", "factor analytic", "Factor analytic", 
        "factor Analytic")
    if (sum(type == type.1, type == type.2) != 1) 
        stop("Assign either Factor Analytic or Normal Theory to 'type'.")
    if (sum(type == type.1) == 1) 
        Type.to.Use <- "Normal Theory"
    if (sum(type == type.2) == 1) 
        Type.to.Use <- "Factor Analytic"
    if (sum(model == model.type1) == 1) {
        if (!is.null(cor.est)) {
            if (!is.null(lambda)) 
                stop("Please enter either cor.est or lambda, but not both.")
            if (is.null(psi.square)) 
                stop("Problem: please enter all of the necessary information: `i', `cor.est' or `lambda', and `psi.square.'")
            if (psi.square <= 0) 
                stop("Problem: `psi.square' must be greater than zero")
            if (length(psi.square) > 1) 
                stop("Problem: 'psi.square' must be one number for the Parallel model. If you want to enter multiple 'psi.square' values, you must use either `True Score' or `Congeneric.'")
            if (is.null(i)) 
                stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
            if (i <= 0) 
                stop("Problem: `i' must be greater than zero")
            if (length(cor.est) >= 2) 
                stop("Problem: you can only enter one 'cor.est' value.")
            lambda.1 <- sqrt(cor.est)
            v.lambda <- matrix(data = lambda.1, nrow = 1, ncol = i)
            v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
            temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                tol.det = 1e-05)
            cov.mat <- Sigma <- temp.mat$Population.Covariance
        }
        if (!is.null(lambda)) {
            if (!is.null(cor.est)) 
                stop("Problem: please enter either cor.est or lambda, but not both.")
            if (is.null(psi.square)) 
                stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
            if (psi.square <= 0) 
                stop("Problem: `psi.square' must be greater than zero.")
            if (length(psi.square) > 1) 
                stop("Problem: 'psi.square' must be one number for the Parallel model. If you want to enter multiple 'psi.square' values, you must use either the True Score or Congeneric model.")
            if (is.null(i)) 
                stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
            if (i <= 0) 
                stop("Problem: `i' must be greater than zero.")
            v.lambda <- matrix(data = lambda, nrow = 1, ncol = i)
            v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
            temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                tol.det = 1e-05)
            cov.mat <- Sigma <- temp.mat$Population.Covariance
        }
        if (!is.null(S)) {
            if (!isSymmetric(S, tol = 1e-05)) 
                stop("Input a symmetric covariance or correlation matrix, 'S'")
            cov.mat <- S
            i <- dim(S)[1]
        }
        if (!is.null(data)) {
            data <- na.omit(data)
            cov.mat <- var(data, y = NULL, na.rm = TRUE)
        }
        sigma.jj <- sum(diag(cov.mat))
        sigma2.Y <- sum(cov.mat)
        p <- ncol(cov.mat)
        alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
        k <- 1 - alpha
        Crit.Value <- qnorm(1 - (1 - conf.level)/2)
        Nec.N <- as.numeric(((Crit.Value^2) * 8 * (k^2) * p)/((width^2) * 
            (p - 1)) + 1)
    }
    if (sum(model == model.type2) == 1) {
        if (!is.null(cor.est)) {
            if (!is.null(lambda)) 
                stop("Problem: please enter either `cor.est' or `lambda', but not both")
            if (is.null(psi.square)) 
                stop("Problem: please enter a vector of `psi.square' values")
            if (is.null(i)) 
                stop("Problem: please enter a value for i")
            if ((i) != length(psi.square)) 
                stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
            if (i <= 0) 
                stop("Problem: `i' must be greater than zero")
            if (length(cor.est) >= 2) 
                stop("You can only enter one 'cor.est' value.")
            lambda.1 <- sqrt(cor.est)
            lambda.vector <- rep(lambda.1, times = i)
            Population.Cov <- covmat.from.cfm(Lambda = lambda.vector, 
                Psi.Square = psi.square)$Population.Covariance
            cor.mat <- cov2cor(Population.Cov)
            cov.mat <- cor.mat
        }
        if (!is.null(lambda)) {
            if (!is.null(cor.est)) 
                stop("Please enter either cor.est or lambda, but not both")
            if (is.null(psi.square)) 
                stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
            if (length(psi.square) != (i)) 
                stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
            if (is.null(i)) 
                stop("You need to enter a quantity for i if you enter lambda and psi.square")
            if (i <= 0) 
                stop("i must be greater than zero")
            if (length(lambda) > 1) 
                stop("'lambda' must be one number for the True Score model. If you want to enter multiple 'lambda' values, you must use the Congeneric model.")
            lambda.vector <- rep(lambda, times = i)
            Population.Cov <- covmat.from.cfm(Lambda = lambda.vector, 
                Psi.Square = psi.square)$Population.Covariance
            cor.mat <- cov2cor(Population.Cov)
            cov.mat <- cor.mat
        }
        if (!is.null(S)) {
            if (!isSymmetric(S, tol = 1e-05)) 
                stop("Input a symmetric covariance or correlation matrix 'S'")
            cor.mat <- cov2cor(S)
            cov.mat <- cor.mat
            i <- dim(S)[1]
        }
        if (!is.null(data)) {
            data <- na.omit(data)
            cor.mat <- cor(data)
            cov.mat <- cor.mat
        }
        p <- ncol(cor.mat)
        j <- cbind(rep(1, times = p))
        Crit.Value <- qnorm(1 - (1 - conf.level)/2)
        step.1 <- (p^2/(p - 1)^2)
        gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
        gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
        gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
        gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
        gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% 
            cor.mat) %*% j)
        gamma.2 <- gamma.2.1 - gamma.2.2
        gamma.final <- gamma.1 * gamma.2
        Nec.N <- as.numeric(((step.1 * gamma.final)/((width/(2 * 
            Crit.Value))^2)) + 1)
    }
    if (sum(model == model.type3) == 1) {
        if (!is.null(cor.est)) {
            if (!is.null(lambda)) 
                stop("Please enter either cor.est or lambda, but not both")
            if (is.null(psi.square)) 
                stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
            if (is.null(i)) 
                stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
            if (i <= 0) 
                stop("i must be greater than zero")
            if ((i) != length(psi.square)) 
                stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
            if (length(cor.est) >= 2) 
                stop("If you have multiple values for 'cor.est', please put them as 'lambda' values. The square root of cor.est equals lambda.")
            print("You entered one value for 'cor.est' with the Congeneric model. This model allows for multiple 'lambda' values. If you have only one 'cor.est' or 'lambda' value, you might want to use the True Score model.")
            lambda.1 <- sqrt(cor.est)
            v.lambda <- matrix(data = lambda.1, nrow = 1, ncol = i)
            v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
            temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                tol.det = 1e-05)
            cov.mat <- Sigma <- temp.mat$Population.Covariance
            cor.mat <- cov2cor(cov.mat)
        }
        if (!is.null(lambda)) {
            if (!is.null(cor.est)) 
                stop("Please enter either cor.est or lambda, but not both")
            if (is.null(psi.square)) 
                stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
            if (is.null(i)) 
                stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
            if (i <= 0) 
                stop("i must be greater than zero")
            if ((i) != length(lambda)) 
                stop("The number of values entered in the lambda vector should be the same as the quantity entered for i")
            if ((i) != length(psi.square)) 
                stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
            v.lambda <- matrix(data = lambda, nrow = 1, ncol = i)
            v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
            temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                tol.det = 1e-05)
            cov.mat <- Sigma <- temp.mat$Population.Covariance
            cor.mat <- cov2cor(cov.mat)
        }
        if (!is.null(data)) {
            data <- na.omit(data)
            cov.mat <- var(data, y = NULL, na.rm = TRUE)
            cor.mat <- cov2cor(cov.mat)
        }
        if (!is.null(S)) {
            if (!isSymmetric(S, tol = 1e-05)) 
                stop("Input a symmetric covariance or correlation matrix 'S'")
            cor.mat <- cov2cor(S)
            cov.mat <- cor.mat
            i <- dim(S)[1]
        }
        p <- ncol(cor.mat)
        j <- cbind(rep(1, times = p))
        Crit.Value <- qnorm(1 - (1 - conf.level)/2)
        step.1 <- (p^2/(p - 1)^2)
        gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
        gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
        gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
        gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
        gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% 
            cor.mat) %*% j)
        gamma.2 <- gamma.2.1 - gamma.2.2
        gamma.final <- gamma.1 * gamma.2
        Nec.N <- as.numeric(((step.1 * gamma.final)/((width/(2 * 
            Crit.Value))^2)) + 1)
    }
    if (!is.null(assurance)) {
        if (assurance > 1) 
            assurance <- assurance/100
        print("An a priori Monte Carlo simulation study has been started so that the exact sample size for the requested condition can be determined. Please be patient, as this process may take several minutes (or longer given the computer and condition).")
        initial.assurance.N <- ceiling(Nec.N) + 1
        if (sum(model == model.type1) == 1) 
            Model.to.Use <- "Parallel"
        if (sum(model == model.type2) == 1) 
            Model.to.Use <- "True Score"
        if (sum(model == model.type3) == 1) 
            Model.to.Use <- "Congeneric"
        n.i <- initial.assurance.N
        if (is.null(start.ss)) {
            Difference <- -1
            while (Difference < 0) {
                CI.Result <- rep(NA, initial.iter)
                for (iters in 1:initial.iter) {
                  sim.data <- var(mvrnorm(n = n.i, mu = rep(0, 
                    i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
                  if (sum(model == model.type1) == 1) 
                    CI.Result.raw <- try(ci.reliability(S = sim.data, 
                      N = n.i, interval.type = "Parallel", analysis.type = Type.to.Use, 
                      conf.level = conf.level), silent = TRUE)
                  if (sum(model == model.type2) == 1) 
                    CI.Result.raw <- try(ci.reliability(S = sim.data, 
                      N = n.i, interval.type = "True Score", analysis.type = Type.to.Use, 
                      conf.level = conf.level), silent = TRUE)
                  if (sum(model == model.type3) == 1) 
                    CI.Result.raw <- try(ci.reliability(S = sim.data, 
                      N = n.i, interval.type = "Congeneric", analysis.type = Type.to.Use, 
                      conf.level = conf.level), silent = TRUE)
                  if (class(CI.Result.raw) == "try-error") {
                    CI.Result[iters] <- NA
                  }
                  if (class(CI.Result.raw) != "try-error") {
                    CI.Result[iters] <- CI.Result.raw$CI.upper - 
                      CI.Result.raw$CI.lower
                  }
                  iters.to.go<-initial.iter-iters
                  max.possible.assurance <- (sum(na.omit(CI.Result)<width) + iters.to.go)/initial.iter
                  if (max.possible.assurance < assurance) {break}
                }
                Difference <- mean(na.omit(CI.Result) < width) - 
                  assurance
                if (Difference < 0) {
                  n.i <- n.i + 1
                }
            }
            initial.n.i <- n.i
        }
        if (!is.null(start.ss)) {
            CI.Result <- NA
            n.i <- start.ss
            Difference <- 0
            for (iters in 1:final.iter) {
                sim.data <- var(mvrnorm(n = n.i, mu = rep(0, 
                  i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
                CI.Result.raw <- try(ci.reliability(S = sim.data, 
                  N = n.i, type = Model.to.Use, analysis.type = Type.to.Use, 
                  conf.level = conf.level), silent = TRUE)
                if (class(CI.Result.raw) == "try-error") {
                  CI.Result[iters] <- NA
                }
                if (class(CI.Result.raw) != "try-error") {
                  CI.Result[iters] <- CI.Result.raw$CI.upper - 
                    CI.Result.raw$CI.lower
                }                                               
            }
            Difference <- mean(na.omit(CI.Result) < width) - 
                assurance
            if (Difference > 0) 
                while (Difference > 0) {
                  CI.Result <- rep(NA, final.iter)
                  for (iters in 1:final.iter) {
                    sim.data <- var(mvrnorm(n = n.i, mu = rep(0, 
                      i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
                    CI.Result.raw <- try(ci.reliability(S = sim.data, 
                      N = n.i, type = Model.to.Use, analysis.type = Type.to.Use, 
                      conf.level = conf.level), silent = TRUE)
                    if (class(CI.Result.raw) == "try-error") {
                      CI.Result[iters] <- NA
                    }
                    if (class(CI.Result.raw) != "try-error") {
                      CI.Result[iters] <- CI.Result.raw$CI.upper - 
                        CI.Result.raw$CI.lower
                    }
                    iters.to.go<-final.iter-iters
                    min.possible.assurance <- (sum(na.omit(CI.Result)<width))/final.iter 
                    if (min.possible.assurance > assurance) {break}
                  }
                  Difference <- mean(na.omit(CI.Result) < width) - 
                    assurance
                  if (Difference > 0) {
                    n.i <- n.i - 1
                  }
                }
        }
        Difference <- -1
        while (Difference < 0) {
            CI.Result <- rep(NA, final.iter)
            for (iters in 1:final.iter) {
                sim.data <- var(mvrnorm(n = n.i, mu = rep(0, 
                  i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
                CI.Result.raw <- try(ci.reliability(S = sim.data, 
                  N = n.i, type = Model.to.Use, analysis.type = Type.to.Use, 
                  conf.level = conf.level), silent = TRUE)
                if (class(CI.Result.raw) == "try-error") {
                  CI.Result[iters] <- NA
                }
                if (class(CI.Result.raw) != "try-error") {
                  CI.Result[iters] <- CI.Result.raw$CI.upper - 
                    CI.Result.raw$CI.lower
                }
              
                iters.to.go<-final.iter-iters
                max.possible.assurance <- (sum(na.omit(CI.Result)<width) + iters.to.go)/final.iter 
                if (max.possible.assurance < assurance) {break}

            }
            Difference <- mean(na.omit(CI.Result) < width) - 
                assurance
            if (Difference < 0) {
                n.i <- n.i + 1
            }
        }
        Nec.N.assurance <- n.i
        empirical.assurance <- round(mean(na.omit(CI.Result) < 
            width), 3)
        print(paste("A sample size of", n.i, "leads to an empirical assurance of", 
            round(mean(na.omit(CI.Result) < width), 3)))
    }
    if (is.null(assurance)) 
        return(list(Required.Sample.Size = ceiling(Nec.N)))
    if (!is.null(assurance)) 
        return(list(Required.Sample.Size = ceiling(Nec.N.assurance), 
            width = width, specified.assurance = assurance, empirical.assurance = empirical.assurance, 
            final.iter = final.iter))
}
