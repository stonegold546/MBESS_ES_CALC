mediation <- 
function (x, mediator, dv, S = NULL, N = NULL, x.location.S = NULL, 
    mediator.location.S = NULL, dv.location.S = NULL, mean.x = NULL, 
    mean.m = NULL, mean.dv = NULL, conf.level = 0.95, bootstrap = FALSE, 
    B = 1000) 
{
    if (bootstrap == TRUE) {
        if (!require(boot)) 
            stop("This function depends on the 'boot' package. Please install the 'boot' package first")
    }
    .mediation <- function(x = x, mediator = mediator, dv = dv, 
        S = S, N = N, x.location.S = x.location.S, mediator.location.S = mediator.location.S, 
        dv.location.S = dv.location.S, mean.x = mean.x, mean.m = mean.m, 
        mean.dv = mean.dv, conf.level = conf.level) {
        if (!is.null(S)) {
            These <- c(x.location.S, mediator.location.S, dv.location.S)
            Cov.Matrix <- as.matrix(S[These, These])
        }
        if (is.null(S)) {
            Data <- na.omit(cbind(x, mediator, dv))
            Cov.Matrix <- var(Data)
            N <- dim(Data)[1]
            mean.dv <- mean(dv)
            mean.x <- mean(x)
            mean.m <- mean(mediator)
            s.dv <- scale(dv)
            s.x <- scale(x)
            s.mediator <- scale(mediator)
            Y.on.X <- lm(dv ~ x)
            resid.Y.on.X <- resid(Y.on.X)
            standardized.Y.on.X <- lm(s.dv ~ s.x)
            standardized.resid.Y.on.X <- resid(standardized.Y.on.X)
            M.on.X <- lm(mediator ~ x)
            resid.M.on.X <- resid(M.on.X)
            standardized.M.on.X <- lm(s.mediator ~ s.x)
            standardized.resid.M.on.X <- resid(standardized.M.on.X)
            Y.on.X.and.M <- lm(dv ~ x + mediator)
            resid.Y.on.X.and.M <- resid(Y.on.X.and.M)
            standardized.Y.on.X.and.M <- lm(s.dv ~ s.x + s.mediator)
            standardized.resid.Y.on.X.and.M <- resid(standardized.Y.on.X.and.M)
            Y.on.M <- lm(dv ~ mediator)
            resid.Y.on.M <- resid(Y.on.M)
            standardized.Y.on.M <- lm(s.dv ~ s.mediator)
            standardized.resid.Y.on.M <- resid(standardized.Y.on.M)
            e.1M <- resid.M.on.X
            e.1Y <- resid.Y.on.X + resid.Y.on.M - resid.Y.on.X.and.M
            standardized.e.1M <- standardized.resid.M.on.X
            standardized.e.1Y <- standardized.resid.Y.on.X + 
                standardized.resid.Y.on.M - standardized.resid.Y.on.X.and.M
            e.0M <- mediator - mean.m
            e.0Y <- dv - mean.dv
            standardized.e.0M <- s.mediator - 0
            standardized.e.0Y <- s.dv - 0
			Residual.Based_Gamma <- as.numeric(1 -(sum(abs(e.1M) + 
                abs(e.1Y)))/(sum(abs(e.0M) + abs(e.0Y))))
                                
			Residual.Based.Standardized_gamma <- as.numeric(1 - 
              (sum(abs(standardized.e.1M) + abs(standardized.e.1Y)))/(sum(abs(standardized.e.0M) + 
              abs(standardized.e.0Y))))
        }
        Cor.Matrix <- cov2cor(Cov.Matrix)
        Dim.Cov.Matrix <- dim(Cov.Matrix)[1]
        s.XY <- Cov.Matrix[Dim.Cov.Matrix, -3]
        S.XX <- Cov.Matrix[1:(Dim.Cov.Matrix - 1), 1:(Dim.Cov.Matrix - 
            1)]
        B.Y_X <- solve(S.XX[1:1, 1:1]) %*% s.XY[1]
        B.Y_X <- cbind(mean.dv - mean.x * B.Y_X, B.Y_X)
        colnames(B.Y_X) <- c("Intercept.Y_X", "c (Regressor)")
        path.c <- B.Y_X[2]
        R2.Y_X <- (t(s.XY[1]) %*% solve(S.XX[1:1, 1:1]) %*% s.XY[1])/Cov.Matrix[Dim.Cov.Matrix, 
            Dim.Cov.Matrix]
        R2.Y_X.Adj <- 1 - ((1 - R2.Y_X) * ((N - 1)/(N - 1 - 1)))
        CI.R2.Y_X <- ci.R2(R2 = R2.Y_X, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 1)
        Model.F.Y_X <- (R2.Y_X/1)/((1 - R2.Y_X)/(N - 1 - 1))
        RMSE.Y_X <- sqrt((1 - R2.Y_X) * ((N - 1) * Cov.Matrix[3, 
            3]/(N - 2)))
        SE.Y_X <- cbind(c(sqrt((1 - R2.Y_X) * ((N - 1) * Cov.Matrix[3, 
            3]/(N - 2))) * sqrt(1/N + mean.x^2/((N - 1) * S.XX[1, 
            1])), sqrt((1 - R2.Y_X)/(N - 1 - 1)) * sqrt(Cov.Matrix[3, 
            3]/S.XX[1, 1])))
        t.Y_X <- t(B.Y_X)/SE.Y_X
        p.Y_X <- 2 * (pt(-1 * abs(t.Y_X), df = N - 1 - 1))
        CL.Low.Y_X <- t(B.Y_X) - qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.Y_X
        CL.Up.Y_X <- t(B.Y_X) + qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.Y_X
        Values.Y_X <- cbind(t(B.Y_X), SE.Y_X, t.Y_X, p.Y_X, CL.Low.Y_X, 
            CL.Up.Y_X)
        colnames(Values.Y_X) <- c("Estimate", "Std. Error", "t value", 
            "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.Y_X <- cbind(RMSE.Y_X, 1, N - 1 - 1, Model.F.Y_X, 
            1 - pf(Model.F.Y_X, 1, N - 1 - 1), R2.Y_X, R2.Y_X.Adj, 
            CI.R2.Y_X$Lower, CI.R2.Y_X$Upper)
        colnames(Model.Fit.Y_X) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.Y_X) <- "Values"
        Regression.of.Y.on.X <- list(Regression.Table = Values.Y_X, 
            Model.Fit = Model.Fit.Y_X)
        B.M_X <- solve(S.XX[1:1, 1:1]) %*% S.XX[2, 1]
        B.M_X <- cbind(mean.m - mean.x * B.M_X, B.M_X)
        colnames(B.M_X) <- c("Intercept.M_X", "a (Regressor)")
        path.a <- B.M_X[2]
        R2.M_X <- (t(S.XX[2, 1]) %*% solve(S.XX[1:1, 1:1]) %*% 
            S.XX[2, 1])/S.XX[2, 2]
        R2.M_X.Adj <- 1 - ((1 - R2.M_X) * ((N - 1)/(N - 1 - 1)))
        CI.R2.M_X <- ci.R2(R2 = R2.M_X, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 1)
        Model.F.M_X <- (R2.M_X/1)/((1 - R2.M_X)/(N - 1 - 1))
        RMSE.M_X <- sqrt((1 - R2.M_X) * ((N - 1) * S.XX[2, 2]/(N - 
            2)))
        SE.M_X <- cbind(c(sqrt((1 - R2.M_X) * ((N - 1) * S.XX[2, 
            2]/(N - 2))) * sqrt(1/N + mean.x^2/((N - 1) * S.XX[1, 
            1])), sqrt((1 - R2.M_X)/(N - 1 - 1)) * sqrt(S.XX[2, 
            2]/S.XX[1, 1])))
        t.M_X <- t(B.M_X)/SE.M_X
        p.M_X <- 2 * (pt(-1 * abs(t.M_X), df = N - 1 - 1))
        CL.Low.M_X <- t(B.M_X) - qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.M_X
        CL.Up.M_X <- t(B.M_X) + qt(1 - (1 - conf.level)/2, df = N - 
            1 - 1) * SE.M_X
        Values.M_X <- cbind(t(B.M_X), SE.M_X, t.M_X, p.M_X, CL.Low.M_X, 
            CL.Up.M_X)
        colnames(Values.M_X) <- c("Estimate", "Std. Error", "t value", 
            "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.M_X <- cbind(RMSE.M_X, 1, N - 1 - 1, Model.F.M_X, 
            1 - pf(Model.F.M_X, 1, N - 1 - 1), R2.M_X, R2.M_X.Adj, 
            CI.R2.M_X$Lower, CI.R2.M_X$Upper)
        colnames(Model.Fit.M_X) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.M_X) <- "Values"
        Regression.of.M.on.X <- list(Regression.Table = Values.M_X, 
            Model.Fit = Model.Fit.M_X)
        B.Y_XM <- solve(S.XX) %*% s.XY
        B.Y_XM <- t(cbind(c(mean.dv - (mean.x * B.Y_XM[1] + mean.m * 
            B.Y_XM[2]), cbind(B.Y_XM))))
        colnames(B.Y_XM) <- c("Intercept.Y_XM", "c.prime (Regressor)", 
            "b (Mediator)")
        path.c.prime <- B.Y_XM[2]
        path.b <- B.Y_XM[3]
        R2.Y_XM <- (t(s.XY) %*% solve(S.XX) %*% s.XY)/Cov.Matrix[Dim.Cov.Matrix, 
            Dim.Cov.Matrix]
        R2.Y_XM.Adj <- 1 - ((1 - R2.Y_XM) * ((N - 1)/(N - 2 - 
            1)))
        CI.R2.Y_XM <- ci.R2(R2 = R2.Y_XM, conf.level = conf.level, 
            Random.Predictors = TRUE, N = N, p = 2)
        Model.F.Y_XM <- (R2.Y_XM/2)/((1 - R2.Y_XM)/(N - 2 - 1))
        RMSE.Y_XM <- sqrt((1 - R2.Y_XM) * ((N - 1) * Cov.Matrix[3, 
            3]/(N - 3)))
        x.prime.x <- cbind(c(N, mean.x * N, mean.m * N), c(mean.x * 
            N, (S.XX[1, 1] * (N - 1) + N * mean.x^2), S.XX[1, 
            2] * (N - 1) + N * mean.x * mean.m), c(mean.m * N, 
            S.XX[1, 2] * (N - 1) + N * mean.x * mean.m, (S.XX[2, 
                2] * (N - 1) + N * mean.m^2)))
        SE.Y_XM <- cbind(sqrt(diag(solve(x.prime.x))) * RMSE.Y_XM)
        t.Y_XM <- t(B.Y_XM)/SE.Y_XM
        p.Y_XM <- 2 * (pt(-1 * abs(t.Y_XM), df = N - 2 - 1))
        CL.Low.Y_XM <- t(B.Y_XM) - qt(1 - (1 - conf.level)/2, 
            df = N - 2 - 1) * SE.Y_XM
        CL.Up.Y_XM <- t(B.Y_XM) + qt(1 - (1 - conf.level)/2, 
            df = N - 2 - 1) * SE.Y_XM
        Values.Y_XM <- cbind(t(B.Y_XM), SE.Y_XM, t.Y_XM, p.Y_XM, 
            CL.Low.Y_XM, CL.Up.Y_XM)
        colnames(Values.Y_XM) <- c("Estimate", "Std. Error", 
            "t value", "p(>|t|)", "Low Conf Limit", "Up Conf Limit")
        Model.Fit.Y_XM <- cbind(RMSE.Y_XM, 2, N - 2 - 1, Model.F.Y_XM, 
            1 - pf(Model.F.Y_XM, 2, N - 2 - 1), R2.Y_XM, R2.Y_XM.Adj, 
            CI.R2.Y_XM$Lower, CI.R2.Y_XM$Upper)
        colnames(Model.Fit.Y_XM) <- c("Residual standard error (RMSE)", 
            "numerator df", "denomenator df", "F-Statistic", 
            "p-value (F)", "R^2", "Adj R^2", "Low Conf Limit", 
            "Up Conf Limit")
        rownames(Model.Fit.Y_XM) <- "Values"
        Regression.of.Y.on.X.and.M <- list(Regression.Table = Values.Y_XM, 
            Model.Fit = Model.Fit.Y_XM)
        s2.X <- Cov.Matrix[1, 1]
        s2.M <- Cov.Matrix[2, 2]
        s2.Y <- Cov.Matrix[3, 3]
        s.YX <- Cov.Matrix[1, 3]
        s.XM <- Cov.Matrix[2, 1]
        s.YM <- Cov.Matrix[3, 2]
        ab <- path.a * path.b
        a.contained <- c((s.YM * s.YX + sqrt(s2.M * s2.Y - s.YM^2) * 
            sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y), (s.YM * 
            s.YX - sqrt(s2.M * s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - 
            s.YX^2))/(s2.X * s2.Y))
        if (path.a > 0) {
            a.contained <- a.contained[a.contained > 0]
            a.contained <- a.contained[abs(a.contained) == max(abs(a.contained))]
        }
        if (path.a < 0) {
            a.contained <- a.contained[a.contained < 0]
            a.contained <- a.contained[abs(a.contained) == max(abs(a.contained))]
        }
        b.contained <- c(sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
            s2.M - s.XM^2), -sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
            s2.M - s.XM^2))
        if (path.b > 0) {
            b.contained <- b.contained[b.contained > 0]
            b.contained <- b.contained[abs(b.contained) == max(abs(b.contained))]
        }
        if (path.b < 0) {
            b.contained <- b.contained[b.contained < 0]
            b.contained <- b.contained[abs(b.contained) == max(abs(b.contained))]
        }
        if (ab > 0) {
            From.a <- a.contained * c(sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
                s2.M - s.XM^2), -sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
                s2.M - s.XM^2))
            From.a <- From.a[From.a > 0]
            From.b <- b.contained * c((s.YM * s.YX + sqrt(s2.M * 
                s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * 
                s2.Y), (s.YM * s.YX - sqrt(s2.M * s2.Y - s.YM^2) * 
                sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y))
            From.b <- From.b[From.b > 0]
            From.b <- From.b[abs(From.b) == max(abs(From.b))]
            From.a <- From.a[abs(From.a) == max(abs(From.a))]
        }
        if (ab < 0) {
            From.a <- a.contained * c(sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
                s2.M - s.XM^2), -sqrt(s2.X * s2.Y - s.YX^2)/sqrt(s2.X * 
                s2.M - s.XM^2))
            From.a <- From.a[From.a < 0]
            From.b <- b.contained * c((s.YM * s.YX + sqrt(s2.M * 
                s2.Y - s.YM^2) * sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * 
                s2.Y), (s.YM * s.YX - sqrt(s2.M * s2.Y - s.YM^2) * 
                sqrt(s2.X * s2.Y - s.YX^2))/(s2.X * s2.Y))
            From.b <- From.b[From.b < 0]
            From.b <- From.b[abs(From.b) == max(abs(From.b))]
            From.a <- From.a[abs(From.a) == max(abs(From.a))]
        }
        Indirect.Effect <- c(Estimate = path.a * path.b)
        Indirect.Effect.Partially.Standardized <- c(Estimate = path.a * 
            path.b/sqrt(Cov.Matrix[3, 3]))
        Index.of.Mediation <- c(Estimate = path.a * path.b * 
            (sqrt(Cov.Matrix[1, 1])/sqrt(Cov.Matrix[3, 3])))
        R2_4.5 <- c(Estimate = (Cor.Matrix[3, 2]^2) - (R2.Y_XM - 
            R2.Y_X))
        R2_4.6 <- c(Estimate = R2.M_X * (R2.Y_XM - R2.Y_X)/(1 - 
            R2.Y_X))
        R2_4.7 <- c(Estimate = (R2.M_X * (R2.Y_XM - R2.Y_X)/(1 - 
            R2.Y_X))/R2.Y_XM)
        Maximum.Possible.Mediation.Effect <- c(Estimate = From.a)
        ab.to.Maximum.Possible.Mediation.Effect_kappa.squared <- c(Estimate = ab/From.a)
        Ratio.of.Indirect.to.Total.Effect <- c(Estimate = 1 - (path.c.prime/path.c))
        Ratio.of.Indirect.to.Direct.Effect <- c(Estimate = path.a * path.b/path.c.prime)
        Success.of.Surrogate.Endpoint <- c(Estimate = path.c/path.a)
        SOS <- c(Estimate = R2_4.5/R2.Y_X)
        if (!is.null(S)) {
            Results.mediation <- list(Y.on.X = Regression.of.Y.on.X, 
                M.on.X = Regression.of.M.on.X, Y.on.X.and.M = Regression.of.Y.on.X.and.M, 
                Indirect.Effect = Indirect.Effect, Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                Index.of.Mediation = Index.of.Mediation, R2_4.5 = R2_4.5, 
                R2_4.6 = R2_4.6, R2_4.7 = R2_4.7, Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                SOS = SOS)
        }
        if (is.null(S)) {
            if (sum(x == 0 | x == 1) != N) {
                Results.mediation <- list(Y.on.X = Regression.of.Y.on.X, 
                  M.on.X = Regression.of.M.on.X, Y.on.X.and.M = Regression.of.Y.on.X.and.M, 
                  Indirect.Effect = Indirect.Effect, Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                  Index.of.Mediation = Index.of.Mediation, R2_4.5 = R2_4.5, 
                  R2_4.6 = R2_4.6, R2_4.7 = R2_4.7, Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                  ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                  Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                  Residual.Based_Gamma = Residual.Based_Gamma, 
                  Residual.Based.Standardized_gamma = Residual.Based.Standardized_gamma, 
                  SOS = SOS)
            }
            if (sum(x == 0 | x == 1) == N) {
                ES <- c(Estimate = ((path.a * path.b)/(sqrt(SE.Y_XM[3]^2 * 
                  path.a^2 + SE.M_X[2]^2 * path.b^2))) * (sqrt(1/sum(x == 
                  0) + 1/sum(x == 1))))
                Results.mediation <- list(Y.on.X = Regression.of.Y.on.X, 
                  M.on.X = Regression.of.M.on.X, Y.on.X.and.M = Regression.of.Y.on.X.and.M, 
                  Indirect.Effect = Indirect.Effect, Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                  Index.of.Mediation = Index.of.Mediation, R2_4.5 = R2_4.5, 
                  R2_4.6 = R2_4.6, R2_4.7 = R2_4.7, Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                  ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                  Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                  Residual.Based_Gamma = Residual.Based_Gamma, 
                  Residual.Based.Standardized_gamma = Residual.Based.Standardized_gamma, 
                  ES.for.two.groups = ES, SOS = SOS)
            }
        }
        return(Results.mediation)
    }
    if (bootstrap == FALSE) {
        return(.mediation(x = x, mediator = mediator, dv = dv, 
            S = S, N = N, x.location.S = x.location.S, mediator.location.S = mediator.location.S, 
            dv.location.S = dv.location.S, mean.x = mean.x, mean.m = mean.m, 
            mean.dv = mean.dv, conf.level = conf.level))
    }
    if (bootstrap == TRUE) {
        if (!require(boot)) 
            stop("This function depends on the package 'boot'. Please install the 'boot' package first.")
        if (!is.null(S)) 
            stop("For the bootstrap procedures to be implemented, you must supply raw data (i.e., not a covariance matrix).")
        Data <- na.omit(cbind(x, mediator, dv))
        N <- dim(Data)[1]
        Boot.This <- function(Data, g) {
            values <- .mediation(x = Data[g, 1], mediator = Data[g, 
                2], dv = Data[g, 3], S = NULL, conf.level = conf.level)
            if (sum(x == 0 | x == 1) != N) {
                Values.to.boot <- c(values$Indirect.Effect, values$Indirect.Effect.Partially.Standardized, 
                  values$Index.of.Mediation, values$R2_4.5, values$R2_4.6, 
                  values$R2_4.7, values$Maximum.Possible.Mediation.Effect, 
                  values$ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  values$Ratio.of.Indirect.to.Total.Effect, values$Ratio.of.Indirect.to.Direct.Effect, 
                  values$Success.of.Surrogate.Endpoint, values$Residual.Based_Gamma, 
                  values$Residual.Based.Standardized_gamma, values$SOS)
            }
            if (sum(x == 0 | x == 1) == N) {
                Values.to.boot <- c(values$Indirect.Effect, values$Indirect.Effect.Partially.Standardized, 
                  values$Index.of.Mediation, values$R2_4.5, values$R2_4.6, 
                  values$R2_4.7, values$Maximum.Possible.Mediation.Effect, 
                  values$ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                  values$Ratio.of.Indirect.to.Total.Effect, values$Ratio.of.Indirect.to.Direct.Effect, 
                  values$Success.of.Surrogate.Endpoint, values$Residual.Based_Gamma, 
                  values$Kris.Standardized, values$ES.for.two.groups, 
                  values$SOS)
            }
            return(Values.to.boot)
        }
        print("Bootstrap resampling has begun. This process may take a considerable amount of time if the number of replications is large, which is optimal for the bootstrap procedure.")
        boot.out <- boot(data = Data, statistic = Boot.This, 
            R = B, stype = "i")
        Mediation.Results <- .mediation(x = x, mediator = mediator, 
            dv = dv, S = NULL, N = N, x.location.S = x.location.S, 
            mediator.location.S = mediator.location.S, dv.location.S = dv.location.S, 
            mean.x = mean.x, mean.m = mean.m, mean.dv = mean.dv, 
            conf.level = conf.level)
        Indirect.Effect <- rbind(c(Mediation.Results$Indirect.Effect, 
            boot.ci(boot.out = boot.out, index = c(1), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(1), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Indirect.Effect) <- c("Estimate", "CI.Low.Perc", 
            "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Indirect.Effect.Partially.Standardized <- rbind(c(Mediation.Results$Indirect.Effect.Partially.Standardized, 
            boot.ci(boot.out = boot.out, index = c(2), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(2), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Indirect.Effect.Partially.Standardized) <- c("Estimate", 
            "CI.Low.Perc", "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Index.of.Mediation <- rbind(c(Mediation.Results$Index.of.Mediation, 
            boot.ci(boot.out = boot.out, index = c(3), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(3), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Index.of.Mediation) <- c("Estimate", "CI.Low.Perc", 
            "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        R2_4.5 <- rbind(c(Mediation.Results$R2_4.5, boot.ci(boot.out = boot.out, 
            index = c(4), conf = conf.level, type = c("perc"))$percent[4:5], 
            boot.ci(boot.out = boot.out, index = c(4), conf = conf.level, 
                type = c("bca"))$bca[4:5]))
        colnames(R2_4.5) <- c("Estimate", "CI.Low.Perc", "CI.Up.Perc", 
            "CI.Low.Bca", "CI.Up.BCa")
        R2_4.6 <- rbind(c(Mediation.Results$R2_4.6, boot.ci(boot.out = boot.out, 
            index = c(5), conf = conf.level, type = c("perc"))$percent[4:5], 
            boot.ci(boot.out = boot.out, index = c(5), conf = conf.level, 
                type = c("bca"))$bca[4:5]))
        colnames(R2_4.6) <- c("Estimate", "CI.Low.Perc", "CI.Up.Perc", 
            "CI.Low.Bca", "CI.Up.BCa")
        R2_4.7 <- rbind(c(Mediation.Results$R2_4.7, boot.ci(boot.out = boot.out, 
            index = c(6), conf = conf.level, type = c("perc"))$percent[4:5], 
            boot.ci(boot.out = boot.out, index = c(6), conf = conf.level, 
                type = c("bca"))$bca[4:5]))
        colnames(R2_4.7) <- c("Estimate", "CI.Low.Perc", "CI.Up.Perc", 
            "CI.Low.Bca", "CI.Up.BCa")
        Maximum.Possible.Mediation.Effect <- rbind(c(Mediation.Results$Maximum.Possible.Mediation.Effect, 
            boot.ci(boot.out = boot.out, index = c(7), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(7), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Maximum.Possible.Mediation.Effect) <- c("Estimate", 
            "CI.Low.Perc", "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        ab.to.Maximum.Possible.Mediation.Effect_kappa.squared <- rbind(c(Mediation.Results$ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
            boot.ci(boot.out = boot.out, index = c(8), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(8), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(ab.to.Maximum.Possible.Mediation.Effect_kappa.squared) <- c("Estimate", 
            "CI.Low.Perc", "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Ratio.of.Indirect.to.Total.Effect <- rbind(c(Mediation.Results$Ratio.of.Indirect.to.Total.Effect, 
            boot.ci(boot.out = boot.out, index = c(9), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(9), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Ratio.of.Indirect.to.Total.Effect) <- c("Estimate", "CI.Low.Perc", 
            "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Ratio.of.Indirect.to.Direct.Effect <- rbind(c(Mediation.Results$Ratio.of.Indirect.to.Direct.Effect, 
            boot.ci(boot.out = boot.out, index = c(10), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(10), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Ratio.of.Indirect.to.Direct.Effect) <- c("Estimate", "CI.Low.Perc", 
            "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Success.of.Surrogate.Endpoint <- rbind(c(Mediation.Results$Success.of.Surrogate.Endpoint, 
            boot.ci(boot.out = boot.out, index = c(11), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(11), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Success.of.Surrogate.Endpoint) <- c("Estimate", 
            "CI.Low.Perc", "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Residual.Based_Gamma <- rbind(c(Mediation.Results$Residual.Based_Gamma, 
            boot.ci(boot.out = boot.out, index = c(12), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(12), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Residual.Based_Gamma) <- c("Estimate", "CI.Low.Perc", 
            "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        Residual.Based.Standardized_gamma <- rbind(c(Mediation.Results$Residual.Based.Standardized_gamma, 
            boot.ci(boot.out = boot.out, index = c(13), conf = conf.level, 
                type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                index = c(13), conf = conf.level, type = c("bca"))$bca[4:5]))
        colnames(Residual.Based.Standardized_gamma) <- c("Estimate", 
            "CI.Low.Perc", "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
        SOS <- rbind(c(Mediation.Results$SOS, boot.ci(boot.out = boot.out, 
            index = c(14), conf = conf.level, type = c("perc"))$percent[4:5], 
            boot.ci(boot.out = boot.out, index = c(14), conf = conf.level, 
                type = c("bca"))$bca[4:5]))
        colnames(SOS) <- c("Estimate", "CI.Low.Perc", "CI.Up.Perc", 
            "CI.Low.Bca", "CI.Up.BCa")
        if (sum(x == 0 | x == 1) != N) {
            Results.mediation <- list(Y.on.X = Mediation.Results$Y.on.X, 
                M.on.X = Mediation.Results$M.on.X, Y.on.X.and.M = Mediation.Results$Y.on.X.and.M, 
                Indirect.Effect = Indirect.Effect, Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                Index.of.Mediation = Index.of.Mediation, R2_4.5 = R2_4.5, 
                R2_4.6 = R2_4.6, R2_4.7 = R2_4.7, Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                Residual.Based_Gamma = Residual.Based_Gamma, 
                Residual.Based.Standardized_gamma = Residual.Based.Standardized_gamma, 
                SOS = SOS)
        }
        if (sum(x == 0 | x == 1) == N) {
            ES.for.two.groups <- rbind(c(Mediation.Results$ES.for.two.groups, 
                boot.ci(boot.out = boot.out, index = c(12), conf = conf.level, 
                  type = c("perc"))$percent[4:5], boot.ci(boot.out = boot.out, 
                  index = c(12), conf = conf.level, type = c("bca"))$bca[4:5]))
            colnames(ES.for.two.groups) <- c("Estimate", "CI.Low.Perc", 
                "CI.Up.Perc", "CI.Low.Bca", "CI.Up.BCa")
            Results.mediation <- list(Y.on.X = Mediation.Results$Y.on.X, 
                M.on.X = Mediation.Results$M.on.X, Y.on.X.and.M = Mediation.Results$Y.on.X.and.M, 
                Indirect.Effect = Indirect.Effect, Indirect.Effect.Partially.Standardized = Indirect.Effect.Partially.Standardized, 
                Index.of.Mediation = Index.of.Mediation, R2_4.5 = R2_4.5, 
                R2_4.6 = R2_4.6, R2_4.7 = R2_4.7, Maximum.Possible.Mediation.Effect = Maximum.Possible.Mediation.Effect, 
                ab.to.Maximum.Possible.Mediation.Effect_kappa.squared = ab.to.Maximum.Possible.Mediation.Effect_kappa.squared, 
                Ratio.of.Indirect.to.Total.Effect = Ratio.of.Indirect.to.Total.Effect, Ratio.of.Indirect.to.Direct.Effect = Ratio.of.Indirect.to.Direct.Effect, 
                Success.of.Surrogate.Endpoint = Success.of.Surrogate.Endpoint, 
                Residual.Based_Gamma = Residual.Based_Gamma, 
                Residual.Based.Standardized_gamma = Residual.Based.Standardized_gamma, 
                ES.for.two.groups = ES.for.two.groups, SOS = SOS)
        }
        return(Results.mediation)
    }
}
