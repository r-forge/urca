cajorls <- function(z, r = 1, reg.number = NULL){
    if ((!inherits(z, "ca.jo")) && (!inherits(z, "cajo.test"))) {
        stop("\nPlease, provide object of class 'ca.jo' or 'cajo.test' as 'z'.\n")
    }
    P <- ncol(z@Z0)
    r <- as.integer(r)
    if((r < 1) || (r > P - 1)){
      stop(paste("Please, provide a cointegration rank 'r' in the interval of 1 to", P - 1, ", in accordance with 'z@V'.\n", sep = " "))
    }
    beta <- matrix(z@V[, 1:r], ncol = r)
    C1 <- diag(r)
    C2 <- matrix(0, nrow = nrow(beta) - r, ncol = r)
    C <- rbind(C1, C2)
    betanorm <- beta %*% solve(t(C) %*% beta)
    ECT <- z@ZK %*% betanorm
    colnames(ECT) <- paste("ect", 1:r, sep = "")
    colnames(betanorm) <- colnames(ECT)
    rownames(betanorm) <- colnames(z@ZK)
    data.mat <- data.frame(z@Z0, ECT, z@Z1)
    text <- colnames(data.mat)[-c(1:P)]
    text1 <- paste(text, "", sep = "+", collapse = "")
    text2 <- paste("~", substr(text1, 1, nchar(text1) - 1))
    if (!is.null(reg.number)) {
        reg.number <- as.integer(reg.number)
        if (reg.number > ncol(z@Z0) || reg.number < 1) {
            stop("\nPlease, provide a valid number of the regression within \n the VECM, numbering from 1 to ", P, ".\n")
        }
        form1 <- formula(paste("z@Z0[, reg.number]", text2, "-1"))
        rlm <- lm(substitute(form1), data = data.mat)
		return(list(rlm = rlm, beta = betanorm))
    }
    else if (is.null(reg.number)) {
        form1 <- formula(paste("z@Z0", text2, "-1"))
        rlm <- lm(substitute(form1), data = data.mat)
		# compute beta statistics
		ZK2 <- z@ZK[, -r] # Y(-1) (2)
		M <- diag(nrow(z@ZK)) - z@Z1 %*% solve(crossprod(z@Z1)) %*% t(z@Z1)
		alpha <- matrix(coef(rlm)[1:r, ], ncol = r, byrow = T) # the coefficients on ecm
		rownames(alpha) <- colnames(coef(rlm)[1:r, ])
		SIGMA <- crossprod(resid(rlm)) / nrow(z@ZK)
		beta.se <- sqrt(diag(kronecker(
		  solve(t(ZK2) %*% M %*% ZK2),
		  solve(t(alpha) %*% solve(SIGMA) %*% alpha)
		)))
		beta.se <- matrix(c(rep(NA, r), beta.se), ncol = r, byrow = T)
		for (i in 1:r) {
		  beta.se[i, ] <- NA
		}
		rownames(beta.se) <- rownames(beta)
		colnames(beta.se) <- colnames(beta)
		beta.t <- betanorm / beta.se
		beta.pval <- 2 * pt(abs(beta.t), df = rlm$df.residual, lower.tail = F)
		if (r == 1) {
		  beta.stats <- data.frame(cbind(betanorm, beta.se, beta.t, beta.pval))
		  colnames(beta.stats) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
		} else {
		  beta.stats <- list()
		  for (ect in 1:ncol(beta)) {
			beta.ect <- data.frame(cbind(round(betanorm[, ect], digits = 15), beta.se[, ect], beta.t[, ect], beta.pval[, ect]))
			colnames(beta.ect) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
			beta.stats[[ect]] <- beta.ect
		  }
		}
		return(list(rlm = rlm, beta = betanorm, beta.stats = beta.stats))
    }
}
